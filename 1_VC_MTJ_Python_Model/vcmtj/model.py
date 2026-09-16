"""Public API for write switching and resistance analysis."""
from dataclasses import dataclass
import numpy as np
from .parameters import ModelParameters
from . import _numba, _numpy

P = 0
AP = 1

def _finite(value, name, dtype=np.float64):
    array = np.asarray(value, dtype=dtype)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array

def _validate_states(value, name="mtj_state"):
    array = _finite(value, name)
    if not np.all((array == P) | (array == AP)):
        raise ValueError(f"{name} must contain only P (0) or AP (1)")
    return array.astype(np.int8)

def _normalize_analysis(state, population, *values):
    # Only state establishes a device axis. Scalar state leaves all inputs
    # as sweep coordinates, even when a sweep length equals n_devices.
    if population is None:
        arrays = np.broadcast_arrays(*values, state)
        return arrays[:-1], arrays[-1], np.zeros(6)
    if not isinstance(population, DevicePopulation):
        raise TypeError("population must be a DevicePopulation")
    n_devices = len(population)
    if state.ndim and state.shape[-1] == n_devices:
        # Device-aware inputs follow ordinary NumPy broadcasting. A sweep
        # against per-device states must explicitly use (..., 1).
        arrays = np.broadcast_arrays(*values, state)
    else:
        arrays = [value[..., None] for value in np.broadcast_arrays(*values, state)]
        arrays = np.broadcast_arrays(*arrays, np.empty(n_devices))[:-1]
    return arrays[:-1], arrays[-1], population._matrix()

@dataclass(frozen=True)
class DevicePopulation:
    """Represent fixed device-to-device Monte-Carlo variation.

    Each device stores six independent, dimensionless standard-normal
    variables: z_rp, z_tmr, z_t50, z_a, z_w, and z_eta. Every array has shape
    (n_devices,) and is read-only. Reusing a population across pulses or
    simulations therefore represents the same physical devices; switching
    randomness is controlled separately by switching_seed in MTJModel.simulate.
    """
    z_rp: np.ndarray
    z_tmr: np.ndarray
    z_t50: np.ndarray
    z_a: np.ndarray
    z_w: np.ndarray
    z_eta: np.ndarray

    def __post_init__(self):
        arrays = []
        for name in ("z_rp", "z_tmr", "z_t50", "z_a", "z_w", "z_eta"):
            value = _finite(getattr(self, name), name).copy()
            if value.ndim != 1 or value.size < 1:
                raise ValueError(f"{name} must be a nonempty one-dimensional array")
            value.setflags(write=False)
            arrays.append(value)
        if len({value.size for value in arrays}) != 1:
            raise ValueError("all device variation arrays must have equal length")
        for name, value in zip(("z_rp", "z_tmr", "z_t50", "z_a", "z_w", "z_eta"), arrays):
            object.__setattr__(self, name, value)

    def __len__(self):
        """Return the number of devices in the population."""
        return self.z_rp.size

    @classmethod
    def monte_carlo(cls, n_devices, device_seed=None):
        """Create a reproducible Monte-Carlo device population.

        Parameters
        ----------
        n_devices : int
            Number of physical device realizations. Must be positive.
        device_seed : int or None, optional
            Seed for the six independent standard-normal device variables.
            It affects device-to-device variation only. The default None
            requests a non-deterministic seed from NumPy.

        Returns
        -------
        DevicePopulation
            Fixed read-only variation arrays, each with shape (n_devices,).
        """
        if isinstance(n_devices, bool) or not isinstance(n_devices, (int, np.integer)) or n_devices < 1:
            raise ValueError("n_devices must be a positive integer")
        rng = np.random.default_rng(device_seed)
        values = rng.standard_normal((6, int(n_devices)))
        return cls(*values)

    @classmethod
    def nominal(cls, n_devices=1):
        """Create one or more identical nominal devices with all z values zero."""
        if isinstance(n_devices, bool) or not isinstance(n_devices, (int, np.integer)) or n_devices < 1:
            raise ValueError("n_devices must be a positive integer")
        zeros = [np.zeros(int(n_devices)) for _ in range(6)]
        return cls(*zeros)

    def _matrix(self):
        return np.column_stack((self.z_rp, self.z_tmr, self.z_t50, self.z_a, self.z_w, self.z_eta))

@dataclass(frozen=True)
class SimulationResult:
    """Store the complete result of a write simulation.

    Attributes
    ----------
    probabilities : numpy.ndarray
        Switching probabilities with shape (n_pulses, n_devices).
    switched : numpy.ndarray
        Boolean switching outcomes with shape (n_pulses, n_devices).
    states : numpy.ndarray
        P/AP trajectory with shape (n_pulses + 1, n_devices). Row zero is the
        initial state; each subsequent row follows one write pulse.
    resistance_ohm : numpy.ndarray or None
        Resistance aligned one-to-one with states, in ohms. It is None when no
        read voltage was requested.
    tmr : numpy.ndarray or None
        Dimensionless device TMR with shape (n_devices,) for a scalar read
        voltage, or None when no read voltage was requested.

    Accessing these stored attributes never reruns the simulation.
    """
    probabilities: np.ndarray
    switched: np.ndarray
    states: np.ndarray
    resistance_ohm: np.ndarray | None
    tmr: np.ndarray | None

@dataclass(frozen=True)
class ResistanceResult:
    """Store resistance and TMR evaluated at an actual MTJ read voltage.

    Attributes
    ----------
    resistance_ohm : numpy.ndarray
        State-dependent resistance in ohms.
    tmr : numpy.ndarray
        Dimensionless tunnel magnetoresistance.

    Shapes follow the documented sweep/device broadcasting rules. A voltage
    sweep of length M with N Monte-Carlo devices returns (M, N) arrays; scalar
    nominal inputs return zero-dimensional arrays.
    """
    resistance_ohm: np.ndarray
    tmr: np.ndarray

class MTJModel:
    """Evaluate the calibrated VC-MTJ compact model.

    The model provides stochastic write simulation, probability-only analysis,
    and static resistance/TMR analysis. With no explicit ModelParameters, it
    uses the calibrated defaults corresponding to the reference Verilog-A
    implementation. NumPy handles vectorized analysis, while stateful
    simulation automatically uses optional Numba acceleration when available.
    """

    def __init__(self, parameters=None):
        """Create a model using calibrated defaults or supplied parameters."""
        self.parameters = ModelParameters() if parameters is None else parameters
        if not isinstance(self.parameters, ModelParameters):
            raise TypeError("parameters must be a ModelParameters instance")

    def evaluate_psw(self, write_voltage_v, pulse_width_ns, initial_state=AP, population=None):
        """Evaluate switching probability without performing switching.

        Parameters
        ----------
        write_voltage_v : float or array-like
            Actual voltage across the MTJ terminals during writing, in volts.
            It is mapped internally to Vfit using the current P/AP state. With
            scalar state, an array is a voltage sweep.
        pulse_width_ns : float or array-like
            Effective write-pulse width in nanoseconds. Arrays broadcast with
            voltage and state. Widths below t_valid_min_ns give zero switching
            probability.
        initial_state : int or array-like, optional
            State at which probability is evaluated. Use P or AP; the default
            is AP. With an N-device population, shape (..., N) establishes a
            per-device state axis.
        population : DevicePopulation or None, optional
            Fixed Monte-Carlo devices. If omitted, nominal parameters are
            evaluated. Device variation is not resampled by this method.

        Returns
        -------
        numpy.ndarray
            Switching probability. Scalar nominal inputs return a
            zero-dimensional array. A sweep of length M with N devices returns
            shape (M, N); per-device state (..., N) retains its broadcast shape.

        Notes
        -----
        This auxiliary method performs no Bernoulli trial, state update, or
        resistance evaluation and is not a prerequisite for simulate().
        """
        voltage = _finite(write_voltage_v, "write_voltage_v")
        width = _finite(pulse_width_ns, "pulse_width_ns")
        state = _validate_states(initial_state, "initial_state")
        if np.any(width < 0.0):
            raise ValueError("pulse_width_ns must be nonnegative")
        (voltage, width), state, z = _normalize_analysis(state, population, voltage, width)
        return _numpy.evaluate_psw(voltage, width, state, z, self.parameters)

    def resistance(self, read_voltage_v, mtj_state=AP, population=None):
        """Evaluate static resistance and TMR without a write simulation.

        Parameters
        ----------
        read_voltage_v : float or array-like
            Actual voltage across the MTJ terminals during reading, in volts.
            Rp and TMR use this voltage directly without Vfit mapping. With
            scalar state, an array is an R-V sweep.
        mtj_state : int or array-like, optional
            State for resistance selection. Use P or AP; the default is AP.
            With N devices, shape (..., N) represents per-device states,
            including a trajectory returned by simulate().
        population : DevicePopulation or None, optional
            Fixed Monte-Carlo devices. If omitted, nominal Rp/TMR parameters
            are used.

        Returns
        -------
        ResistanceResult
            Resistance in ohms and dimensionless TMR. A voltage sweep of length
            M with N devices returns shape (M, N). For scalar voltage and a
            state trajectory (K, N), resistance is (K, N) and TMR is (N,).
        """
        voltage = _finite(read_voltage_v, "read_voltage_v")
        state = _validate_states(mtj_state)
        (v,), state, z = _normalize_analysis(state, population, voltage)
        resistance, tmr = _numpy.resistance(v, state, z, self.parameters)
        if population is not None and voltage.ndim == 0:
            _, tmr = _numpy.resistance(voltage, P, z, self.parameters)
        return ResistanceResult(resistance, tmr)

    def simulate(self, write_voltage_v, pulse_width_ns, initial_state=AP, population=None,
                 switching_seed=None, read_voltage_v=None):
        """Simulate stochastic VC-MTJ switching.

        Parameters
        ----------
        write_voltage_v : float or array-like
            Actual voltage across the MTJ terminals during writing, in volts.
            A scalar represents one pulse; a one-dimensional array represents
            a shared sequence. A two-dimensional array may specify
            (n_pulses, n_devices) values.
        pulse_width_ns : float or array-like
            Effective write-pulse width in nanoseconds. Its pulse/device shape
            broadcasts with write_voltage_v.
        initial_state : int or array-like, optional
            State before the first pulse. Use P or AP. The default is AP; an
            array may provide one state per device.
        population : DevicePopulation or None, optional
            Fixed Monte-Carlo devices. If omitted, one nominal device is used.
            Device variation remains fixed throughout the pulse sequence.
        switching_seed : int or None, optional
            Seed controlling Bernoulli switching trials only, independently of
            device_seed. The default None requests a non-deterministic seed.
        read_voltage_v : float or None, optional
            Actual MTJ terminal voltage for post-switching resistance/TMR
            evaluation, in volts. The default None skips read evaluation. The
            current simulation interface accepts a scalar read voltage.

        Returns
        -------
        SimulationResult
            Probabilities and switching outcomes have shape
            (n_pulses, n_devices); states have shape
            (n_pulses + 1, n_devices). With a read voltage, resistance matches
            the state shape and TMR has shape (n_devices,).
        """
        pop = DevicePopulation.nominal() if population is None else population
        if not isinstance(pop, DevicePopulation):
            raise TypeError("population must be a DevicePopulation")
        n_devices = len(pop)
        voltage = self._pulse_matrix(write_voltage_v, "write_voltage_v", n_devices)
        width = self._pulse_matrix(pulse_width_ns, "pulse_width_ns", n_devices)
        shape = np.broadcast_shapes(voltage.shape, width.shape, (1, n_devices))
        voltage, width = np.broadcast_to(voltage, shape).copy(), np.broadcast_to(width, shape).copy()
        if np.any(width < 0.0):
            raise ValueError("pulse_width_ns must be nonnegative")
        initial = np.broadcast_to(_validate_states(initial_state, "initial_state"), (n_devices,)).copy()
        uniforms = np.random.default_rng(switching_seed).random(shape)
        z = pop._matrix()
        engine = _numba.simulate_sequence if _numba.NUMBA_AVAILABLE else _numpy.simulate_sequence
        probabilities, switched, states = engine(voltage, width, initial, z, uniforms, self.parameters)
        resistance_ohm = tmr = None
        if read_voltage_v is not None:
            read = _finite(read_voltage_v, "read_voltage_v")
            if read.ndim != 0:
                raise ValueError("simulate currently requires scalar read_voltage_v")
            read_result = self.resistance(read, states, pop)
            resistance_ohm, tmr = read_result.resistance_ohm, read_result.tmr
        return SimulationResult(probabilities, switched, states, resistance_ohm, tmr)

    @staticmethod
    def _pulse_matrix(value, name, n_devices):
        array = _finite(value, name)
        if array.ndim == 0:
            return array.reshape(1, 1)
        if array.ndim == 1:
            return array[:, None]
        if array.ndim == 2 and array.shape[1] in (1, n_devices):
            return array
        raise ValueError(f"{name} must be scalar, a sequence, or (n_pulses, n_devices)")
