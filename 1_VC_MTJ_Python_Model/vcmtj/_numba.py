"""Optional Numba acceleration for stateful pulse-sequence evolution."""
import numpy as np
import math
from . import _numpy

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised in installations without Numba
    NUMBA_AVAILABLE = False
    njit = None

_PARAMETER_NAMES = (
    "vratio_p", "vratio_ap", "vfit_min", "vfit_max", "t_valid_min_ns",
    "y0_a", "y0_b", "y0_c", "t50_a", "t50_b", "t50_c",
    "amp_a", "amp_b", "amp_c", "w_a", "w_b", "w_c",
    "st50_a", "st50_b", "st50_c", "samp_a", "samp_b", "samp_c",
    "sw", "eta0", "seta", "amp_min", "w_min_ns", "psw_min", "psw_max",
)

def _evolve_python(voltage, width, initial_state, z, uniforms, c):
    (ratio_p, ratio_ap, vmin, vmax, tmin,
     y_a, y_b, y_c, t_a, t_b, t_c, a_a, a_b, a_c, w_a, w_b, w_c,
     st_a, st_b, st_c, sa_a, sa_b, sa_c, sw, eta0, seta,
     amin, wmin, pmin, pmax) = c
    n_pulses, n_devices = uniforms.shape
    probabilities = np.empty((n_pulses, n_devices))
    switched = np.empty((n_pulses, n_devices), dtype=np.bool_)
    states = np.empty((n_pulses + 1, n_devices), dtype=np.int8)
    states[0] = initial_state
    # Fixed device eta factors mirror VA initial_step; O(N) storage, never
    # pulse-by-device directional probability matrices.
    eta_scales = np.empty((2, n_devices))
    for device in range(n_devices):
        eta = eta0 + seta * z[device, 5]
        eta_scales[0, device] = math.exp(min(-0.5 * eta, 709.0))
        eta_scales[1, device] = math.exp(min(0.5 * eta, 709.0))
    for pulse in range(n_pulses):
        for device in range(n_devices):
            state = states[pulse, device]
            vf = voltage[pulse, device] / (ratio_ap if state == 1 else ratio_p)
            probability = 0.0
            if vf >= vmin and width[pulse, device] >= tmin:
                vf = min(vf, vmax)
                vv = vf * vf
                y0 = y_a * vv + y_b * vf + y_c
                t50 = t_a * vv + t_b * vf + t_c
                t50 += (st_a * vv + st_b * vf + st_c) * z[device, 2]
                amp = a_a * vv + a_b * vf + a_c
                amp = max(amp + (sa_a * vv + sa_b * vf + sa_c) * z[device, 3], amin)
                spread = max(w_a * vv + w_b * vf + w_c + sw * z[device, 4], wmin)
                x = 2.0 * math.log(4.0) / spread * (width[pulse, device] - t50)
                term = math.exp(-abs(x))
                sigmoid = 1.0 / (1.0 + term) if x >= 0.0 else term / (1.0 + term)
                common = min(max(y0 + amp * sigmoid, pmin), pmax)
                scale = eta_scales[state, device]
                probability = min(max(-math.expm1(math.log1p(-common) * scale), pmin), pmax)
            did_switch = uniforms[pulse, device] < probability
            probabilities[pulse, device] = probability
            switched[pulse, device] = did_switch
            if did_switch:
                state = 1 - state
            states[pulse + 1, device] = state
    return probabilities, switched, states

if NUMBA_AVAILABLE:
    _evolve_numba = njit(cache=True)(_evolve_python)
else:
    _evolve_numba = _evolve_python

def simulate_sequence(write_voltage_v, pulse_width_ns, initial_state, z, uniforms, p):
    """Calculate current-direction Psw and evolve states inside the JIT loop.

    Uniforms come from the same NumPy generator used by the fallback, keeping
    switching_seed semantics independent of compilation and device variation.
    """
    if not NUMBA_AVAILABLE:
        return _numpy.simulate_sequence(write_voltage_v, pulse_width_ns, initial_state, z, uniforms, p)
    coefficients = np.array([getattr(p, name) for name in _PARAMETER_NAMES])
    return _evolve_numba(write_voltage_v, pulse_width_ns, initial_state, z, uniforms, coefficients)
