"""Numerical behavior, reproducibility, and VA-parity tests."""
from dataclasses import replace
import math
import numpy as np
import pytest

from vcmtj import AP, P, DevicePopulation, MTJModel, ModelParameters
from vcmtj import _numba, _numpy

def _literal_va_psw(voltage, width, state, z, p):
    ratio = p.vratio_ap if state == AP else p.vratio_p
    vf = voltage / ratio
    if vf < p.vfit_min or width < p.t_valid_min_ns:
        return 0.0
    vf = min(vf, p.vfit_max)
    q = lambda stem: getattr(p, stem + "_a") * vf**2 + getattr(p, stem + "_b") * vf + getattr(p, stem + "_c")
    y0 = q("y0")
    t50 = q("t50") + q("st50") * z[2]
    amp = max(q("amp") + q("samp") * z[3], p.amp_min)
    spread = max(q("w") + p.sw * z[4], p.w_min_ns)
    common = min(max(y0 + amp / (1.0 + math.exp(-2.0 * math.log(4.0) / spread * (width - t50))), p.psw_min), p.psw_max)
    eta_scale = math.exp((0.5 if state == AP else -0.5) * (p.eta0 + p.seta * z[5]))
    return min(max(1.0 - math.exp(math.log1p(-common) * eta_scale), p.psw_min), p.psw_max)

def test_state_constants():
    assert P == 0 and AP == 1

def test_device_seed_reproducibility_and_independence():
    a = DevicePopulation.monte_carlo(n_devices=32, device_seed=123)
    b = DevicePopulation.monte_carlo(n_devices=32, device_seed=123)
    c = DevicePopulation.monte_carlo(n_devices=32, device_seed=124)
    for name in ("z_rp", "z_tmr", "z_t50", "z_a", "z_w", "z_eta"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))
    assert not np.array_equal(a._matrix(), c._matrix())

def test_device_values_are_fixed_across_pulses():
    population = DevicePopulation.monte_carlo(20, 9)
    before = population._matrix().copy()
    MTJModel().simulate([1.4, 1.6, 1.7], [0.9, 1.1, 1.3], population=population, switching_seed=2)
    np.testing.assert_array_equal(population._matrix(), before)
    with pytest.raises(ValueError):
        population.z_rp[0] = 0.0

def test_psw_bounds_width_gate_and_voltage_range():
    model = MTJModel()
    pop = DevicePopulation.monte_carlo(2000, 4)
    values = model.evaluate_psw(1.7, 1.2, AP, pop)
    assert np.all(values >= model.parameters.psw_min)
    assert np.all(values <= model.parameters.psw_max)
    assert np.all(model.evaluate_psw(1.7, 0.699, AP, pop) == 0.0)
    assert model.evaluate_psw(1.27, 1.2, AP) == 0.0
    assert model.evaluate_psw(3.0, 1.2, AP) == model.evaluate_psw(1.92, 1.2, AP)

def test_evaluate_psw_matches_single_pulse_simulation():
    model = MTJModel()
    pop = DevicePopulation.monte_carlo(128, 33)
    expected = model.evaluate_psw(1.6, 1.2, AP, pop)
    run = model.simulate(1.6, 1.2, AP, pop, switching_seed=44)
    np.testing.assert_allclose(run.probabilities[0], expected, rtol=2e-14)

def test_switching_seed_controls_outcomes_only():
    model = MTJModel()
    pop = DevicePopulation.monte_carlo(1000, 81)
    a = model.simulate(1.6, 1.2, AP, pop, switching_seed=7)
    b = model.simulate(1.6, 1.2, AP, pop, switching_seed=7)
    c = model.simulate(1.6, 1.2, AP, pop, switching_seed=8)
    np.testing.assert_array_equal(a.switched, b.switched)
    np.testing.assert_array_equal(a.states, b.states)
    np.testing.assert_array_equal(a.probabilities, c.probabilities)
    assert np.any(a.switched != c.switched)
    np.testing.assert_array_equal(pop._matrix(), DevicePopulation.monte_carlo(1000, 81)._matrix())

def test_numpy_and_numba_trajectory_consistency():
    p = ModelParameters()
    pop = DevicePopulation.monte_carlo(64, 17)
    rng = np.random.default_rng(8)
    voltage = rng.uniform(1.1, 2.1, (12, 64))
    width = rng.uniform(0.4, 1.8, (12, 64))
    uniforms = rng.random((12, 64))
    initial = rng.integers(0, 2, 64, dtype=np.int8)
    expected = _numpy.simulate_sequence(voltage, width, initial, pop._matrix(), uniforms, p)
    actual = _numba.simulate_sequence(voltage, width, initial, pop._matrix(), uniforms, p)
    np.testing.assert_allclose(actual[0], expected[0], rtol=2e-14, atol=2e-15)
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(actual[2], expected[2])

def test_resistance_equations_and_actual_voltage():
    p = ModelParameters(vratio_p=0.5, vratio_ap=0.5)
    model = MTJModel(p)
    pop = DevicePopulation.monte_carlo(10, 21)
    voltage = 0.3
    result = model.resistance(voltage, AP, pop)
    vv = voltage**2
    rp = 1000 * (p.rp_mu_a / (p.rp_mu_b * vv + p.rp_mu_c) + (p.rp_sig_a + p.rp_sig_b / (1 + p.rp_sig_c * vv)) * pop.z_rp)
    tmr = p.tmr_mu_a / (1 + p.tmr_mu_b * vv) + (p.tmr_sig_a + p.tmr_sig_b / (1 + p.tmr_sig_c * vv)) * pop.z_tmr
    np.testing.assert_allclose(result.tmr, np.maximum(tmr, p.tmr_min), rtol=2e-15)
    np.testing.assert_allclose(result.resistance_ohm, np.maximum(rp, p.rp_min_ohm) * (1 + np.maximum(tmr, p.tmr_min)), rtol=2e-15)
    assert model.resistance(0.3, P).resistance_ohm != model.resistance(0.6, P).resistance_ohm

def test_trajectory_and_resistance_indexing():
    model = MTJModel()
    pop = DevicePopulation.monte_carlo(50, 5)
    run = model.simulate([1.6, 1.7, 1.5], [1.0, 1.1, 1.2], AP, pop, 6, 0.1)
    assert run.probabilities.shape == run.switched.shape == (3, 50)
    assert run.states.shape == run.resistance_ohm.shape == (4, 50)
    assert run.tmr.shape == (50,)
    np.testing.assert_array_equal(run.switched, run.states[1:] != run.states[:-1])
    direct = model.resistance(0.1, run.states, pop)
    np.testing.assert_allclose(run.resistance_ohm, direct.resistance_ohm)
    rp = model.resistance(0.1, P, pop).resistance_ohm
    rap = model.resistance(0.1, AP, pop).resistance_ohm
    np.testing.assert_allclose(run.resistance_ohm[-1], np.where(run.states[-1] == AP, rap, rp))

def test_resistance_voltage_array_shapes():
    model = MTJModel()
    voltage = np.linspace(-0.5, 0.5, 11)
    nominal = model.resistance(voltage, AP)
    mc = model.resistance(voltage, AP, DevicePopulation.monte_carlo(7, 2))
    assert nominal.resistance_ohm.shape == nominal.tmr.shape == (11,)
    assert mc.resistance_ohm.shape == mc.tmr.shape == (11, 7)

@pytest.mark.parametrize("voltage,width,state,z", [
    (1.6, 1.2, AP, np.zeros(6)),
    (1.5, 0.95, P, np.array([0.2, -0.4, 0.5, -0.7, 0.3, 1.1])),
    (2.2, 1.7, AP, np.array([-0.8, 0.6, -0.2, 0.4, -0.5, -1.0])),
])
def test_fixed_va_probability_parity_points(voltage, width, state, z):
    p = ModelParameters()
    pop = DevicePopulation(*(np.array([value]) for value in z))
    actual = MTJModel(p).evaluate_psw(voltage, width, state, pop)[0]
    assert actual == pytest.approx(_literal_va_psw(voltage, width, state, z, p), rel=2e-14, abs=2e-15)

def test_fixed_va_rp_tmr_rap_parity_point():
    p = ModelParameters()
    z = np.array([0.35, -0.7, 0, 0, 0, 0])
    pop = DevicePopulation(*(np.array([value]) for value in z))
    voltage = 0.42
    vv = voltage**2
    rp = max(1000 * (p.rp_mu_a / (p.rp_mu_b * vv + p.rp_mu_c) + (p.rp_sig_a + p.rp_sig_b / (1 + p.rp_sig_c * vv)) * z[0]), p.rp_min_ohm)
    tmr = max(p.tmr_mu_a / (1 + p.tmr_mu_b * vv) + (p.tmr_sig_a + p.tmr_sig_b / (1 + p.tmr_sig_c * vv)) * z[1], p.tmr_min)
    assert MTJModel(p).resistance(voltage, P, pop).resistance_ohm[0] == pytest.approx(rp, rel=2e-15)
    result = MTJModel(p).resistance(voltage, AP, pop)
    assert result.tmr[0] == pytest.approx(tmr, rel=2e-15)
    assert result.resistance_ohm[0] == pytest.approx(rp * (1 + tmr), rel=2e-15)

def test_optional_read_and_input_validation():
    run = MTJModel().simulate(1.6, 1.2, switching_seed=1)
    assert run.resistance_ohm is None and run.tmr is None
    with pytest.raises(ValueError):
        MTJModel().evaluate_psw(1.6, -0.1)
    with pytest.raises(ValueError):
        MTJModel().simulate(1.6, 1.2, initial_state=2)
    with pytest.raises(ValueError):
        DevicePopulation.monte_carlo(0, 1)
    assert replace(ModelParameters(), eta0=0.0).eta0 == 0.0

@pytest.mark.parametrize("m,n", [(5, 3), (3, 3), (1, 1)])
def test_scalar_state_sweep_keeps_separate_device_axis(m, n):
    model = MTJModel()
    pop = DevicePopulation.monte_carlo(n_devices=n, device_seed=10)
    voltage = np.linspace(1.4, 1.8, m)
    r = model.resistance(read_voltage_v=voltage, mtj_state=AP, population=pop)
    p = model.evaluate_psw(write_voltage_v=voltage, pulse_width_ns=1.2,
                           initial_state=AP, population=pop)
    assert r.resistance_ohm.shape == r.tmr.shape == p.shape == (m, n)
    for j, v in enumerate(voltage):
        np.testing.assert_allclose(r.resistance_ohm[j], model.resistance(v, AP, pop).resistance_ohm)
        np.testing.assert_allclose(p[j], model.evaluate_psw(v, 1.2, AP, pop))

def test_per_device_probability_and_single_pulse():
    model = MTJModel()
    pop = DevicePopulation.monte_carlo(n_devices=3, device_seed=123)
    state = np.array([P, AP, P])
    actual = model.evaluate_psw(1.6, 1.2, state, pop)
    assert actual.shape == (3,)
    for i in range(3):
        expected = _literal_va_psw(1.6, 1.2, state[i], pop._matrix()[i], model.parameters)
        assert actual[i] == pytest.approx(expected, abs=2e-15)
    run = model.simulate(1.6, 1.2, state, pop, switching_seed=99)
    np.testing.assert_allclose(actual, run.probabilities[0], rtol=2e-14, atol=2e-15)

@pytest.mark.parametrize("voltage_shape", [(), (3,), (1, 3), (4, 1), (4, 3)])
def test_device_aware_resistance_and_probability_broadcast(voltage_shape):
    model = MTJModel()
    pop = DevicePopulation.monte_carlo(n_devices=3, device_seed=4)
    state = np.array([[P, AP, P], [AP, P, AP], [AP, AP, P], [P, P, AP]])
    voltage = np.linspace(1.4, 1.8, int(np.prod(voltage_shape))).reshape(voltage_shape)
    r = model.resistance(voltage, state, pop)
    ps = model.evaluate_psw(voltage, 1.2, state, pop)
    assert r.resistance_ohm.shape == ps.shape == (4, 3)
    assert r.tmr.shape == ((3,) if voltage.ndim == 0 else (4, 3))
    broadcast_voltage = np.broadcast_to(voltage, state.shape)
    for k in range(4):
        for i in range(3):
            device = DevicePopulation(*(np.array([value]) for value in pop._matrix()[i]))
            expected = model.resistance(broadcast_voltage[k, i], state[k, i], device)
            assert r.resistance_ohm[k, i] == pytest.approx(expected.resistance_ohm[0])
            assert ps[k, i] == pytest.approx(
                _literal_va_psw(broadcast_voltage[k, i], 1.2, state[k, i],
                               pop._matrix()[i], model.parameters), abs=2e-15)

def test_explicit_sweep_over_per_device_states():
    pop = DevicePopulation.monte_carlo(n_devices=3, device_seed=7)
    p = MTJModel().evaluate_psw(np.array([1.4, 1.6])[:, None], 1.2, [P, AP, P], pop)
    assert p.shape == (2, 3)
    with pytest.raises(ValueError):
        MTJModel().resistance(np.array([0.1, 0.2]), [P, AP, P], pop)

@pytest.mark.skipif(not _numba.NUMBA_AVAILABLE, reason="requires actual Numba JIT")
def test_real_jit_public_seed_consistency(monkeypatch):
    pop = DevicePopulation.monte_carlo(n_devices=128, device_seed=4)
    model = MTJModel(replace(ModelParameters(), vratio_p=0.7, vratio_ap=0.85))
    arguments = dict(write_voltage_v=[1.0, 1.6, 2.2, 1.5],
                     pulse_width_ns=[1.1, 0.699, 1.2, 1.4],
                     initial_state=np.arange(128) % 2, population=pop,
                     switching_seed=987, read_voltage_v=0.1)
    jit = model.simulate(**arguments)
    assert _numba._evolve_numba.nopython_signatures
    monkeypatch.setattr(_numba, "NUMBA_AVAILABLE", False)
    fallback = model.simulate(**arguments)
    np.testing.assert_allclose(jit.probabilities, fallback.probabilities, rtol=2e-14, atol=2e-15)
    np.testing.assert_array_equal(jit.switched, fallback.switched)
    np.testing.assert_array_equal(jit.states, fallback.states)
    np.testing.assert_array_equal(jit.resistance_ohm, fallback.resistance_ohm)

def test_unavailable_numba_uses_numpy(monkeypatch):
    monkeypatch.setattr(_numba, "NUMBA_AVAILABLE", False)
    test_numpy_and_numba_trajectory_consistency()
    run = MTJModel().simulate(1.6, 1.2, switching_seed=5)
    assert run.states.shape == (2, 1)

@pytest.mark.parametrize("parameters", [
    ModelParameters(),
    replace(ModelParameters(), vratio_p=0.7, vratio_ap=0.85, psw_min=0.01, psw_max=0.95),
])
def test_sequence_gates_clamps_and_mc_guards(parameters):
    z = np.array([[0, 0, 0, 0, -10, 0], [0, 0, 2, -100, 0, 2],
                  [0, 0, -2, 1, 0, -2]], dtype=float)
    voltage = np.array([[1.0]*3, [1.6]*3, [1.6]*3, [4.0]*3])
    widths = np.array([[1.2]*3, [0.699]*3, [0.7]*3, [2.0]*3])
    initial = np.array([P, AP, P], dtype=np.int8)
    uniforms = np.random.default_rng(56).random(voltage.shape)
    expected = _numpy.simulate_sequence(voltage, widths, initial, z, uniforms, parameters)
    actual = _numba.simulate_sequence(voltage, widths, initial, z, uniforms, parameters)
    np.testing.assert_allclose(actual[0], expected[0], rtol=2e-14, atol=2e-15)
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(actual[2], expected[2])
    assert np.all(actual[0][:2] == 0)

def _benchmark():
    # Keep the optional developer benchmark here to preserve package structure.
    import platform
    import statistics
    import time
    if not _numba.NUMBA_AVAILABLE:
        raise SystemExit("Numba is unavailable; no JIT benchmark was performed.")
    import numba
    print(f"Python {platform.python_version()}, NumPy {np.__version__}, Numba {numba.__version__}")
    print("End-to-end simulate, no read; cold median of 3 fresh cache-disabled dispatchers; warm median of 7.")
    print("pulses devices numpy_ms cold_ms warm_ms speedup max_psw_error states_equal")
    original_dispatcher = _numba._evolve_numba
    try:
        for pulses, devices in [(20, 100), (100, 100), (100, 1000),
                                (100, 10000), (1000, 1000)]:
            pop = DevicePopulation.monte_carlo(n_devices=devices, device_seed=123)
            model = MTJModel()
            args = dict(write_voltage_v=np.linspace(1.4, 1.85, pulses),
                        pulse_width_ns=np.linspace(0.7, 1.5, pulses),
                        initial_state=AP, population=pop, switching_seed=456)
            def measure():
                start = time.perf_counter()
                result = model.simulate(**args)
                return time.perf_counter() - start, result
            _numba.NUMBA_AVAILABLE = False
            measure()
            numpy_times = [measure()[0] for _ in range(7)]
            _, expected = measure()
            _numba.NUMBA_AVAILABLE = True
            cold_times = []
            for _ in range(3):
                _numba._evolve_numba = numba.njit(cache=False)(_numba._evolve_python)
                elapsed, actual = measure()
                cold_times.append(elapsed)
                assert _numba._evolve_numba.nopython_signatures
            measure()
            warm_times = [measure()[0] for _ in range(7)]
            np.testing.assert_allclose(actual.probabilities, expected.probabilities, rtol=2e-14, atol=2e-15)
            np.testing.assert_array_equal(actual.switched, expected.switched)
            np.testing.assert_array_equal(actual.states, expected.states)
            numpy_t, cold_t, warm_t = map(statistics.median, (numpy_times, cold_times, warm_times))
            error = np.max(np.abs(actual.probabilities - expected.probabilities))
            print(f"{pulses} {devices} {numpy_t*1000:.3f} {cold_t*1000:.3f} "
                  f"{warm_t*1000:.3f} {numpy_t/warm_t:.2f} {error:.3g} True", flush=True)
    finally:
        _numba._evolve_numba = original_dispatcher
        _numba.NUMBA_AVAILABLE = True

if __name__ == "__main__":
    import sys
    if sys.argv[1:] == ["--benchmark"]:
        _benchmark()
    else:
        raise SystemExit("Usage: python tests/test_model.py --benchmark")
