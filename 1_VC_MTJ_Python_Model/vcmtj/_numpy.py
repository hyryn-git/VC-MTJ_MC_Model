"""Vectorized NumPy equations for the VC-MTJ compact model."""
import numpy as np

def evaluate_psw(write_voltage_v, pulse_width_ns, mtj_state, z, p):
    """Evaluate VA-equivalent switching probabilities on broadcast arrays."""
    ratio = np.where(mtj_state == 1, p.vratio_ap, p.vratio_p)
    vfit_raw = write_voltage_v / ratio
    vfit = np.minimum(vfit_raw, p.vfit_max)
    vv = vfit * vfit
    y0 = p.y0_a * vv + p.y0_b * vfit + p.y0_c
    t50 = p.t50_a * vv + p.t50_b * vfit + p.t50_c
    amp = p.amp_a * vv + p.amp_b * vfit + p.amp_c
    width = p.w_a * vv + p.w_b * vfit + p.w_c
    t50 = t50 + (p.st50_a * vv + p.st50_b * vfit + p.st50_c) * z[..., 2]
    amp = amp + (p.samp_a * vv + p.samp_b * vfit + p.samp_c) * z[..., 3]
    width = width + p.sw * z[..., 4]
    amp = np.maximum(amp, p.amp_min)
    width = np.maximum(width, p.w_min_ns)
    x = 2.0 * np.log(4.0) / width * (pulse_width_ns - t50)
    term = np.exp(-np.abs(x))
    sigmoid = np.where(x >= 0.0, 1.0 / (1.0 + term), term / (1.0 + term))
    common = np.clip(y0 + amp * sigmoid, p.psw_min, p.psw_max)
    eta = p.eta0 + p.seta * z[..., 5]
    scale = np.exp(np.minimum(np.where(mtj_state == 1, 0.5, -0.5) * eta, 709.0))
    with np.errstate(over="ignore", invalid="ignore"):
        result = -np.expm1(np.log1p(-common) * scale)
    result = np.clip(result, p.psw_min, p.psw_max)
    invalid = (vfit_raw < p.vfit_min) | (pulse_width_ns < p.t_valid_min_ns)
    return np.where(invalid, 0.0, result)

def resistance(read_voltage_v, mtj_state, z, p):
    """Evaluate Rp, TMR, and state-dependent resistance at actual voltage."""
    vv = read_voltage_v * read_voltage_v
    rp_mu = p.rp_mu_a / (p.rp_mu_b * vv + p.rp_mu_c)
    rp_sig = p.rp_sig_a + p.rp_sig_b / (1.0 + p.rp_sig_c * vv)
    rp = np.maximum(1000.0 * (rp_mu + rp_sig * z[..., 0]), p.rp_min_ohm)
    tmr_mu = p.tmr_mu_a / (1.0 + p.tmr_mu_b * vv)
    tmr_sig = p.tmr_sig_a + p.tmr_sig_b / (1.0 + p.tmr_sig_c * vv)
    tmr = np.maximum(tmr_mu + tmr_sig * z[..., 1], p.tmr_min)
    return np.where(mtj_state == 1, rp * (1.0 + tmr), rp), tmr

def simulate_sequence(write_voltage_v, pulse_width_ns, initial_state, z, uniforms, p):
    """Evolve a pulse-by-device trajectory using NumPy calculations."""
    n_pulses, n_devices = write_voltage_v.shape
    probabilities = np.empty((n_pulses, n_devices))
    switched = np.empty((n_pulses, n_devices), dtype=np.bool_)
    states = np.empty((n_pulses + 1, n_devices), dtype=np.int8)
    states[0] = initial_state
    current = initial_state.copy()
    for pulse in range(n_pulses):
        probability = evaluate_psw(write_voltage_v[pulse], pulse_width_ns[pulse], current, z, p)
        did_switch = uniforms[pulse] < probability
        current = np.where(did_switch, 1 - current, current).astype(np.int8)
        probabilities[pulse], switched[pulse], states[pulse + 1] = probability, did_switch, current
    return probabilities, switched, states
