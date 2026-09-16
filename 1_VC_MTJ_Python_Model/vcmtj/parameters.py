"""Calibrated parameters for the VC-MTJ compact model."""
from dataclasses import asdict, dataclass
import math

@dataclass(frozen=True)
class ModelParameters:
    """Physical fit coefficients and numerical guards used by the model.

    Voltage coefficients use volts, time coefficients use nanoseconds, and
    resistance coefficients produce kOhm before conversion to ohms. Defaults
    reproduce ``vc_mtj_core.va``; Monte-Carlo variation is held separately.
    """
    # Actual terminal write voltage to calibration-domain voltage mapping.
    vratio_p: float = 0.8
    vratio_ap: float = 0.8
    vfit_min: float = 1.6
    vfit_max: float = 2.4
    # P-state resistance mean and standard-deviation model.
    rp_mu_a: float = 114.1122817
    rp_mu_b: float = 0.180017789
    rp_mu_c: float = 1.108333982
    rp_sig_a: float = 2.45705674
    rp_sig_b: float = 3.93397863
    rp_sig_c: float = 0.916399026
    # Tunnel-magnetoresistance mean and standard-deviation model.
    tmr_mu_a: float = 0.572198649
    tmr_mu_b: float = 0.505290180
    tmr_sig_a: float = 0.001562837
    tmr_sig_b: float = 0.091027265
    tmr_sig_c: float = 0.261932936
    # Nominal switching-probability sigmoid parameters.
    y0_a: float = 0.1005543629
    y0_b: float = -0.3642518154
    y0_c: float = 0.3375512292
    t50_a: float = -0.0268571501
    t50_b: float = -0.2023076311
    t50_c: float = 1.7689730294
    amp_a: float = -1.0679461125
    amp_b: float = 5.0983018830
    amp_c: float = -5.2246142286
    w_a: float = -0.2484296130
    w_b: float = 0.8697977039
    w_c: float = -0.4387266590
    eta0: float = 0.2030650335

    # Device-to-device variation applied to the Psw model.
    st50_a: float = 0.1544664595
    st50_b: float = -0.5386098902
    st50_c: float = 0.5807208485
    samp_a: float = 0.0814842091
    samp_b: float = -0.2770523218
    samp_c: float = 0.3186725567
    sw: float = 0.0920203869
    seta: float = 0.6685692209
    # Validity limits and numerical guards.
    t_valid_min_ns: float = 0.7
    rp_min_ohm: float = 1.0
    tmr_min: float = 0.0
    w_min_ns: float = 0.001
    amp_min: float = 1e-6
    psw_min: float = 0.0
    psw_max: float = 0.999999

    def __post_init__(self) -> None:
        """Validate that a parameter set is finite and numerically usable."""
        values = asdict(self)
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError("all model parameters must be finite")
        for name in ("vratio_p", "vratio_ap", "rp_min_ohm", "w_min_ns", "amp_min"):
            if values[name] <= 0.0:
                raise ValueError(f"{name} must be positive")
        if self.vfit_min <= 0.0 or self.vfit_max < self.vfit_min:
            raise ValueError("require 0 < vfit_min <= vfit_max")
        if self.t_valid_min_ns < 0.0 or self.tmr_min < 0.0:
            raise ValueError("validity and clipping minima must be nonnegative")
        if not 0.0 <= self.psw_min <= self.psw_max < 1.0:
            raise ValueError("require 0 <= psw_min <= psw_max < 1")
