"""One write pulse applied to a fixed Monte-Carlo population."""
import numpy as np
from vcmtj import AP, DevicePopulation, MTJModel

model = MTJModel()
devices = DevicePopulation.monte_carlo(
    n_devices=10_000,
    device_seed=123,
)
run = model.simulate(
    write_voltage_v=1.6,
    pulse_width_ns=1.2,
    initial_state=AP,
    population=devices,
    switching_seed=456,
    read_voltage_v=0.1,
)
print("mean Psw:", run.probabilities[0].mean())
print("switched fraction:", run.switched[0].mean())
print("final AP fraction:", np.mean(run.states[-1] == AP))
print("mean final resistance (ohm):", run.resistance_ohm[-1].mean())
