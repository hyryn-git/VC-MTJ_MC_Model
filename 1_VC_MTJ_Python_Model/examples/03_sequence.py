"""A pulse sequence and its state/resistance trajectory."""
from vcmtj import AP, DevicePopulation, MTJModel

model = MTJModel()
devices = DevicePopulation.monte_carlo(
    n_devices=8,
    device_seed=123,
)
run = model.simulate(
    write_voltage_v=[1.6, 1.7, 1.6],
    pulse_width_ns=[1.0, 1.1, 1.2],
    initial_state=AP,
    population=devices,
    switching_seed=456,
    read_voltage_v=0.1,
)
print("probabilities (pulse x device):\n", run.probabilities)
print("switched (pulse x device):\n", run.switched)
print("states[0] is initial; later rows follow each pulse:\n", run.states)
print("resistance aligned with every state row (ohm):\n", run.resistance_ohm)
