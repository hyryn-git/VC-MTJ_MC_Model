"""One nominal VC-MTJ write/read operation."""
from vcmtj import AP, MTJModel

model = MTJModel()
run = model.simulate(
    write_voltage_v=1.6,
    pulse_width_ns=1.2,
    initial_state=AP,
    switching_seed=456,
    read_voltage_v=0.1,
)
print("Psw:", run.probabilities[0, 0])
print("switched:", run.switched[0, 0])
print("final state:", run.states[-1, 0])
print("final resistance (ohm):", run.resistance_ohm[-1, 0])
print("TMR:", run.tmr[0])
