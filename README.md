# Field-Free VC-MTJ Monte Carlo Model

A field-free VC-MTJ compact modeling framework supporting Monte Carlo
simulation, with a Verilog-A implementation for Cadence/Spectre circuit
simulation and a Python implementation for fast standalone statistical
simulation.

## Repository contents

- [0_VC_MTJ_VerilogA_Model](0_VC_MTJ_VerilogA_Model/README.md)
  provides the Verilog-A implementation for Cadence Virtuoso and Spectre,
  including schematic-based and netlist-based testbenches.
- [1_VC_MTJ_Python_Model](1_VC_MTJ_Python_Model/README.md)
  provides standalone Monte-Carlo switching, pulse-sequence simulation, and
  R-V/TMR analysis using NumPy with optional Numba acceleration.
- [Manual_Field_Free_VC-MTJ_v1.0.1.pdf](Manual_Field_Free_VC-MTJ_v1.0.1.pdf)
  is the user manual for the Verilog-A/Cadence implementation.
