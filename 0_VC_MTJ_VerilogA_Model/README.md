# Field-Free VC-MTJ Compact Model

A Verilog-A model capturing hybrid dynamics and supporting Monte Carlo simulation in Cadence Virtuoso and Spectre.

## Quick Start & More Information

Two testbench implementations are provided:

- [Schematic-based testbench](testbench/schematic_based/README.md): Virtuoso schematic-based testbench for typical-corner and Monte Carlo simulations.
- [Netlist-based testbench](testbench/netlist_based/README.md): Standalone Spectre testbenches for typical-corner and Monte Carlo simulations.

## File Organization

- `model_file/`
  - `vc_mtj_core.lib`: Model library and simulation corners.
  - `vc_mtj_core.mdl`: Spectre subcircuit wrapper.
  - `vc_mtj_core.va`: Verilog-A compact model.

- `testbench/`
  - `schematic_based/`
    - `tb/`: Schematic testbench containing the schematic view and the saved ADE XL setup for both typical transient and mismatch Monte Carlo simulations.
    - `vc_mtj_3t/`: VC-MTJ symbol and Spectre interface.
  - `netlist_based/`
    - `tb_typical.scs`: Typical-corner transient testbench.
    - `run_tb_typical.sh`: Typical-corner simulation script.
    - `tb_mc.scs`: One-run mismatch Monte Carlo testbench.
    - `run_tb_mc.sh`: Monte Carlo simulation script.
