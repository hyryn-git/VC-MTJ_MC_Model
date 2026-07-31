# Schematic-Based Testbench

A Cadence Virtuoso schematic-based testbench for the VC-MTJ compact model. The saved ADE XL setup supports both typical transient and mismatch Monte Carlo simulations.

## Quick Start

Add this directory as a Virtuoso library and open the saved `tb/adexl` view.

Update the model-library path to the local location of `model_file/vc_mtj_core.lib` when necessary.

### Typical Transient Simulation

Select **Single Run, Sweeps and Corners**, enable the typical corner, and run the configured transient analysis.

### Monte Carlo Simulation

Select **Monte Carlo Sampling**, enable the MC corner, and run the simulation with mismatch variation enabled.

## Model Corners

- `typical`: Uses the `mtj_t` model section for nominal transient simulation.
- `MC`: Uses the `mtj_mc` model section for mismatch Monte Carlo simulation.

## Library Contents

- `tb/schematic`: Schematic testbench.
- `tb/adexl`: Saved ADE XL setup supporting both typical transient and mismatch Monte Carlo simulations.
- `vc_mtj_3t`: VC-MTJ symbol and Spectre interface.
