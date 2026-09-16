# Field-Free VC-MTJ Python Model

## Overview

This package is the Python implementation of a field-free VC-MTJ probabilistic
compact model for standalone Monte Carlo and stochastic switching simulation.
It supports nominal simulation, device-level Monte Carlo simulation,
stochastic single-pulse switching, pulse-sequence simulation, resistance,
R-V, and TMR analysis, and optional Numba acceleration.

## File Organization

- `vcmtj/`
  - `__init__.py`: Package entry point and public exports.
  - `model.py`: Public API for VC-MTJ simulation.
  - `parameters.py`: Calibrated model parameters.
  - `_numpy.py`: NumPy numerical implementation.
  - `_numba.py`: Optional Numba-accelerated simulation backend.
- `examples/`
  - Example scripts demonstrating typical model usage.
- `tests/`
  - Regression and model-validation tests.
- `pyproject.toml`
  - Python package and dependency configuration.

## Installation

From the Python model directory:

```bash
pip install .
```

For optional Numba acceleration:

```bash
pip install ".[numba]"
```

## Quick Start & Examples

Three example scripts are provided:

- `01_basic_simulation.py`: One nominal VC-MTJ write/read simulation,
  including switching probability, stochastic switching result, final state,
  and resistance.
- `02_monte_carlo.py`: Device-level Monte Carlo simulation with fixed device
  variations, stochastic switching, and final resistance distribution.
- `03_sequence.py`: Multi-pulse stochastic simulation with probability,
  switching, state, and resistance trajectories.

Open and run the corresponding example for the desired simulation flow.

The main user interface is:

- `model.simulate(...)`: Stochastic write simulation with optional resistance
  evaluation.

Two auxiliary interfaces are provided:

- `model.evaluate_psw(...)`: Switching-probability evaluation without
  stochastic state update.
- `model.resistance(...)`: Standalone resistance and TMR evaluation,
  including R-V analysis.

`DevicePopulation.monte_carlo(...)` creates a fixed device-level Monte Carlo
population. `device_seed` controls device-to-device variation, while
`switching_seed` controls stochastic switching trials.

## Notes

- `P = 0` and `AP = 1`.
- `write_voltage_v` and `read_voltage_v` represent actual voltages across
  the MTJ terminals.
- The switching model internally applies the calibrated voltage mapping
  consistent with the reference Verilog-A model.
- The model is calibrated over a fitting-voltage range of 1.6-2.4 V.
- The minimum valid effective write-pulse width is 0.7 ns.
- Numba acceleration is optional and workload-dependent.
- Detailed model definitions and advanced usage will be provided separately
  in the model manual.
