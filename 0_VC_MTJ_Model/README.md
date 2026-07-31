# Field-Free VC-MTJ Compact Model
	A Verilog-A model capturing hybrid dynamics and supporting Monte Carlo simulation in Cadence Virtuoso and Spectre.

	## Quick Start & More Info
		Two testbench implementations are provided:
		- Schematic-based testbench: Virtuoso schematic-based testbench for typical-corner and Monte Carlo simulations, read 'testbench/schematic_based/README.md'
		- Netlist-based testbench: standalone Spectre testbenches for typical-corner and Monte Carlo simulations, read 'testbench/netlist_based/README.md'

	## Files Organization

	- model_file
	  - vc_mtj_core.lib: model library and simulation corners
	  - vc_mtj_core.mdl: Spectre subcircuit wrapper
	  - vc_mtj_core.va: Verilog-A compact model

	- testbench
	  - schematic_based
		- tb: schematic testbench containing the schematic view and the saved ADE XL setup for both typical transient and mismatch Monte Carlo simulations
		- vc_mtj_3t: VC-MTJ symbol and Spectre interface
	  - netlist_based
		- tb_typical.scs: typical-corner transient testbench
		- run_tb_typical.sh: typical-corner simulation script
		- tb_mc.scs: one-run mismatch Monte Carlo testbench
		- run_tb_mc.sh: Monte Carlo simulation script
