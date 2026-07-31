# Netlist-Based Testbench
	Standalone Spectre testbenches for the VC-MTJ compact model.


	## Quick Start

		Run the typical-corner testbench with:
		./run_tb_typical.sh

		Run the Monte Carlo mismatch testbench with:
		./run_tb_mc.sh

		Each script removes the previous results, runs the corresponding Spectre netlist, and opens ViVA.
		Use the Results Browser in ViVA to open the generated psf directory and select the desired waveforms.

	## Model Corners
		tb_typical.scs: uses the mtj_t corner for a standard transient simulation.
		tb_mc.scs: uses the mtj_mc corner for one mismatch Monte Carlo run with an embedded transient simulation.

	## Files
		- tb_typical.scs: typical-corner transient testbench
		- run_tb_typical.sh: typical-corner simulation script
		- tb_mc.scs: one-run mismatch Monte Carlo testbench
		- run_tb_mc.sh: Monte Carlo simulation script
