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
- [Manual_Field_Free_VC-MTJ_v1.1.pdf](Manual_Field_Free_VC-MTJ_v1.1.pdf)
  is the user manual for the Verilog-A/Cadence/Python implementation.

## Associated Publication

This software is associated with the following paper, accepted for presentation at the IEEE International Electron Devices Meeting (IEDM) 2026:

**Monte-Carlo Modeling and Experimental Characterization of Hybrid Dynamics for Field-Free Voltage-Controlled Magnetic Tunnel Junctions**

Yaoru Hou, Yantong Di, Xuezhao Wu, Hongjie Yin, Haoran Du, Zeying Ding, Albert Lee, Di Wu, Jayong Wang, Zheng Zhu, Zhihua Xiao, Hao Cai, and Qiming Shao.

*IEEE International Electron Devices Meeting (IEDM), 2026.*

Yaoru Hou and Yantong Di contributed equally. The corresponding authors are Zheng Zhu, Hao Cai, and Qiming Shao.

The paper has been accepted but has not yet been presented or published online. DOI, page numbers, and final proceedings information will be added after official publication.

## Citation

If you use this software in your research, please cite this repository using the metadata in [CITATION.cff](CITATION.cff). Please also cite the associated IEDM 2026 paper when its final bibliographic information becomes available.

## License

This project is licensed under the BSD 3-Clause License. See [LICENSE](LICENSE) for the full license text.

Copyright (c) 2026 HKUST, SEU and Quantsyn company.

