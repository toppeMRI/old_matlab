## TOPPE Matlab code

Current version is toppe v2b

+toppe: package (namespace) with basic toppe functions such as writemod.m, plotmod.m, playseq.m.

   +toppe/+utils/: scripts for loading raw data from P-files, and other scripts.

   +toppe/+utils/+spiral/: create stack-of-spirals readout module, and reconstruct with reconSoS.m

   +toppe/+utils/+rf/: make slice-selective SLR pulse with makeslr.m

./pulseq/: convert to/from Pulseq file format

./examples/: complete TOPPE sequence examples [TODO]

All code in this repository is **open-source** under the terms of the GNU LGPL 2.0 license (see LICENSE for details).

For an overview of the TOPPE platform, including the 'official' TOPPE user guide, see [https://toppemri.github.io/](https://toppemri.github.io/).

To scan, you'll also need the 'EPIC' source code for the TOPPE binary executable; see [https://toppemri.github.io/](https://toppemri.github.io/) for details.

For questions about TOPPE, contact Jon-Fredrik Nielsen at jfnielse@umich.edu.


## toppe v2 (old) files:
The 'lib' subfolder in this respository contains basic MATLAB functions for creating and viewing TOPPE sequence files.

The 'examples' folder contains several tested sequences and is a good place to start.

