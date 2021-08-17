# caraslab-spikesortingKS2

Pipeline for spike sorting multi-channel data acquired with TDT and Intan hardware.

Preprocessing steps convert data from native formats (Synapse or OpenEphysGUI) to .mat and/or .dat files. Files are then common median referenced and high-pass filtered, and run through Kilosort2 for sorting. 
This version runs with the modified Kilosort2 code also present in this repository. 

Required before running for the first time:
- npy-matlab (https://github.com/kwikteam/npy-matlab)
- Kilosort2: install CUDA Kilosort2 according to the developers' instructions (https://github.com/MouseLand/Kilosort) using the caraslab-Kilosort2-master.
