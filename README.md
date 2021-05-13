## Thermosetting supramolecular polymerisation of compartmentalised DNA fibres with stereo-sequence and length control
*Michael Dore, Tuan Trinh, Marlo Zorman, Donatien de Rochambeau, Pengfei Xu, Xin Luo, Jacob Remington, Violeta Toader, Jianing Li, Hanadi Sleiman*

DNA nanostructures are highly addressable and compatible with biological systems, but often require hundreds of unique strands for their assembly. On the other hand, nature can assemble complex structures from identical building blocks by compartmentalising the process: assembling molecules into sub-components and bringing these together across multiple length scales. Inspired by this process, we report DNA-polymer conjugates that assemble through a unique heat-driven hierarchical mechanism to form fibres displaying blocks of different DNA sequences along their axis of polymerisation . Ions are manipulated to pre-assemble short cylindrical micelle ‘monomers’ that then non-covalently fuse at 90 °C. These one-dimensional ‘thermoset’ DNA-fibres retain the compartmentalisation programmed in the pre-assembled segments, allowing nanoscale organisation of gold nanoparticles. Length control over fibre segments is also achieved through gel purification of pre-assembled cylindrical micelles. Importantly, the DNA-polymers are sequence-defined, allowing us to show that stereochemical sequence of the hydrophobic core can be amplified into distinctive morphological traits in the DNA-fibres. Molecular dynamics simulations are used to model the structure of the supramolecular fibres and elucidate how stereochemical sequence manifests itself in different fibre-fibre interactions.


### Description
This repository is a supplement to the methods and SI of (DOI) and allows the user to replicate all molecular dynamics work. It includes simulation input structures, figure scripts/structures, and analysis scripts/data. Trajectory coordinates are not included here but are available upon request. Simulations and analysis were done using the Schrodinger Suite, and figures were made with Pymol and Matplotlib.

### Files
*analysis_utils.py* -- Contains analysis utility scripts. Functions used to fit a circle to backbone atoms of the hydrophobic core were adapted from [Fitting a Circle to Cluster of 3D Points](https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/). Contains some analysis code that was not used in the final manuscript but which may be useful.

*circle_fitting.py* -- Used to find the conformations from the trajectory of a single amphiphile that best represent a monomeric unit. It produces a time-plot of scores and a csv file listing top conformation candidates (frames) and scores.

*sasa_analysis.py* -- Calculates and plots the non-DNA solvent accessible surface area (SASA) of amphiphiles. Used for both single molecules and tubes

*replica_sasa_analysis.py* -- Calculates and plots the non-DNA solvent accessible surface area (SASA) range of all amphiphile replicas.

*pymol_fig.txt* -- Used to set representations, colors, lighting, etc, of structures in Pymol.

*structs/* -- Contains all input structures used for simulations

### Dependencies
* Schrodinger
* Pymol
* Numpy
* Matplotlib
