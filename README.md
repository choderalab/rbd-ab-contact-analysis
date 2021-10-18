# rbd-ab-contact-analysis
Code and workflow for running MD simulations on Folding@home to analyze RBD:Ab contacts

Publication: https://doi.org/10.1038/s41586-021-03807-6

## License
* This software is licensed under the [MIT license](https://opensource.org/licenses/MIT) - a copy of this license is provided as `SOFTWARE_LICENSE`
* The data in this repository is made available under the Creative Commons [CC0 (“No Rights Reserved”) License](https://creativecommons.org/share-your-work/public-domain/cc0/) - a copy of this license is provided as `DATA_LICENSE`

## Manifest

* `01_system_preparation` - Contains all scripts and PDBs before and after structure preparation, system parameterization, and equilibration. Also contains the representative glycan structures added to RBD:S309 and RBD:S2H97
* `02_analysis_scripts` - Contains the scripts used to analyze RBD:S309 and RBD:S2H97 contacts in the trajectories

## Contributors

* Ivy Zhang
* William G. Glass
* Tristan Croll
* Aoife M. Harbison
* Elisa Fadda
* John D. Chodera

# Simulation details

**Structure preparation**
The RBD:S309 complex was constructed from PDB ID 7JX3 (Chains A, B, and R). 7JX3 was first refined using [ISOLDE](https://isolde.cimr.cam.ac.uk/) to better fit the experimental electron density using detailed manual inspection. Refinement included adjusting several rotamers, flipping several peptide bonds, fixing several weakly resolved waters, and building in a missing four-residue-long loop. Though the N343 glycan N-Acetylglucosamine (NAG) was present in 7JX3, ISOLDE was used to construct a complex glycan at N343. The full glycosylation pattern was determined from [Shajahan et al.](http://doi.org/10.1093/glycob/cwaa101) and [Watanabe et al.](http://doi.org/10.1126/science.abb9983) The glycan structure used for N343 (FA2G2) corresponds to the most stable conformer obtained from multi microsecond molecular dynamics (MD) simulations of cumulative sampling [Harbison et al.](http://doi.org/10.1093/glycob/cwy097). The base NAG residue in FA2G2 was aligned to the corresponding NAG stub in the RBD:S309 model and any resulting clashes were refined in ISOLDE. The same process was repeated for the RBD:S2H97 crystal structure.

**System solvation and parametrization**
The refined glycosylated complexes were prepared for simulation using the [AmberTools20](https://ambermd.org/AmberTools.php) tleap suite. All relevant disulfide bridges were specified as well as covalent connectivity within each glycan structure. The glycosylated protein was parameterized with the Amber ff14SB ([Maier et al., 2015](https://doi.org/10.1021/acs.jctc.5b00255)) and GLYCAM_06j-1 ([Kirschner et al., 2008](https://doi.org/10.1002/jcc.20820)) force fields. The system was solvated using the TIP3P rigid water model ([Jorgensen et al., 1983](https://doi.org/10.1063/1.445869)) in a cubic box with 1.5 nm solvent padding on all sides. The solvated system was then minimally neutralized with 0.15 M NaCl using the Li/Merz ion parameters of monovalent ions for the TIP3P water model (12-6 normal usage set) (Li et al., 2015).

**System equilibration**
Each system was energy-minimized with an energy tolerance of 10 kJ mol−1 and and equilibrated five times independently using the OpenMMTools 0.20.0 (https://github.com/choderalab/openmmtools) BAOAB Langevin integrator ([Leimkuhler and Matthews, 2013](https://doi.org/10.1063/1.4802990)) for 20 ns in the NPT (p=1 atm, T = 310 K) ensemble with a timestep of 4.0 femtoseconds, a collision rate of 1.0 picoseconds-1, and a relative constraint tolerance of 1 ✕ 10−5. Hydrogen atom masses were set to 4.0 amu by transferring mass from connected heavy atoms, bonds to hydrogen were constrained, and center of mass motion was not removed. Pressure was controlled by a molecular-scaling Monte Carlo barostat with an update interval of 25 steps. Non-bonded interactions were treated with the Particle Mesh Ewald method ([Darden et al., 1993](https://doi.org/10.1063/1.464397)) using a real-space cutoff of 1.0 nm and the OpenMM default relative error tolerance of 0.0005, with grid spacing selected automatically. The simulations were subsequently packaged to seed for production simulation on Folding@home ([Shirts and Pande, 2000](https://science.sciencemag.org/content/290/5498/1903.full), [Zimmerman et al., 2020](https://doi.org/10.1101/2020.06.27.175430)). Default parameters were used unless noted otherwise.

**Folding@home simulation**
The equilibrated structures (5 for each complex) were then used to initiate parallel distributed MD simulations on Folding@home ([Shirts and Pande, 2000](https://science.sciencemag.org/content/290/5498/1903.full), [Zimmerman et al., 2020](https://doi.org/10.1101/2020.06.27.175430)). Simulations were run with OpenMM 7.4.2 (Folding@home core22 0.0.13). Production simulations used the same Langevin integrator as the NpT equilibration described above. 5000 and 4985 independent MD simulations were generated on Folding@home for the complexes RBD:S309 and RBD:S2H97, respectively. Conformational snapshots (frames) were stored at an interval of 1 ns/frame for subsequent analysis. The final datasets contained 1.1 ms and 623.7 µs of aggregate simulation time for RBD:S309 and RBD:S2H97, respectively. This trajectory dataset (without solvent) is be available at the MolSSI COVID-19 Molecular Structure and Therapeutics Hub:
- RBD:S309: https://covid.molssi.org//simulations/#foldinghome-simulations-of-the-sars-cov-2-spike-rbd-bound-to-monoclonal-antibody-s309
- RBD:S2H97: https://covid.molssi.org//simulations/#foldinghome-simulations-of-the-sars-cov-2-spike-rbd-bound-to-monoclonal-antibody-s2h97 
