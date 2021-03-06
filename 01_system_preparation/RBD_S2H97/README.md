### Author: Ivy Zhang  

*Summary:* The RBD:S2H97 complex (fully glycosylated) was saved as two structures (see `./PDBs/5_rbd_renumbered_final.pdb` `./PDBs/5_s2h97_renumbered_final.pdb`). These structures were then read into `tleap` for solvation and parametrization, because the [GLYCAM_06j-1 forcefield](https://pubmed.ncbi.nlm.nih.gov/17849372/) is not available in `OpenMM 7.4.2`. The `.inpcrd` and `.prmtop` files generated by `tleap` were then read into `OpenMM 7.4.2` and equilibrated for 20 ns. The output files from equilibration were then used to seed Folding@home simulations.

## Glycan structures

Glycosylation sites and patterns were determined by [Elisa Fadda](https://www.maynoothuniversity.ie/people/elisa-fadda). The recent [Science paper](https://science.sciencemag.org/content/early/2020/05/01/science.abb9983) by Watanabe *et al.* was used to propose the most likely glycosylation pattern.

Pre-equilibrated glycan structures are stored in the `01_system_preparation/representative_glycan_structures` directory. 

# System Preparation

## Obtain starting structure
This pipeline uses [7M7W](https://www.rcsb.org/structure/7M7W) as the starting structure. Chains S, C, and D (and waters within 3 angstroms of the complex) were kept, everything else was removed. 

After visual inspection, we discovered a missing loop in the constant domain of S2H97. We built the missing loop in (and add the N343 glycan) using ISOLDE.

Location: `./PDBs/0_glycos.pdb`, `./PDBs/1_glycos_renamed.pdb` (same as `0_glycos.pdb` but with residues renamed according to AMBER residue naming convention)

Glycosylation patterns were kept the same as in the original structures, as a reminder these are:

*RBD*
* N343 (FA2G2)

Note: We did not build in the glycan at N331 because the residues 324-331 are missing in 7M7W (meaning the there is no GlcNac at N331), therefore, we wouldn't be able to model the N331 glycan in properly.

## Capping the PDB file

We then took the structure and capped termini with missing residues. If the chain terminus matched the SEQRES (i.e. had no missing terminal residues), we did not cap it. When selecting rotamers for the cap, we chose the rotamer that most closely matched the rotamer already present.

Location (capped): `./PDBs/2_glycos_renamed_noh_capped.pdb`

## Cleaning the PDB file

We first split the structure such that the RBD and S2H97 antibody were separate PDB files. These are located at `.PDBs/3_rbd.pdb` and `.PDBs/3_s2h97.pdb`.

We then renumbered each PDB file using `.PDBs/renumber.py`. These are located at: `.PDBs/4_rbd_renumbered.pdb` and `.PDBs/4_s2h97_renumbered.pdb`.

Next, we  manually edited the spacing of lines containing `CCYX` to have one less space before the `CCYX` and one more space after the `CCYX`.
- Originally, it looked like this: `N[space][space][space]CCYXL[space]215` (this was throwing an error when we tried to parametrize in `tleap`)
- We modified the spacing to look like this: `N[space][space]CCYX[space]L[space]215`

Finally, we added `TER` cards between the glycan residues and moved TER cards for NME residues to be _after_ the last NME residue line (necessary for `tleap`) using `./PDBs/correct_TER_NME_CYS_CONECT.py`. We also moved water lines to the end of the PDB.

The final structures are saved in `./PDBs` as:

* `5_rbd_renumbered_final.pdb`
* `5_s2h97_renumbered_final.pdb`

## Running tleap

Location: `./run_tleap/`

We used `tleap` to solvate and parametrize the system. 

### Addition of solvent and ions

Within the `rbd_s2h97_tleap.in` file, the `solvateoct` command solvates the system using a truncated octahedron (with 22 angstrom padding) and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined we use the methodology from OpenMM:

```python
from math import floor
numWaters = 77458
numPositive = 0
numNegative = 0 
totalCharge = 11
ionicStrength = 0.15

if totalCharge > 0:
    numNegative += totalCharge
else:
    numPositive -= totalCharge

numIons = (numWaters - numPositive - numNegative) * ionicStrength / (55.4)  # Pure water is about 55.4 molar (depending on temperature)
numPairs = int(floor(numIons + 0.5))
numPositive += numPairs
numNegative += numPairs
```

The `numWaters` variable was determined by running the commands in the `rbd_s2h97_tleap.in` script before adding ions. The above method gave `numPositive` was 210 and `numNegative` was 221.

### Running tleap

The `tleap` commands we ran are present in `rbd_s2h97_tleap.in`.

This produced fully solvated (with ions) systems. Check the `leap.log` files for full descriptions of the output. The `.prmtop` and `.inpcrd` files could now be used for minimisation and equilibration. The generated system can be viewed here: `./run_tleap/rbd_s2h97_tleap_out.pdb`, though note that `tleap` combines all chains into 1, so it is a bit difficult to use PyMOL to manipulate the system.

## Equilibration

Once prepared the system was equilibrated here: `./equil`
The equilibrated PDB were renumbered such that chains and residue ids are numbered according to the original PDB structures using `.equil/correct_equilibrated_s2h97.py`

