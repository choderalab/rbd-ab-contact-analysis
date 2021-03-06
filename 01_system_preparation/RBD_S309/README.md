### Author: Ivy Zhang  

*Summary:* The RBD:S309 complex (fully glycosylated) was saved as three monomeric structures (see `./PDBs/7jx3_rbd_renumbered_final.pdb` `./PDBs/7jx3_s309_A_final.pdb` and `./PDBs/7jx3_s309_B_final.pdb`). These monomeric structures were then read into `tleap` for solvation and parametrization, because the [GLYCAM_06j-1 forcefield](https://pubmed.ncbi.nlm.nih.gov/17849372/) is not available in `OpenMM 7.4.2`. The `.inpcrd` and `.prmtop` files generated by `tleap` were then read into `OpenMM 7.4.2` and equilibrated for 20 ns. The output files from equilibration were then used to seed Folding@home simulations.

## Glycan structures

Glycosylation sites and patterns were determined by [Elisa Fadda](https://www.maynoothuniversity.ie/people/elisa-fadda). The recent [Science paper](https://science.sciencemag.org/content/early/2020/05/01/science.abb9983) by Watanabe *et al.* was used to propose the most likely glycosylation pattern.

Pre-equilibrated glycan structures are stored in the `01_system_preparation/representative_glycan_structures` directory. 

# System Preparation

## Obtain starting structure
This pipeline uses a version of the [7jx3 structure](https://www.rcsb.org/structure/7JX3) rebuilt (flipped some rotamers and fixed poorly resolved waters) by Tristan Croll.

Location: `./PDBs/01_7jx3_rebuild.pdb`

Glycosylation patterns were kept the same as in the original structures, as a reminder these are:

*RBD*
* N343 (FA2G2)

Note: We did not build in the glycan at N331 because the residues 328-331 are missing in 7jx3 (meaning the there is no GlcNac at N331), therefore, we wouldn't be able to model the N331 glycan in properly.

## ISOLDE

After visual inspection, we discovered a missing loop in the constant domain of S309. In order to build in the missing loop, we used ISOLDE to copy and paste the residues from another human IgG antibody [4jhw](https://www.rcsb.org/structure/4jhw), which has the same sequence for the chain with the missing loop.

Next, we used ISOLDE to build in the glycan, which involved:

1) aligning the FA2G2 (open) structure (`systems/representative_glycan_structures/fa2g2_n_open.pdb`) to the glycan stub present in 7jx3, 
2) adding missing bonds in the sugar and between the sugar and N343, 
3) deleting the glycan stub in 7jx3. 

Finally, there were severe clashes introduced when building in the glycan, so we used ISOLDE to remove these clashes and refine non-ideal rotamers.

Location (final output from ISOLDE): `./PDBs/02_7jx3_rebuild_isolde_without_H.pdb`

For an in-depth discussion of steps taken within the ISOLDE package check `./PDBs/README.md`

## Capping the PDB file

We then took the structure outputed from ISOLDE and capped termini with missing residues. If the chain terminus matched the SEQRES (i.e. had no missing terminal residues), we did not cap it. When selecting rotamers for the cap, we chose the rotamer that most closely matched the rotamer already present.

Location (capped): `./PDBs/03_7jx3_rebuild_isolde_without_H_capped.pdb`

We then deleted the S2H14 antibody (chains C and D) from the structure.

Location: `./PDBs/04_7jx3_rebuild_isolde_without_H_capped_without_S2H14.pdb`

## Cleaning the PDB file

We first split the structure such that each chain is a separate PDB file. These are located at `.PDB/05_*.pdb`. 

We then renumbered the RBD PDB file using `.PDBs/renumber_rbd.py`. This is located at: `.PDBs/06_7jx3_rbd_renumbered.pdb`.

Next, we  manually edited the spacing of lines containing `CCYX` to have one less space before the `CCYX` and one more space after the `CCYX`.
- Originally, it looked like this: `N[space][space][space]CCYXL[space]215` (this was throwing an error when we tried to parametrize in `tleap`)
- We modified the spacing to look like this: `N[space][space]CCYX[space]L[space]215`

Finally, we added `TER` cards between the glycan residues and moved TER cards for NME residues to be _after_ the last NME residue line (necessary for `tleap`) using `./PDBs/correct_TER_NME_CYS_CONECT.py`. We also moved water lines to the end of the PDB.

The final structures are saved in `./PDBs` as:

* `7jx3_rbd_renumbered_final.pdb`
* `7jx3_s309_A_final.pdb`
* `7jx3_s309_B_final.pdb`

## Running tleap

Location: `./run_tleap/`

We used `tleap` to solvate and parametrize the system. 

### Addition of solvent and ions

Within the `7jx3_s309_tleap.in` files, the `solvateoct` command solvates the system using a truncated octahedron (with 22 angstrom padding) and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined we use the methodology from OpenMM:

```python
from math import floor
numWaters = 82743
numPositive = 0
numNegative = 0 
totalCharge = 12
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

The `numWaters` variable was determined by running the commands in the `7jx3_s309_tleap.in` script before adding ions. The above method gave `numPositive` was 224 and `numNegative` was 236.

### Running tleap

The `tleap` commands we ran are present in `7jx3_s309_tleap.in`.

This produced fully solvated (with ions) systems. Check the `leap.log` files for full descriptions of the output. The `.prmtop` and `.inpcrd` files could now be used for minimisation and equilibration. The generated system can be viewed here: `./run_tleap/7jx3_s309_out.pdb`, though note that `tleap` combines all chains into 1, so it is a bit difficult to use PyMOL to manipulate the system.

## Equilibration

Once prepared the system was equilibrated here: `./equil`
The equilibrated PDB were renumbered such that chains and residue ids are numbered according to the original PDB structures using `.equil/correct_equilibrated_s309.py`
To render the equilibrated PDB properly for analysis, we ran `.equil/correct_equilibrated_s309_mda.py`.
