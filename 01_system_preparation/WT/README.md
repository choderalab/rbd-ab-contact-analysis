### Author: Ivy Zhang

The RBD:S309:S304:S2H14 complex (glycosylated at N343 in the RBD) was prepped as a whole, then split up into separate structures (RBD:S309, RBD:S304, RBD:S309:S304) at the end.

Upon analyzing the trajectories, we realized there was an imaging issue due to the solvent box being too small, so we modified the box to a truncated octahedron with 22 angstrom padding (RBD:S309 and RBD:S304). Later, we realized that that there was an issue with glycan residue naming (we had named a residue 0FA, but it should be 0fA, which caused the stereochemistry of the methyl group to be incorrect). We fixed the naming for RBD:S309 only.  

## Glycan structures

Glycosylation sites and patterns were determined by [Elisa Fadda](https://www.maynoothuniversity.ie/people/elisa-fadda). The recent [Science paper](https://science.sciencemag.org/content/early/2020/05/01/science.abb9983) by Watanabe *et al.* was used to propose the most likely glycosylation pattern.

Pre-equilibrated glycan structures are stored in the `systems/representative_glycan_structures` directory. 

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

We then took the structure outputted from ISOLDE and capped temini with missing residues. If the chain terminus matched the SEQRES (i.e. had no missing terminal residues), we did not cap it. When selecting rotamers for the cap, we chose the rotamer that most closely matched the rotamer already present.

Location (capped): `./PDBs/03_7jx3_rebuild_isolde_without_H_capped.pdb`

We then deleted the S2H14 antibody (chains C and D) from the structure.

Location: `./PDBs/04_7jx3_rebuild_isolde_without_H_capped_without_S2H14.pdb`

## Cleaning the PDB file

We first split the structure such that each chain is a separate PDB file (except S304 chains H and L, which must be kept as one PDB because there is a disulfide bond between the C-termini of chains H and L). These are located at `.PDB/05_*.pdb`. 

We then renumbered the RBD PDB file using `.PDBs/renumber_rbd.py`. This is located at: `.PDB/06_7jx3_rbd_renumbered.pdb`.

Next, we  manually edited the spacing of lines containing `CCYX` to have one less space before the `CCYX` and one more space after the `CCYX`.
- Originally, it looked like this: `N[space][space][space]CCYXL[space]215` (this was throwing an error when we tried to parametrize in `tleap`)
- We modified the spacing to look like this: `N[space][space]CCYX[space]L[space]215`

Finally, we added `TER` cards between the glycan residues and moved TER cards for NME residues to be _after_ the last NME residue line (necessary for `tleap`) using `./PDBs/correct_TER_NME_CYS_CONECT.py`. We also moved water lines to the end of the PDB.

The final structures are saved in `./PDBs` as:

* `7jx3_rbd_renumbered_final.pdb`
* `7jx3_s309_A_final.pdb`
* `7jx3_s309_B_final.pdb`
* `7jx3_s304_H_L_final.pdb`

## Running tleap

Location: `./run_tleap/`

We used `tleap` to solvate the systems and parametrize them. We did this separately for each of the 3 systems, each of which has a separate directory in `./run_tleap`. Below, we use `system_name` in place of `7jx3_s309`, `7jx3_s304`, and `7jx3_s309_s304` to outline the relevant input/output files for running `tleap`.

### Addition of solvent and ions

Within the `[system_name]_tleap.in` files, the `solvateBox` command solvates the system and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined we use the methodology from OpenMM (using `7jx3_s309_s304` as an example):

```python
from math import floor
numWaters = 70255
numPositive = 0
numNegative = 0 
totalCharge = 14
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

The `numWaters` variable was determined by running the commands in the `[system_name]_tleap.in` script before adding ions. The above method gave `numPositive` (i.e. Na+) as 190 and `numNegative` (i.e. Cl-) as 204.

For `7jx3_s309`: `numWaters` was 41624 and `totalCharge` was 12. From the computation, `numPositive` was 224 and `numNegative` was 236.

For `7jx3_s304`: `numWaters` was 54631 and `totalCharge` was 6. From the computation, `numPositive` was 148 and `numNegative` was 162.

Note that we solvated `7jx3_s309` and `7jx3_s304` using 22 angstroms padding with a truncated octahderal box. This was to ensure no overlap between neighbouring periodic images.

## Running tleap

The `tleap` commands we ran are present in `[system_name]_tleap.in`.

This produced fully solvated (with ions) systems. Check the `[system_name]/leap.log` files for full descriptions of the output. The `.prmtop` and `.inpcrd` files could now be used for minimisation and equilibration. The generated system can be viewed here: `./run_tleap/[system_name]_out.pdb`, though note that `tleap` combines all chains into 1, so it is a bit difficult to use PyMOL to manipulate the system.

## Equilibration

Once prepared the system was equilibrated here: `./equil`
The equilibrated PDBs were renumbered such that chains and residue ids are numbered according to the original PDB structures using `.equil/*/correct_equilibrated*.py`
To render the equilibrated PDB properly in pymol, it must be read into `MDAnalysis` and written out: 
```python
import MDAnalysis as mda
ref = mda.Universe("equilibrated.pdb")
protein = ref.select_atoms("all")
protein.write("protein.pdb")
```
