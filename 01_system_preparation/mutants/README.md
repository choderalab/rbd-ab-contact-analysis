### Author: Ivy Zhang

The mutant RBDs (glycosylated at N343) were prepped from the WT RBD using `mutate.py` and then each of the mutant RBDs was loaded into tleap with WT S309 to generate the RBD:S309 mutant complexes.

# System Preparation

After mutating WT RBD at P337 to all other 19 amino acids using `mutate.py`,  we added `TER` cards between the glycan residues and moved TER cards for NME residues to be _after_ the last NME residue line (necessary for `tleap`) using `./PDBs/correct_TER_NME_CYS_CONECT_waters.py`. We also moved water lines to the end of the PDB.

The final structures are saved in `./PDBs` as:

- `7jx3_rbd_{mutant}_final.pdb`

## Running tleap

Location: `./run_tleap/`

The mutant PDBs were copied into numbered dirs `1-19` in alphabetical order by amino acid mutant. E.g. 1 is P337A, 2 is P337C, 19 is P337Y, etc. 

{'P337A': 1, 'P337C': 2, 'P337D': 3, 'P337E': 4, 'P337F': 5, 'P337G': 6, 'P337H': 7, 'P337I': 8, 'P337K': 9, 'P337L': 10, 'P337M': 11, 'P337N': 12, 'P337Q': 13, 'P337R': 14, 'P337S': 15, 'P337T': 16, 'P337V': 17, 'P337W': 18, 'P337Y': 19}

Each of the copied PDBs was named `7jx3_rbd_mutant_final.pdb` in the corresponding numbered dir.

We used `tleap` to solvate the systems and parametrize them. We did this separately for each of the systems by running tleap with `7jx3_s309_tleap_template.in` first to determine the number of solvent molecules and charge and then used this info to edit the tleap file using `edit_tleap.py` to use the right number of positive/negative ions. The edited `{mutant}_tleap.in` files were saved in the numbered dirs and then used to run tleap again. 

### Addition of solvent and ions

Within the `[mutant]_tleap.in` files, the `solvateoct` command solvates the system (with 22 angstrom padding) and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined we use the methodology from OpenMM:

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

The `numWaters` variable was determined by running the commands in the `[mutant]_tleap.in` script before adding ions.

We solvated with 22 angstroms padding because we were previously observing the complex turning in its rectangular solvent box such that the long axis of the complex would align with the short axis of the box, causing the protein to interact with the image of itself. 

### Running tleap

The `tleap` commands we ran are present in `[mutant]_tleap.in`.

This produced fully solvated (with ions) systems. Check the `[mutant]/leap.log` files for full descriptions of the output. The `.prmtop` and `.inpcrd` files could now be used for minimisation and equilibration. The generated system can be viewed here: `./run_tleap/[mutant]_out.pdb`, though note that `tleap` combines all chains into 1, so it is a bit difficult to use PyMOL to manipulate the system.

## Equilibration

Once prepared the system was equilibrated here: `./equil`
The equilibrated PDBs were renumbered such that chains and residue ids are numbered according to the original PDB structures using `.equil/*/correct_equilibrated*.py`
To render the equilibrated PDB properly for analysis, we ran `.equil/correct_equilibrated_s309_mda.py`.
