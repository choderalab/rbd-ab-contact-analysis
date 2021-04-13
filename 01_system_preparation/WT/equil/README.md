### Author: Ivy Zhang

# Equilibration Protocol

The RBD:S309 (`7jx3_s309.inpcrd`, `7jx3_s309.prmtop`), RBD:S304 (`7jx3_s304.inpcrd`, `7jx3_s304.prmtop`), RBD:S309:S304 (`7jx3_s309_s304.inpcrd`, `7jx3_s309_s304.prmtop`) systems were equilibrated for 20 ns. Check `equil_NPT_20ns.py` for complete simulation details.

* Virtual bonds were added across chains that should be kept together to prevent imaging issues
* Performed in the NPT ensemble using the `Openmmtools` Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.

Outputs are in `output/`

Note: `.dcd` not included since it is large and unnecessary

Since `tleap` merges all chains into one and renumbers residues starting from 1, `output/equilibrated.pdb` was renumbered using `correct_equilibrated_*.py` and
 saved as `output/equilibrated.pdb`. The original, non-renumbered structure was moved
to `output/equilibrated_old.pdb`.

Writing PDB files with `OpenMM` in `correct_equilibrated_*.py` doesn't retain glycan `CONECT` records. To retain these (useful for imaging after a production run) `correct_equilibrated_*_mda.py` was used to create `output/equilibrated_for_imaging.pdb`. This file was used for imaging / analysing runs returned from Folding@Home trajectories.

# Info for Fah 

| System name     | Number of atoms | Time (s) per 10 ns | 
| --------------- | --------------- | ------------------ |
| RBD:S309        | 143492          | 24295.882          |
| RBD:S304        | 142497          | 24299.824          |
| RBD:S309:S304   | 284612          | 52308.055          |


### with octahedron
| System name     | Number of atoms | Time (s) per 20 ns |
| --------------- | --------------- | ------------------ |
| RBD:S309        | 257386          | 67636.275          |
| RBD:S304        | 244411          | 62544.364          |
