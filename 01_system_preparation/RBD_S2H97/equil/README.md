### Author: Ivy Zhang

# Equilibration Protocol

The RBD:S2H97 (`rbd_s2h97.inpcrd`, `rbd_s2h97.prmtop`) system were equilibrated for 20 ns. Check `equil_NPT_20ns.py` for complete simulation details.

* Virtual bonds were added across chains that should be kept together to prevent imaging issues
* Performed in the NPT ensemble using the `Openmmtools` Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.

Outputs are in `output/`

Note: `.dcd` not included since it is large and unnecessary

Since `tleap` merges all chains into one and renumbers residues starting from 1, `output/equilibrated.pdb` was renumbered using `correct_equilibrated_s309.py` and
 saved as `output/equilibrated.pdb`. The original, non-renumbered structure was moved
to `output/X/equilibrated_old.pdb` (where `X` is a particular RUN number).

Writing PDB files with `OpenMM` in `correct_equilibrated_s309.py` doesn't retain glycan `CONECT` records.

