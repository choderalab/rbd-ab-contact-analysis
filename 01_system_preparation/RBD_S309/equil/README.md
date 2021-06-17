### Author: Ivy Zhang

# Equilibration Protocol

The RBD:S309 (`7jx3_s309.inpcrd`, `7jx3_s309.prmtop`) system were equilibrated for 20 ns. Check `equil_NPT_20ns.py` for complete simulation details.

* Virtual bonds were added across chains that should be kept together to prevent imaging issues
* Performed in the NPT ensemble using the `Openmmtools` Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.

Outputs are in `output/`

Note: `.dcd` not included since it is large and unnecessary

Since `tleap` merges all chains into one and renumbers residues starting from 1, `output/equilibrated.pdb` was renumbered using `correct_equilibrated_s309.py` and
 saved as `output/equilibrated.pdb`. The original, non-renumbered structure was moved
to `output/equilibrated_old.pdb`.

Writing PDB files with `OpenMM` in `correct_equilibrated_s309.py` doesn't retain glycan `CONECT` records. To retain these (useful for imaging after a production run) `correct_equilibrated_s309_mda.py` was used to create `output/X/equilibrated_for_imaging.pdb` (where `X` is a particular RUN number). This file was used for imaging / analysing runs returned from Folding@Home trajectories.
