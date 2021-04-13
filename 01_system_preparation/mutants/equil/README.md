### Author: Ivy Zhang

# Equilibration Protocol

The RBD:S309 mutant (`7jx3_s309_mutant.inpcrd`, `7jx3_s309_mutant.prmtop`) systems were equilibrated for 20 ns. Check `equil_NPT_20ns.py` for complete simulation details.

* Virtual bonds were added across chains that should be kept together to prevent imaging issues 
* Performed in the NPT ensemble using the `Openmmtools` Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.

Output directories are in zipped.

Note: `.dcd` not included since it is large and unnecessary

Since `tleap` merges all chains into one and renumbers residues starting from 1, `*/equilibrated.pdb` was renumbered using `correct_equilibrated_*.py` and
 saved as `*/equilibrated.pdb`. The original, non-renumbered structure was moved
to `*/equilibrated_old.pdb`.

