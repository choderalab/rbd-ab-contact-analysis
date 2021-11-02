import MDAnalysis as mda
import argparse

# Read filenames
parser = argparse.ArgumentParser(description='correct equilibrated.pdb')
parser.add_argument('-mutant_id', dest='mutant_id', type=int, help='the mutant index')
parser.add_argument('-job_id', dest='job_id', type=int, help="the job index")
parser.add_argument('-prmtop_file', dest='prmtop_file', type=str, help='the prmtop file created from tleap')
parser.add_argument('-inpcrd_file', dest='inpcrd_file', type=str, help='the inprcd file created from tleap')
parser.add_argument('-ref_file', dest='ref_file', type=str, help='a reference PDB to get the correct unit cell dimensions (e.g. ./output/equilibrated.pdb)')

args = parser.parse_args()

mutant_id = args.mutant_id
job_id = args.job_id
prmtop_file = args.prmtop_file
inpcrd_file = args.inpcrd_file
ref_file = args.ref_file

# Load in the topology from tleap output files 
u = mda.Universe(prmtop_file, inpcrd_file)

u_dim = mda.Universe(ref_file)
dimensions = u_dim.dimensions

# RBD
# find the end of the RBD, this will change with each mutated system
rbd_end_sel = u.select_atoms("resname NME and name HH33")
print(rbd_end_sel)
# since there are two NME residues, we take the first one as this will be the end of the RBD chain
rbd_end_index = rbd_end_sel[0].index 
print(rbd_end_index)

#rbd = u.select_atoms("index 0-3077")
rbd = u.select_atoms(f"index 0-{str(rbd_end_index)}")
new_rbd_resids = [i for i in range(331, 531)]
rbd.residues.resids = new_rbd_resids

#rbd_glycans = u.select_atoms("index 3078-3311")
rbd_glycans_start = rbd_end_index + 1
rbd_glycans_end = rbd_glycans_start + 233 
rbd_glycans = u.select_atoms(f"index {rbd_glycans_start}-{rbd_glycans_end}")
new_rbd_glycan_resids = [i for i in range(1, len(rbd_glycans.residues.resids) + 1)]
rbd_glycans.residues.resids = new_rbd_glycan_resids

# S309 Chain A
s309_a_start = rbd_glycans_end + 1
s309_a_end = s309_a_start + 3393# + 1
s309_a = u.select_atoms(f"index {str(s309_a_start)}-{str(s309_a_end)}")
#s309_a = u.select_atoms("index 3312-6705")
new_s309_a_resids = [i for i in range(1, len(s309_a.residues.resids) + 1)]
s309_a.residues.resids = new_s309_a_resids

# S309 Chain B
s309_b_start = s309_a_end + 1
s309_b_end = s309_b_start + 3211# + 1
s309_b = u.select_atoms(f"index {str(s309_b_start)}-{str(s309_b_end)}")
#s309_b = u.select_atoms("index 6706-9917")
new_s309_b_resids = [i for i in range(1, len(s309_b.residues.resids) + 1)]
s309_b.residues.resids = new_s309_b_resids

# Solvent
solvent = u.select_atoms("resname Na+ or resname Cl-")
new_solvent_ion_resids = [i for i in range(1, len(solvent.residues.resids) + 1)]
solvent.residues.resids = new_solvent_ion_resids

# Create the new system by merging each universe
new_system = mda.Merge(rbd, rbd_glycans, s309_a, s309_b, solvent)

# Name each chain
new_system.segments.segids = ['R', 'X', 'A', 'B', 'Y']

new_system.dimensions = dimensions

# Write out the new system
#new_system.atoms.write("./output/equilibrated_for_imaging.pdb")
new_system.atoms.write(f"./{mutant_id}/{job_id}/equilibrated_for_imaging.pdb")
