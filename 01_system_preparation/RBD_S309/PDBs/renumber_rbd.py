from simtk.openmm import app
import argparse

# Read filename
parser = argparse.ArgumentParser(description='renumber')
parser.add_argument('name', type=str, help='name of file for which to renumber')
parser.add_argument('--split', default=False, action='store_true')
args = parser.parse_args()
pdb = app.PDBFile(args.name)

# Rename residues
current_start = 530
for residue in pdb.topology.residues():
    if residue.chain.id == 'X':
        current_start += 1
        residue.id = str(current_start) # Update the residue id

app.PDBFile.writeFile(pdb.topology, pdb.positions, open(f"{args.name[:-4]}_renumbered.pdb", "w"), keepIds=True)
