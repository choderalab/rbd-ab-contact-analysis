from simtk.openmm import app
import argparse

# Read filename
parser = argparse.ArgumentParser(description='renumber')
parser.add_argument('name', type=str, help='name of file for which to renumber')
args = parser.parse_args()
pdb = app.PDBFile(args.name)

# Renumber glycan residues
if "rbd" in args.name:
    current_start = {'X': 0, 'S': 330, 'Y': 0}
elif "s2h97" in args.name:
    current_start = {'C': 0, 'D': 0, 'Z': 0}
for residue in pdb.topology.residues():
    if residue.chain.id == 'X':
        current_start['X'] += 1
        residue.id = str(current_start['X']) # Update the residue id
    elif residue.chain.id == 'S':
        current_start['S'] += 1
        residue.id = str(current_start['S']) # Update the residue id
    elif residue.chain.id == 'Y':
        current_start['Y'] += 1
        residue.id = str(current_start['Y']) # Update the residue id
    elif residue.chain.id == 'C':
        current_start['C'] += 1
        residue.id = str(current_start['C']) # Update the residue id
    elif residue.chain.id == 'D':
        current_start['D'] += 1
        residue.id = str(current_start['D']) # Update the residue id
    elif residue.chain.id == 'Z':
        current_start['Z'] += 1
        residue.id = str(current_start['Z']) # Update the residue id


app.PDBFile.writeFile(pdb.topology, pdb.positions, open(f"{args.name[:-4]}_renumbered.pdb", "w"), keepIds=True)