import pymol
from pymol import cmd
import sys
import argparse

# https://sourceforge.net/p/pymol/mailman/message/27979284/
# https://sourceforge.net/p/pymol/mailman/message/11671708/

# Read in arguments
parser = argparse.ArgumentParser(description='run pymol mutagenesis')
parser.add_argument('pdb', type=str, help='path to input file')
parser.add_argument('selection', type=str, help='e.g. "R/337/"')
parser.add_argument('mutant', type=str, help='amino acid three letter code, e.g. ALA')
parser.add_argument('wt', type=str, help='for the last part of the outfile e.g. "P337"')
args = parser.parse_args()
# Beware, the error message you get is because pymol unsuccessfully tries 
#to open the command line arguments as files.

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S'}

# for amino_acid in d.keys():

print(f"Mutating {args.mutant}...")

# Launch pymol session
pymol.pymol_argv = ["pymol", "-qc"] + sys.argv[1:]
pymol.finish_launching()

# Load RBD (no solvent)
cmd.load(args.pdb)

# Mutate
cmd.wizard("mutagenesis")
cmd.do("refresh_wizard")
cmd.get_wizard().set_mode(args.mutant)
cmd.get_wizard().do_select(args.selection) 

# Select the rotamer
cmd.frame(1)

# Apply the mutation
cmd.get_wizard().apply()
cmd.set_wizard()

# Save
cmd.save(f"7jx3_rbd_{args.wt}{d[args.mutant]}.pdb")

# Wait for everything to be done
cmd.refresh() 

# # Remove object
# cmd.remove(args.pdb[:-4])
