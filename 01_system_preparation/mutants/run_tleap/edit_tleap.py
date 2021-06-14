import argparse
import re
from math import floor

# Read filename
parser = argparse.ArgumentParser(description='edit tleap file')
parser.add_argument('out_name', type=str, help='name of out file to check')
parser.add_argument('in_name', type=str, help='name of in file to edit')
args = parser.parse_args()

# Retrieve charge and num of waters
with open(args.out_name, "r") as f:
    lines_out = f.readlines()

for line in lines_out:
    if "Total unperturbed charge" in line:
        charge = float(line.split(":")[1].strip('\n'))
    if "residues" in line:
        result = re.findall(r"\d*", line)
        result_filtered = [r for r in result if r]
        num_waters = int(result_filtered[0])

# Compute number of ions
numWaters = num_waters
numPositive = 0
numNegative = 0 
totalCharge = charge
ionicStrength = 0.15

if totalCharge > 0:
    numNegative += totalCharge
else:
    numPositive -= totalCharge

numIons = (numWaters - numPositive - numNegative) * ionicStrength / (55.4)  # Pure water is about 55.4 molar (depending on temperature)
numPairs = int(floor(numIons + 0.5))
numPositive += numPairs
numNegative += numPairs
print(f"num positive: {numPositive}")
print(f"num negative: {numNegative}")

# Edit tleap file
with open(args.in_name, "r") as f:
    lines_in = f.readlines()

new_lines = []
for line in lines_in:
    if "addionsrand complex" in line:
        line = f"addionsrand complex Na+ {int(numPositive)} Cl- {int(numNegative)}\n"
    new_lines.append(line)

with open("7jx3_s309_tleap.in", 'w') as f:
    f.writelines(new_lines)
