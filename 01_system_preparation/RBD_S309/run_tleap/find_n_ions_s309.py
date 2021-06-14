from math import floor
numWaters = 82796
numPositive = 0
numNegative = 0 
totalCharge = 12
ionicStrength = 0.15

if totalCharge > 0:
    numNegative += totalCharge
else:
    numPositive -= totalCharge

numIons = (numWaters - numPositive - numNegative) * ionicStrength / (55.4)  # Pure water is about 55.4 molar (depending on temperature)
numPairs = int(floor(numIons + 0.5))
numPositive += numPairs
numNegative += numPairs

print(f"numPositive: {numPositive}")
print(f"numNegative: {numNegative}")
