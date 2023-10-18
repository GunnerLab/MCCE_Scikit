"""Module `constants` holds physical and thermodynamical constants."""

PH2KCAL = 1.364  # pH "unit" to Kcal/mol
KCAL2KT = 1.688  # Kcal/mol to kT
KJ2KCAL = 0.239  # Kj/mol to Kcal/mol
ROOMT = 298.15  # default simulation temp [CC: mismatched precision]

VALID_MC_METHODS = [
    "MONTERUNS",
]
