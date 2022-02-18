from pyopenms import *
import numpy as np
# https://pyopenms.readthedocs.io/en/latest/aasequences.html


# Get fragments and their average weights
formula_str = "DFPIANGER"
all_fragments_str = sum([[formula_str[i:i+j] for i in range(len(formula_str)-j+1)] for j in range(1,len(formula_str)+1)], [])
all_fragments = [AASequence.fromString(f) for f in all_fragments_str]

all_weights = [f.getAverageWeight() for f in all_fragments]
y_ions = [[f.getAverageWeight(Residue.ResidueType.YIon, 0) for f in all_fragments]] +\
    [[f.getAverageWeight(Residue.ResidueType.YIon, j) for f in all_fragments] for j in [1, 2, 3]]
x_ions = [[f.getAverageWeight(Residue.ResidueType.XIon, 0) for f in all_fragments]] +\
    [[f.getAverageWeight(Residue.ResidueType.XIon, j) for f in all_fragments] for j in [1, 2, 3]]
a_ions = [[f.getAverageWeight(Residue.ResidueType.AIon, 0) for f in all_fragments]] +\
    [[f.getAverageWeight(Residue.ResidueType.AIon, j) for f in all_fragments] for j in [1, 2, 3]]


seq = AASequence.fromString("C378H630N105O118S1NH3")

suffix = seq.getSuffix(3) # y3 ion "GER"
print("="*35)
print("y3 ion sequence:", suffix)
y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2) # y3++ ion
suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 # CORRECT
suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0 # CORRECT
suffix.getMonoWeight(Residue.ResidueType.BIon, 2) / 2.0 # INCORRECT

print("y3 mz:", suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 )
print("y3 molecular formula:", y3_formula)