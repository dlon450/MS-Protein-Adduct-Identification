from pyomo.environ import *
from isotope_pattern import objective_func
import numpy as np
import idaes

from file_io import read
from isotope_pattern import average_mass

PROTON_MASS = 1.007825
            

def solve_nlp(compounds, compound_masses, compound_maximum, compound_minimum, exp_masses, exp_abundance, peak_mass, tolerance=4):
    '''
    Nonlinear Program to find optimal match for peak
    '''
    model = ConcreteModel()
    model.x = Var(compounds, initialize=0, domain=NonNegativeIntegers)
    model.z = Var(domain=Binary)

    ############ needs to be an expression ############
    # model.obj = Objective(objective_func(''.join(model.x[c].value*c for c in compounds), \
    #     exp_masses, exp_abundance, peak_mass), sense=minimize)
    ###################################################
    model.obj = Objective(expr = 0)

    model.constraints = ConstraintList()
    for c in compounds:
        model.constraints.add(model.z * compound_minimum[c] <= model.x[c])
        model.constraints.add(model.x[c] <= model.z * compound_maximum[c])
    model.constraints.add(sum([model.x[c] * compound_masses[c] for c in compounds]) <= peak_mass + tolerance)
    model.constraints.add(sum([model.x[c] * compound_masses[c] for c in compounds]) >= peak_mass - tolerance)
    results = SolverFactory('ipopt').solve(model, tee=True)
    # solver.options['timelimit'] = 100
    return results, ' + '.join(model.x[c].value*c for c in compounds)


if __name__ == '__main__':

    base_path = "Data/"
    fn = "Ubiquitin_plusC_1in100_000001"
    unbound_file_path = base_path + "Deconvoluted Spectra/Ubi_1in100_broad band_000001.xlsx"
    compounds_file_path = base_path + "Compound Constraints/Compounds_CisOxTrans_nlp.xlsx"
    bound_file_path = base_path + "Deconvoluted Spectra/" + fn + ".xlsx"

    bound_df, unbound_df, compounds = read(bound_file_path, unbound_file_path, compounds_file_path)
    exp_masses = bound_df['m/z'].to_numpy()
    exp_abundance = bound_df['I'].to_numpy()
    formulas = compounds["Formula"].to_numpy()
    charges = compounds["Charge"].to_numpy()
    max_amount = compounds["Max"].to_numpy()
    min_amount = compounds["Min"].to_numpy()
    masses = np.vectorize(average_mass)(formulas) - np.dot(PROTON_MASS, charges)

    masses_dict = dict(zip(formulas, masses))
    max_amount_dict = dict(zip(formulas, max_amount))
    min_amount_dict = dict(zip(formulas, min_amount))
    
    results, species = solve_nlp(formulas, masses_dict, max_amount_dict, min_amount_dict, \
        exp_masses, exp_abundance, peak_mass=8775, tolerance=4)