from ortools.sat.python import cp_model
import numpy as np
import time
from utils import read
from isotope_pattern import peak_isotope

PROTON_MASS = 1.0078250319


class VarArraySolutionPrinter(cp_model.CpSolverSolutionCallback):
    '''
    Print intermediate solutions
    '''

    def __init__(self, variables):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__all_variables = []
        self.__variables = variables
        self.__solution_count = 0

    def on_solution_callback(self):
        self.__solution_count += 1
        self.__all_variables.append([[v._IntVar__var.name, self.Value(v)] for v in self.__variables])
        # for v in self.__variables:
        #     print('%s=%i' % (v, self.Value(v)), end=' ')
        # print()

    def solution_count(self):
        return self.__solution_count
    
    def get_variables(self):
        return self.__all_variables, self.__solution_count

def feasible_set_search(compounds, compound_masses, compound_maximum, compound_minimum, \
        peak_mass, tolerance, multi_protein=False, primaries=None, min_primaries=None, \
            max_primaries=None, metal_idx=None, max_per_metal=None, adducts=None, max_adducts=None, \
                valence=None):
    '''
    Returns feasible integer solutions https://developers.google.com/optimization/cp/cp_solver#all_solutions
    '''
    # Creates the model
    model = cp_model.CpModel()

    # Creates the variables
    x = [model.NewIntVar(int(compound_minimum[c]), int(compound_maximum[c]), c) for c in compounds]
    
    # Add multi-protein constraints
    if multi_protein:
        z = [model.NewBoolVar('used_protein_{c}') for c in compounds if primaries[c]]
        j = 0
        for i, c in enumerate(compounds):
            if primaries[c]:
                model.Add(x[i] == 0).OnlyEnforceIf(z[j].Not())
                model.Add(x[i] > 0).OnlyEnforceIf(z[j])
                j += 1
        model.Add(sum(z) >= min_primaries)
        model.Add(sum(z) <= max_primaries)

    # Maximum adducts in potential species
    z2 = [model.NewBoolVar('used_adduct_{c}') for c in compounds if adducts[c]]
    j = 0
    for i, c in enumerate(compounds):
        if adducts[c]:
            model.Add(x[i] == 0).OnlyEnforceIf(z2[j].Not())
            model.Add(x[i] > 0).OnlyEnforceIf(z2[j])
            j += 1
    model.Add(sum(z2) <= max_adducts)

    # Add ratio constraints (max per metal)
    for j, c in enumerate(compounds):
        if not np.isnan(max_per_metal[c]):
            model.Add(x[j] <= int(max_per_metal[c])*x[metal_idx])

    # Valence
    model.Add(sum([x[i] for i, c in enumerate(compounds) if not np.isnan(max_per_metal[c])]) <= valence*x[metal_idx])

    # Add min/max constraints
    model.Add(sum([x[i]*compound_masses[c] for i, c in enumerate(compounds)]) <= peak_mass + tolerance)
    model.Add(sum([x[i]*compound_masses[c] for i, c in enumerate(compounds)]) >= peak_mass - tolerance)

    # Create a solver and solve.
    solver = cp_model.CpSolver()
    solution_printer = VarArraySolutionPrinter(x)
    # Enumerate all solutions.
    solver.parameters.enumerate_all_solutions = True
    # Solve.
    status = solver.Solve(model, solution_printer)
    # print('Status = %s' % solver.StatusName(status))
    print('Number of solutions found: %i' % solution_printer.solution_count())
    
    return solution_printer.get_variables()


# if __name__ == '__main__':

#     print('Reading data...')

#     base_path = "Data/"
#     fn = "Ubiquitin_plusC_1in100_000001"
#     unbound_file_path = base_path + "Deconvoluted Spectra/Ubi_1in100_broad band_000001.xlsx"
#     compounds_file_path = base_path + "Compound Constraints/Compounds_CisOxTrans_nlp.xlsx"
#     bound_file_path = base_path + "Deconvoluted Spectra/" + fn + ".xlsx"

#     bound_df, unbound_df, compounds = read(bound_file_path, unbound_file_path, compounds_file_path)

#     exp_masses = bound_df['m/z'].to_numpy()
#     exp_abundance = bound_df['I'].to_numpy()
#     formulas = compounds["Formula"].to_numpy()
#     charges = compounds["Charge"].to_numpy()
#     max_amount = compounds["Max"].to_numpy()
#     min_amount = compounds["Min"].to_numpy()
#     masses = np.vectorize(peak_isotope)(formulas) - np.dot(PROTON_MASS, charges)

#     masses = (masses*100000).astype(int)

#     masses_dict = dict(zip(formulas, masses))
#     max_amount_dict = dict(zip(formulas, max_amount))
#     min_amount_dict = dict(zip(formulas, min_amount))

#     primaries_dict = dict(zip(formulas, compounds['Primaries'].notna().to_numpy()))
#     min_primaries = 3

#     print('Searching for solutions...')
#     start = time.time()
#     print('Time started:', time.ctime(start))

#     for peak in [8775,8792,8811,8828,9002,9019,9037,9055]:
#         print(f'Peak {peak}: ', end='')
#         solutions = feasible_set_search(formulas, masses_dict, max_amount_dict, min_amount_dict, \
#             peak_mass=peak*100000, tolerance=400000, primaries=primaries_dict, min_primaries=min_primaries)
    
#     end = time.time()
#     print('Time ended:', time.ctime(end))
#     print('Elapsed (minutes):', str((end-start)/60.))