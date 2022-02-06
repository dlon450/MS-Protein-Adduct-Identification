import pandas as pd
import numpy as np
from itertools import combinations  
from isotope_pattern import average_mass
import time
from constraint_optimisation import feasible_set_search

PROTON_MASS = 1.007825


def feasible_set_df(compounds, peaks_mass, tolerance, multi_protein=False, min_primaries=None, precision=5):
    '''
    Returns feasible set of integer solutions 
    '''
    factor = 10.**precision
    solution_dict = {}

    formulas = compounds["Formula"].to_numpy()
    charges = compounds["Charge"].to_numpy()
    masses = np.vectorize(average_mass)(formulas) - np.dot(PROTON_MASS, charges)

    masses = (masses*factor).astype(int)

    masses_dict = dict(zip(formulas, masses))
    charges_dict = dict(zip(formulas, charges))
    max_amount_dict = dict(zip(formulas, compounds["Max"].to_numpy()))
    min_amount_dict = dict(zip(formulas, compounds["Min"].to_numpy()))

    primaries_dict = None
    if multi_protein:
        primaries_dict = dict(zip(formulas, compounds['Primaries'].notna().to_numpy()))

    print('Searching for solutions...')
    start = time.time()
    print('Time started:', time.ctime(start))

    for peak in peaks_mass:
        print(f'Peak {round(peak, 2)} ------ ', end='')
        solutions, count = feasible_set_search(formulas, masses_dict, max_amount_dict, min_amount_dict, \
            peak_mass=int(peak*factor), tolerance=int(tolerance*factor), multi_protein=multi_protein, \
                primaries=primaries_dict, min_primaries=min_primaries)

        if count != 0:
            potential_cmpds = [[[
                c[0]*c[1] for c in soln if c[1] != 0],
                sum([masses_dict[c[0]]*c[1] for c in soln if c[1] != 0])/factor,
                sum([charges_dict[c[0]]*c[1] for c in soln if c[1] != 0])] 
                for soln in solutions]
            output = list(zip(*potential_cmpds))
            solution_dict[peak] = {
                "Compound": output[0], 
                "Protons": output[2],
                "Mass": output[1]         
            }
    
    end = time.time()
    print('Time ended:', time.ctime(end))
    print('Elapsed (seconds):', str((end-start)))

    return solution_dict