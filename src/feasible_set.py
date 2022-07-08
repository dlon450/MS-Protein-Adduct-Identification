# import pandas as pd
import numpy as np
from isotope_pattern import peak_isotope
import time
from constraint_optimisation import feasible_set_search

PROTON_MASS = 1.0078250319


def feasible_set_df(compounds, peaks_mass, tolerance, multi_protein=False, min_primaries=None, max_primaries=None, \
    max_adducts=None, valence=None, precision=5):
    '''
    Returns feasible set of integer solutions 
    '''
    factor = 10.**precision
    solution_dict = {}

    formulas = compounds["Formula"].to_numpy()
    charges = compounds["Charge of compound/fragment"].to_numpy()
    masses = np.vectorize(peak_isotope)(formulas) - np.dot(PROTON_MASS, charges)

    masses = (masses*factor).astype(int)

    masses_dict = dict(zip(formulas, masses))
    charges_dict = dict(zip(formulas, charges))
    max_amount_dict = dict(zip(formulas, compounds["Max"].to_numpy()))
    min_amount_dict = dict(zip(formulas, compounds["Min"].to_numpy()))
    adducts_dict = dict(zip(formulas, compounds['Compound/Fragment Type'] == 'Adducts'))
    metal_idx = compounds.index[compounds['Compound/Fragment Type'] == 'Metal'].tolist()[0]
    max_per_metal = dict(zip(formulas, compounds['Maximum number of the corresponding ligand per metal'].to_numpy()))

    primaries_dict = None
    if multi_protein:
        primaries_dict = dict(zip(formulas, compounds['Compound/Fragment Type'] == 'Protein'))

    print('Searching for solutions...')
    start = time.time()
    print('Time started:', time.ctime(start))

    for peak in peaks_mass:
        print(f'Peak {round(peak, 2)} ------ ', end='')
        solutions, count = feasible_set_search(formulas, masses_dict, max_amount_dict, min_amount_dict,\
            peak_mass=int(peak*factor), tolerance=int(tolerance*factor), multi_protein=multi_protein, \
                primaries=primaries_dict, min_primaries=min_primaries, max_primaries=max_primaries, \
                    metal_idx=metal_idx, max_per_metal=max_per_metal, adducts=adducts_dict, \
                        max_adducts=max_adducts, valence=valence)

        if count != 0:
            potential_cmpds = [[[
                c[0]*c[1] for c in soln if c[1] != 0],
                sum([masses_dict[c[0]]*c[1] for c in soln if c[1] != 0])/factor,
                sum([charges_dict[c[0]]*c[1] for c in soln if c[1] != 0])] 
                for soln in solutions]
            output = list(zip(*potential_cmpds))
            solution_dict[peak] = {
                "Species": output[0], 
                "Proton Offset": output[2],
                "Theoretical Peak Mass": output[1],
                "ppm": abs(output[1] - peak) / output[1] * 1000000
            }
    
    end = time.time()
    print('Time ended:', time.ctime(end))
    print('Elapsed (seconds):', str((end-start)))

    return solution_dict