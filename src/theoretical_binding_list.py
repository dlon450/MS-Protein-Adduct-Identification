import pandas as pd
import numpy as np
from itertools import combinations  
from isotope_pattern import average_mass
import time
from constraint_optimisation import feasible_set_search

PROTON_MASS = 1.007825


def feasible_set_df(compounds, peaks_mass, tolerance, precision=5):
    '''
    Returns feasible set of integer solutions 
    '''
    factor = 10.**precision
    solution_dict = {}

    formulas = compounds["Formula"].to_numpy()
    charges = compounds["Charge"].to_numpy()
    max_amount = compounds["Max"].to_numpy()
    min_amount = compounds["Min"].to_numpy()
    masses = np.vectorize(average_mass)(formulas) - np.dot(PROTON_MASS, charges)

    masses = (masses*factor).astype(int)

    masses_dict = dict(zip(formulas, masses))
    max_amount_dict = dict(zip(formulas, max_amount))
    min_amount_dict = dict(zip(formulas, min_amount))
    charges_dict = dict(zip(formulas, charges))

    print('Searching for solutions...')
    start = time.time()
    print('Time started:', time.ctime(start))

    for peak in peaks_mass:
        print(f'Peak {round(peak, 2)} ------ ', end='')
        solutions, count = feasible_set_search(formulas, masses_dict, max_amount_dict, min_amount_dict, \
            peak_mass=int(peak*factor), tolerance=int(tolerance*factor))

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


def create_binding_list(other_compounds_df, protein_df, secondary_df):

    # species_mass = sum_of_Mass - Charge*PROTON_MASS

    formulas = other_compounds_df["Formula"].to_numpy()
    charges = other_compounds_df["Charge"].to_numpy()
    masses = np.vectorize(average_mass)(formulas)

    base_compound = protein_df["Formula"].values[0]
    base_charge = protein_df["Charge"].values[0]
    base_mass = average_mass(base_compound) - base_charge * PROTON_MASS

    if not secondary_df.empty:
        secondary_compound = secondary_df["Formula"].values[0] # Base Ub+Pt symbol 
        secondary_charge = secondary_df["Charge"].values[0] 
        base_mass += average_mass(secondary_compound) - secondary_charge * PROTON_MASS

    # binding list contains bonded compound and masses
    number_of_compounds = other_compounds_df.shape[0]
    
    potential_bonded_cmpd = [""]*(2**number_of_compounds-1) # 2^n-1 ways from 1 to n
    potential_bonded_mass = [0]*(2**number_of_compounds-1)
    potential_charge = [base_charge + secondary_charge]*(2**number_of_compounds-1)
    i = 0

    for k in range(number_of_compounds):

        # select k items from a list of n (where n=number of compounds, k=length of combinations)
        combs = combinations(np.arange(1,other_compounds_df.shape[0]+1), k+1)

        # create binding elements and total masses
        for bonded_compound_ids in combs:
            bonded_compound = ['']*len(bonded_compound_ids)
            total_mass = base_mass

            for j, idx in enumerate(bonded_compound_ids):
                bonded_compound[j] = formulas[idx-1] # multiply by min 
                total_mass += masses[idx-1] - charges[idx-1] * PROTON_MASS
                potential_charge[i] += charges[idx-1]
            
            potential_bonded_cmpd[i] = [base_compound, secondary_compound] + bonded_compound
            potential_bonded_mass[i] = total_mass
            i += 1

        print(str(k+1) + ' out of ' + str(number_of_compounds))
    
    return pd.DataFrame({
        "Compound": potential_bonded_cmpd, 
        "Protons": potential_charge,
        "Mass": potential_bonded_mass
        })