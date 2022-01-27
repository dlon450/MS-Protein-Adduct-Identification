import pandas as pd
import numpy as np
from itertools import combinations  

def create_binding_list(other_compounds_df, protein_df, secondary_df):

    formulas = other_compounds_df["Formula"].to_numpy()
    masses = other_compounds_df["Mass"].to_numpy()

    base_compound = protein_df["Formula"].values[0]
    base_mass = protein_df["Mass"].values[0]

    if not secondary_df.empty:
        base_compound += secondary_df["Formula"].values[0] # Base Ub+Pt symbol 
        base_mass += secondary_df["Mass"].values[0] # Base Ub+Pt mass

    # binding list contains bonded compound and masses
    binding_df = pd.DataFrame([])
    number_of_compounds = other_compounds_df.shape[0]
    
    potential_bonded_cmpd = [""]*(2**number_of_compounds-1) # 2^n-1 ways from 1 to n
    potential_bonded_mass = [0]*(2**number_of_compounds-1)
    i = 0

    for k in range(number_of_compounds):

        # select k items from a list of n (where n=number of compounds, k=length of combinations)
        combs = combinations(np.arange(1,other_compounds_df.shape[0]+1), k+1)

        # create binding elements and total masses
        for bonded_compound_ids in combs:
            bonded_compound = [base_compound] + ['']*len(bonded_compound_ids)
            total_mass = base_mass

            for j, idx in enumerate(bonded_compound_ids):
                bonded_compound[j+1] = formulas[idx-1]
                total_mass += masses[idx-1]
            
            potential_bonded_cmpd[i] = bonded_compound
            potential_bonded_mass[i] = total_mass
            i += 1
    
    binding_df = pd.DataFrame({
        "Compound": potential_bonded_cmpd, 
        "Mass": potential_bonded_mass
        })

    return binding_df