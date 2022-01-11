import pandas as pd
import numpy as np
from itertools import combinations  

import config


def create_binding_list(other_compounds_df, protein_df):

    formulas = other_compounds_df["Formula"].to_numpy()
    masses = other_compounds_df["Mass"].to_numpy()
    base_compound = protein_df["Formula"].values[0].translate(config.SUB) # Base Ub symbol
    base_mass = protein_df["Mass"].values[0] # Base Ub mass

    # binding list contains bonded compound and masses
    binding_df = pd.DataFrame([])
    number_of_compounds = other_compounds_df.shape[0]
    
    potential_bonded_cmpd = [""]*(2**number_of_compounds-1) # 2^n-1 ways from 1 to n
    potential_bonded_mass = [0]*(2**number_of_compounds-1)
    i = 0

    for k in range(number_of_compounds):

        # select k items from a list of n (where n=number of compounds, k=length of combinations)
        combs = combinations(other_compounds_df.ID, k+1)

        # create binding elements and total masses
        for bonded_compound_ids in combs:
            bonded_compound = base_compound
            total_mass = base_mass

            for idx in bonded_compound_ids:
                bonded_compound += " + " + formulas[idx-1].translate(config.SUB)
                total_mass += masses[idx-1]
            
            potential_bonded_cmpd[i] = bonded_compound
            potential_bonded_mass[i] = total_mass
            i += 1
    
    binding_df = pd.DataFrame({"Compound": potential_bonded_cmpd, "Mass": potential_bonded_mass})

    return binding_df