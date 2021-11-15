import pandas as pd
import numpy as np
from itertools import combinations  

import config


def create_binding_list(other_compounds_df, protein_df):

    #binding list contains bonded compound and masses
    binding_df = pd.DataFrame([])
    number_of_compounds = other_compounds_df.shape[0]

    for k in range(1, number_of_compounds+1):

        #how many ways can you select k items from a list of n (where n=number of compounds, k=length of combinations)
        #use ID's as from ID's we can extract the masses and compound names
        combs = combinations(other_compounds_df.ID,k)

        #create binding elements and total masses
        for bonded_compound_ids in combs:
            bonded_compound = protein_df["Formula"].values[0].translate(config.SUB) #Base Ub symbol
            total_mass = protein_df["Mass"].values[0] #Base Ub mass

            for idx in bonded_compound_ids:
                bonded_compound += " + " + other_compounds_df.loc[other_compounds_df.ID == idx, "Formula"].values[0].translate(config.SUB)
                total_mass += other_compounds_df.loc[other_compounds_df.ID == idx, "Mass"].values[0]
            
            binding_df = binding_df.append(pd.DataFrame({"Compound":bonded_compound, "Mass":total_mass}, index=range(1)))

    binding_df.index = range(binding_df.shape[0])

    return binding_df