import pandas as pd 
import numpy as np

from os import listdir
import config


def read(bound_file_path, compounds_file_path, adducts_file_path):
    '''
    Read in excel files as pd.DataFrame objects
    '''
    bound_df = pd.read_excel(bound_file_path)
    compounds_df = pd.read_excel(compounds_file_path)
    adducts_df = pd.read_excel(adducts_file_path)
    adducts_df = adducts_df[adducts_df['Formula'] != 'H']
    adducts_df.columns = ['Compound/Fragment', 'Formula', 'Min', 'Max', 'Charge of compound/fragment']
    adducts_df['Compound/Fragment Type'] = 'Adducts'
    all_compounds = pd.concat([compounds_df, adducts_df], ignore_index=True, sort=False)

    return bound_df, all_compounds


def normalise(spectrum):
    X = spectrum["I"].to_numpy()
    spectrum["normalised_intensity"] = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    return spectrum