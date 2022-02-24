import pandas as pd 
import numpy as np

from os import listdir
import config


def read(bound_file_path, compounds_file_path):
    '''
    Read in excel files as pd.DataFrame objects
    '''
    bound_df = pd.read_excel(bound_file_path)
    compounds_df = pd.read_excel(compounds_file_path)

    return bound_df, compounds_df


def normalise(spectrum):
    X = spectrum["I"].to_numpy()
    spectrum["normalised_intensity"] = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    return spectrum