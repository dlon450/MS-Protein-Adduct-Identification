import pandas as pd 
import numpy as np

import config
from file_io import read
from normalise import normalise_spectrums
from peak_search import peak_find, match_peaks
from theoretical_binding_list import create_binding_list

import time

def search(bound_file_path, unbound_file_path, compounds_file_path):
    '''
    Search for appropriate binding sites

    Parameters
    ----------
    bound_file_path: str
    unbound_file_path: str
    compounds_file_path: str

    Returns
    -------
    binding_sites_df: pd.DataFrame
    '''
    bound_df, unbound_df, compounds = read(bound_file_path, unbound_file_path, compounds_file_path)

    # extract information about Ub and the other compounds
    other_compounds_df = compounds[compounds.Name != config.protein] # extract all compounds but Ub
    protein_df = compounds[compounds.Name == config.protein] # extract Ub information

    # Perform peak detection
    # scale spectrums between 0 and 1
    unbound_df, bound_df = normalise_spectrums(unbound_df, bound_df)

    # find peaks 
    peaks = peak_find(bound_df)

    # create theoretical binding site list
    binding_df = create_binding_list(other_compounds_df, protein_df)

    # filter list using max value
    max_binding_mass = np.max(binding_df["Mass"]) # maximum possible mass from all theoretical binding sites
    peaks = peaks[peaks["m/z"] < max_binding_mass] # remove incorrect peaks

    # match peaks to theoretical list
    binding_sites_df = match_peaks(peaks, binding_df)

    return binding_sites_df

if __name__ == "__main__":
    base_path = "Data/"
    bound = base_path + "Deconvoluted Spectra/Ubiquitin_plusC_1in100_000001.xlsx"
    unbound = base_path + "Deconvoluted Spectra/Ubi_1in100_broad band_000001.xlsx"
    compounds = base_path + "Compound Constraints/Compounds_CisOxTrans.xlsx"

    binding_sites = search(bound, unbound, compounds)