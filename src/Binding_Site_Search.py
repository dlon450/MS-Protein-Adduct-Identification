import pandas as pd 
import numpy as np

import config
from file_io import read
from normalise import normalise_spectrums
from peak_search import *
from theoretical_binding_list import create_binding_list

import time


def search(bound_file_path, unbound_file_path, compounds_file_path, tolerance=config.tolerance, \
    protein=config.protein, peak_height=config.peak_height, secondary='Platinum', plot_peak_graph=False):
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
    print(tolerance, protein, peak_height)
    bound_df, unbound_df, compounds = read(bound_file_path, unbound_file_path, compounds_file_path)

    # extract information about Ub and the other compounds
    other_compounds_df = compounds[~compounds.Name.isin([protein, secondary])] # extract all compounds but primary/secondary
    secondary_df = compounds[compounds.Name == secondary]
    protein_df = compounds[compounds.Name == protein] # extract Ub information

    # Perform peak detection
    # scale spectrums between 0 and 1
    unbound_df, bound_df = normalise_spectrums(unbound_df, bound_df)

    # find peaks 
    peaks, peaks_idx = peak_find(bound_df, np.float(peak_height))

    # create theoretical binding site list
    binding_df = create_binding_list(other_compounds_df, protein_df, secondary_df)

    # filter list using max value
    max_binding_mass = np.max(binding_df["Mass"]) # maximum possible mass from all theoretical binding sites
    peaks = peaks[peaks["m/z"] < max_binding_mass] # remove incorrect peaks
    
    binding_sites_df = match_peaks(peaks, binding_df, bound_df, np.float(tolerance))

    # format compound names
    binding_sites_df['Compound'] = [' + '.join(cmpd).translate(config.SUB) for cmpd in binding_sites_df['Compound']]

    if plot_peak_graph:
        plot_peaks(bound_df, peaks_idx)

    return binding_sites_df


if __name__ == "__main__":
    base_path = "Data/"
    fns = ["Ubi_O_1in100_broadband_000001", "Ubi_T_1in100_broadband_000001", "Ubiquitin_plusC_1in100_000001"]
    fns = ["Ubiquitin_plusC_1in100_000001"]
    unbound = base_path + "Deconvoluted Spectra/Ubi_1in100_broad band_000001.xlsx"
    compounds = base_path + "Compound Constraints/Compounds_CisOxTrans.xlsx"

    print("Searching...")
    for fn in fns:
        bound = base_path + "Deconvoluted Spectra/" + fn + ".xlsx"
        binding_sites = search(bound, unbound, compounds)
        # print(binding_sites[binding_sites.Best])
        print(binding_sites)
        # print(binding_sites.to_html())