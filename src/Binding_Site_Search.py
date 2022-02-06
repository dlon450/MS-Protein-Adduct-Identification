import pandas as pd 
import numpy as np

import config
from utils import *
from peak_search import *
from theoretical_binding_list import feasible_set_df

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
    print(f'\n-------- CONFIGURATION: Tolerance={tolerance}, Peak_Height={peak_height} --------\n')
    bound_df, unbound_df, compounds = read(bound_file_path, unbound_file_path, compounds_file_path)

    unbound_df, bound_df = normalise_spectrums(unbound_df, bound_df) # scale spectrums between 0 and 1
    peaks, peaks_idx = peak_find(bound_df, float(peak_height)) # find peaks 
    binding_dicts = feasible_set_df(compounds, peaks, float(tolerance)) # feasible set of integer combinations

    # calculate objectives of solns
    print('\nFinding loss...')
    start = time.time()
    full_data = True
    best_compounds = [{}]*len(binding_dicts)
    for i, (peak, binding_dict) in enumerate(binding_dicts.items()):
        print(f'Peak {round(peak, 2)} ------', end=' ')
        binding_sites_record = match_peaks(peak, binding_dict, bound_df, full=full_data)
        best_compounds[i] = binding_sites_record
    print('Elapsed (seconds):', str((time.time()-start)))

    # format compound names
    if full_data:
        print('\nCombining DataFrame...')
        binding_sites_df = pd.concat(best_compounds)
    else:
        binding_sites_df = pd.DataFrame(best_compounds)
    binding_sites_df['Compound'] = [' + '.join(cmpd).translate(config.SUB) for cmpd in binding_sites_df['Compound']]

    if plot_peak_graph:
        plot_peaks(bound_df, peaks_idx)
    return binding_sites_df


if __name__ == "__main__":
    base_path = "Data/"
    fns = ["Ubi_O_1in100_broadband_000001", "Ubi_T_1in100_broadband_000001", "Ubiquitin_plusC_1in100_000001"]
    fns = ["Ubiquitin_plusC_1in100_000001"]
    unbound = base_path + "Deconvoluted Spectra/Ubi_1in100_broad band_000001.xlsx"
    compounds = base_path + "Compound Constraints/Compounds_CisOxTrans_nlp.xlsx"

    for fn in fns:
        bound = base_path + "Deconvoluted Spectra/" + fn + ".xlsx"
        binding_sites = search(bound, unbound, compounds)
        print(binding_sites)