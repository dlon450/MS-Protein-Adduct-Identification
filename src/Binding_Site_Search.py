import pandas as pd 
import numpy as np

import config
from utils import *
from peak_search import *
from feasible_set import feasible_set_df

import time


def search(bound_file_path, compounds_file_path, tolerance=config.tolerance, peak_height=config.peak_height,\
        multi_protein=config.multi_protein, min_primaries=config.min_primaries, max_primaries=config.max_primaries,\
            only_best=config.only_best, plot_peak_graph=False):
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
    print(f'\nCONFIGURATION: Tolerance={tolerance}, Peak_Height={peak_height}, Only_best={only_best}, Multi-protein={multi_protein}, n_primaries=[{min_primaries}, {max_primaries}]\n')
    full_data = only_best != 'on'
    multi_protein = multi_protein == 'on'
    bound_df, compounds = read(bound_file_path, compounds_file_path)

    bound_df = normalise(bound_df) # scale spectrums between 0 and 1
    peaks, peaks_idx, keep = peak_find(bound_df, float(peak_height)) # find peaks 
    binding_dicts = feasible_set_df(compounds, peaks, float(tolerance), multi_protein, int(min_primaries), int(max_primaries)) # feasible set of integer combinations

    # calculate objectives of solns
    print('\nFinding loss...')
    start = time.time()
    best_compounds = [{}]*len(binding_dicts)
    for i, (peak, binding_dict) in enumerate(binding_dicts.items()):
        print(f'Peak {round(peak, 2)} ------', end=' ')
        binding_sites_record = match_peaks(peak, binding_dict, bound_df, full=full_data)
        best_compounds[i] = binding_sites_record
    print('Elapsed (seconds):', str((time.time()-start)))

    if best_compounds == []:
        raise ValueError('No solutions found. Ensure constraints produce feasible solutions.')

    # format compound names
    if full_data:
        print('\nCombining DataFrame...')
        binding_sites_df = pd.concat(best_compounds)
    else:
        binding_sites_df = pd.DataFrame(best_compounds)
    binding_sites_df['Compound'] = [' + '.join(cmpd) for cmpd in binding_sites_df['Compound']]

    if plot_peak_graph:
        plot_peaks(bound_df, peaks_idx, keep)
    return binding_sites_df


if __name__ == "__main__":

    compounds = "Data/Compound Constraints/Compounds_CisOxTrans_nlp.xlsx"
    bound = "Data/Deconvoluted Spectra/Ubiquitin_plusC_1in100_000001.xlsx"

    # compounds = "Data/Other Spectra/compounds_rc.xlsx"
    # bound = "Data/Other Spectra/bound_spectrum_rc.xlsx"

    binding_sites = search(bound, compounds)
    print(binding_sites.head(10))