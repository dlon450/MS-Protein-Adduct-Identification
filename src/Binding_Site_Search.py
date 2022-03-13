import pandas as pd 
import numpy as np
import time
import config

from utils import *
from peak_search import *
from feasible_set import feasible_set_df
from os import listdir
from os.path import isfile, join


def search(bound_file_path, compounds_file_path, adducts_file_path, tolerance=config.tolerance, peak_height=config.peak_height,\
        multi_protein=config.multi_protein, min_primaries=config.min_primaries, max_primaries=config.max_primaries,\
            max_adducts=config.max_adducts, valence=config.valence, only_best=config.only_best, min_dist_between_peaks=5., \
                calibrate=config.calibrate, plot_peak_graph=False):
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
    calibrate = calibrate == 'on'
    bound_df, compounds = read(bound_file_path, compounds_file_path, adducts_file_path)

    bound_df = normalise(bound_df) # scale spectrums between 0 and 1
    protein_str = compounds[compounds['Compound/Fragment Type'] == 'Protein']['Formula'].to_numpy()[0]
    peaks, peaks_idx, keep = peak_find(bound_df, float(peak_height), float(min_dist_between_peaks), calibrate, protein_str) # find peaks 
    binding_dicts = feasible_set_df(compounds, peaks, float(tolerance), multi_protein, \
        int(min_primaries), int(max_primaries), int(max_adducts), int(valence)) # feasible set of integer combinations

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
    binding_sites_df['Species'] = [' + '.join(cmpd) for cmpd in binding_sites_df['Species']]

    if plot_peak_graph:
        plot_peaks(bound_df, peaks_idx, keep)
    return binding_sites_df[[
        'Species', 'Proton Offset', 'Intensity', 'Theoretical Peak Mass', 'Experimental Peak',
        'ppm', 'Closeness of Fit (Loss)', 'Closest Fit']]


def search_all(dirpath, compounds_file, adducts_file):

    calibrated = '_softwareCal' if config.calibrate == 'on' else '_noSWCal'
    bound_files = [f for f in listdir(dirpath) if isfile(join(dirpath, f))]
    for bound_fn in bound_files:
        binding_sites = search(join(dirpath, bound_fn), compounds_file, adducts_file)
        output_fn = bound_fn[:bound_fn.rfind('.')] + f'{calibrated}.csv'
        binding_sites.to_csv(join(dirpath, output_fn), index=False)


if __name__ == "__main__":

    compounds = "Data/Compound Constraints/Compounds_CisOxTrans_latest.xlsx"
    adducts = "Data/Compound Constraints/Standard_Adducts.xlsx"
    bound = "Data/Deconvoluted Spectra/Ubiquitin_plusC_1in100_000001.xlsx"
    
    # binding_sites = search(bound, compounds, adducts)
    # pd.set_option("display.max_rows", None)
    # print(binding_sites)
    # binding_sites.to_csv('Data/cristian_data_1.xlsx', index=False)

    search_all('Data/Deconvoluted Spectra', compounds, adducts)