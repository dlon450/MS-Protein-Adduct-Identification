# from wsgiref import validate
import pandas as pd 
# import numpy as np
import time
import config
from scipy import interpolate

from utils import read, normalise
from peak_search import peak_find, match_peaks, plot_peaks
from feasible_set import feasible_set_df
# from os import listdir
# from os.path import isfile, join


def search(bound_file_path, compounds_file_path, adducts_file_path, tolerance=config.tolerance, peak_height=config.peak_height,\
        multi_protein=config.multi_protein, min_primaries=config.min_primaries, max_primaries=config.max_primaries,\
            max_adducts=config.max_adducts, valence=config.valence, only_best=config.only_best, min_dist_between_peaks=4., \
                calibrate=config.calibrate, manual_calibration=0., plot_peak_graph=False, weight=10., return_peaks=False):
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
    print(f'\nCONFIGURATION: Tolerance={tolerance}, Peak_Height={peak_height}, Multi-protein={multi_protein}, Calibrate={calibrate}\n')
    full_data = only_best != 'on'
    multi_protein = multi_protein == 'on'
    bound_df, compounds = read(bound_file_path, compounds_file_path, adducts_file_path)

    bound_df = normalise(bound_df) # scale spectrums between 0 and 1
    protein_strs = compounds[compounds['Compound/Fragment Type'] == 'Protein']['Formula'].to_numpy()
    if len(protein_strs) > 1:
        multi_protein = True

    peaks, peaks_idx, keep, peak_intensities_dict = peak_find(bound_df, float(peak_height), \
        float(min_dist_between_peaks), calibrate, protein_strs, float(manual_calibration)) # find peaks 

    if plot_peak_graph:
        plot_peaks(bound_df, peaks_idx, keep)

    binding_dicts = feasible_set_df(compounds, peaks, float(tolerance), multi_protein, \
        int(min_primaries), int(max_primaries), int(max_adducts), int(valence)) # feasible set of integer combinations

    if return_peaks:
        return interpolate.interp1d(bound_df['m/z'].to_numpy(), bound_df['normalised_intensity'].to_numpy(), kind='linear'), peaks, list(binding_dicts.keys())

    # calculate objectives of solns
    print('\nFinding loss...')
    start = time.time()
    best_compounds = [{}]*len(binding_dicts)
    for i, (peak, binding_dict) in enumerate(binding_dicts.items()):
        print(f'Peak {round(peak, 2)} ------', end=' ')
        intensity = peak_intensities_dict[peak]
        binding_sites_record = match_peaks(peak, binding_dict, bound_df, intensity, full=full_data, weight=weight)
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

    return binding_sites_df[[
        'Species', 'Proton Offset', 'Intensity', 'Theoretical Peak Mass', 'Experimental Peak',
        'ppm', 'Closeness of Fit (Loss)', 'Closest Fit']]


if __name__ == "__main__":

    # compounds = "Data/Compound Constraints/Compounds_CisOxTrans_latest.xlsx"
    adducts = "Data/Compound Constraints/Standard_Adducts.xlsx"
    # bound = "Data/Deconvoluted Spectra/uc_medres_precal.xlsx"

    bound = r"Data\Input Data\RAPTA-C-20220421T034334Z-001\RAPTA-C\bound_spectrum_mb_rapc_precal.xlsx"
    compounds = r"Data\Input Data\RAPTA-C-20220421T034334Z-001\RAPTA-C\compounds_mb_rapc.xlsx"

    # plt.rcParams["figure.figsize"] = (18,3)
    binding_sites = search(bound, compounds, adducts)
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    print(binding_sites)
    # binding_sites.to_csv('test_output.csv', index=False)