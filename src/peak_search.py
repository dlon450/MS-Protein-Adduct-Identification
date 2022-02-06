import pandas as pd 
import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from isotope_pattern import *
import time


def peak_find(bound_df: pd.DataFrame, peak_height: float):
    '''
    Find peaks above 'peak_height' (default 0.05) intensity
    '''
    peaks_idx, _ = find_peaks(bound_df['normalised_intensity'], height=peak_height)
    peaks = bound_df.loc[peaks_idx]

    return peaks["m/z"].to_numpy(), peaks_idx


def plot_peaks(bound_df: pd.DataFrame, peaks: pd.DataFrame):
    '''
    Plot MS and label peaks
    '''
    plt.plot(bound_df['m/z'], bound_df['normalised_intensity'])
    plt.plot(bound_df['m/z'][peaks], bound_df['normalised_intensity'][peaks], "x")
    plt.show()


def match_peaks(peak, binding_dict, bound_df, full=False):
    '''
    Match peaks to theoretical list
    '''
    start = time.time()
    # make new columns in binding_dict
    n_solns = len(binding_dict['Mass'])
    binding_dict['Peak'] = [peak]*n_solns
    binding_dict['Loss'] = [-1]*n_solns
    binding_dict['Best'] = [False]*n_solns

    # intensities f()
    find_intensities = interpolate.interp1d(bound_df['m/z'].to_numpy(), bound_df['I'].to_numpy(), kind='linear')

    binding_site_record = calculate_score_no_interpolation(peak, binding_dict, bound_df, full) #  allocate DTW scores with no interpolation
    binding_site_record['Intensity'] = find_intensities(binding_site_record['Mass'])
    
    print('Elapsed (seconds):', str((time.time()-start)))
    if full:
        return pd.DataFrame(binding_site_record)
    return binding_site_record