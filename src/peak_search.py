from cProfile import label
import pandas as pd 
import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from isotope_pattern import *
import time


def peak_find(bound_df: pd.DataFrame, peak_height: float, min_dist_between_peaks=5.):
    '''
    Find distinct peaks above 'peak_height' (default 0.05) intensity 
    '''
    peaks_idx, _ = find_peaks(bound_df['normalised_intensity'], height=peak_height)
    peaks = bound_df.loc[peaks_idx]

    peak_masses = peaks["m/z"].to_numpy()
    n = len(peak_masses)
    peak_I = peaks["I"].to_numpy()
    keep = [False]*len(peak_masses)

    # remove peaks close together
    k = 0
    max_I = peak_I[k]
    for i in range(1, n):
        if peak_masses[i] - peak_masses[i-1] <= min_dist_between_peaks:
            current_I = peak_I[i]
            if current_I > max_I:
                max_I = current_I
                k = i
        else:
            keep[k] = True
            k = i
            max_I = peak_I[i]
    keep[k] = True

    return peak_masses[keep], peaks_idx, keep


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


def plot_peaks(bound_df: pd.DataFrame, peaks: pd.DataFrame, keep: np.array):
    '''
    Plot MS and label peaks
    '''
    plt.plot(bound_df['m/z'], bound_df['normalised_intensity'])
    plt.plot(bound_df['m/z'][peaks], bound_df['normalised_intensity'][peaks], "x", label='Identified peaks')
    plt.plot(bound_df['m/z'][peaks[keep]], bound_df['normalised_intensity'][peaks[keep]], "ko", markersize=3, label='Filtered peaks')
    plt.xlabel('m/z')
    plt.ylabel('Relative abundance')
    plt.legend()
    plt.show()