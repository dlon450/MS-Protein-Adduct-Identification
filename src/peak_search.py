import pandas as pd 
import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from isotope_pattern import *


def peak_find(bound_df: pd.DataFrame, peak_height: np.float):
    '''
    Find peaks above 'peak_height' (default 0.05) intensity
    '''
    peaks_idx, _ = find_peaks(bound_df['normalised_intensity'], height=peak_height)
    peaks = bound_df.loc[peaks_idx]

    return peaks, peaks_idx


def plot_peaks(bound_df: pd.DataFrame, peaks: pd.DataFrame):
    '''
    Plot MS and label peaks
    '''
    plt.plot(bound_df['m/z'], bound_df['normalised_intensity'])
    plt.plot(bound_df['m/z'][peaks], bound_df['normalised_intensity'][peaks], "x")
    plt.show()


def match_peaks(peaks, binding_df, bound_df, tolerance):
    '''
    Match peaks to theoretical list
    '''
    peaks_mass = peaks["m/z"].to_numpy()
    peaks_intensity = peaks["I"].to_numpy()
    binding_masses = binding_df["Mass"].to_numpy()
    
    # make new column in binding_df called Loss
    binding_df['Peak'] = None
    binding_df['Loss'] = -1
    binding_df['Best'] = False

    # intensities f()
    find_intensities = interpolate.interp1d(bound_df['m/z'].to_numpy(), bound_df['I'].to_numpy(), kind='linear')
    
    # find sites by matching
    binding_site_df = pd.DataFrame()
    all_binding_df = [pd.DataFrame()] * len(peaks_mass)

    for i, mass in enumerate(peaks_mass):

        lower_bound = np.round(mass - tolerance)
        upper_bound = np.round(mass + tolerance)

        search_mask = (binding_masses > lower_bound) & (binding_masses < upper_bound)
        if any(search_mask):
            # calculate_score(find_intensities, mass, search_mask, binding_df) # allocate DTW scores
            calculate_score_no_interpolation(mass, search_mask, binding_df, bound_df) #  allocate DTW scores with no interpolation

    binding_site_df = binding_df[binding_df['Loss'] != -1].sort_values(by='Mass')
    binding_site_df['Intensity'] = find_intensities(binding_site_df['Mass'].to_numpy())

    return binding_site_df