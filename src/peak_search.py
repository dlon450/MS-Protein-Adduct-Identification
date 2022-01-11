import pandas as pd 
import numpy as np
from scipy.signal import find_peaks
import config
import matplotlib.pyplot as plt

def peak_find(bound_df, peak_height):
    '''
    Find peaks above 'peak_height' (default 0.05) intensity
    '''
    peaks_idx, _ = find_peaks(bound_df['normalised_intensity'], height=peak_height)
    peaks = bound_df.loc[peaks_idx]

    return peaks, peaks_idx

def plot_peaks(bound_df, peaks):
    '''
    Plot MS and label peaks
    '''
    plt.plot(bound_df['m/z'], bound_df['normalised_intensity'])
    plt.plot(bound_df['m/z'][peaks], bound_df['normalised_intensity'][peaks], "x")
    plt.show()

def match_peaks(peaks, binding_df, tolerance):
    '''
    Match peaks to theoretical list
    '''
    peaks_mass = peaks["m/z"].to_numpy()
    peaks_intensity = peaks["I"].to_numpy()
    binding_masses = binding_df["Mass"].to_numpy()

    # find sites by matching
    binding_site_df = pd.DataFrame([])
    mask = np.zeros(binding_masses.size)

    for mass in peaks_mass:

        lower_bound = np.round(mass - tolerance)
        upper_bound = np.round(mass + tolerance)

        search_mask = (binding_masses > lower_bound) & (binding_masses < upper_bound)
        mask = mask + search_mask

    binding_site_df = binding_df[mask > 0.].sort_values(by='Mass').reset_index(drop=True)
    return binding_site_df
    
