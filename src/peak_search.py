import pandas as pd 
import numpy as np
from scipy.signal import find_peaks
import config

def peak_find(bound_df):
    '''
    Find peaks above 'peak_height' (default 0.05) intensity
    '''
    peaks_idx, _ = find_peaks(bound_df['normalised_intensity'], height=config.peak_height)
    peaks = bound_df.loc[peaks_idx]

    return peaks

def match_peaks(peaks, binding_df):
    # find sites by matching
    binding_site_df = pd.DataFrame([])

    for idx in peaks.index:
        mass = peaks.loc[idx, "m/z"]
        intensity = peaks.loc[idx, "I"]

        lower_bound = np.round(mass - config.tolerance)
        upper_bound = np.round(mass + config.tolerance)

        search_mask = (binding_df["Mass"] > lower_bound) & (binding_df["Mass"] < upper_bound)
        number_of_searches = np.sum(search_mask)
        
        if number_of_searches > 0:
            binding_sites = binding_df[search_mask]
            binding_site_df = binding_site_df.append(binding_sites, ignore_index=True)

    return binding_site_df
    
