# from cProfile import label
import pandas as pd 
import numpy as np
# from scipy import interpolate
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from isotope_pattern import *
import time


def peak_find(bound_df: pd.DataFrame, peak_height: float, min_dist_between_peaks=4., \
    calibrate="Automatic", protein_strs=['C378H629N105O118S1'], manual_calibration=0.):
    '''
    Find distinct peaks above 'peak_height' (default 0.01) intensity 
    '''
    peaks_idx, properties = find_peaks(bound_df['normalised_intensity'], height=peak_height)
    peaks = bound_df.loc[peaks_idx]

    peak_masses = peaks["m/z"].to_numpy()
    # n = len(peak_masses)
    peak_I = peaks["I"].to_numpy()
    # keep = [False]*len(peak_masses)

    args = np.argsort(-peak_I)
    keep = get_peaks(peak_masses[args])[np.argsort(args)] # inverse the argsort to get indices of peak_masses

    # remove peaks close together
    # k = 0
    # max_I = peak_I[k]
    # for i in range(1, n):
    #     if peak_masses[i] - peak_masses[i-1] <= min_dist_between_peaks:
    #         current_I = peak_I[i]
    #         if current_I > max_I:
    #             max_I = current_I
    #             k = i
    #     else:
    #         keep[k] = True
    #         k = i
    #         max_I = peak_I[i]
    # keep[k] = True

    if calibrate == "Automatic":
        shift = calibration_shift(peak_masses[keep], protein_formulas=protein_strs)
        bound_df['m/z'] += shift
        peak_masses += shift
    elif calibrate == "Manual":
        bound_df['m/z'] += manual_calibration
        peak_masses += manual_calibration

    return peak_masses[keep], peaks_idx, keep, dict(zip(peak_masses, properties['peak_heights']))


def get_peaks(masses):
    peaks = []
    masses_copy = masses.copy()
    while len(masses_copy) != 0:
        mass_of_greatest_peak = masses_copy[0]
        peaks.append(mass_of_greatest_peak)
        peaks_not_close_to_m = np.where(abs(masses_copy-mass_of_greatest_peak)>15.)
        masses_copy = masses_copy[peaks_not_close_to_m]
    return np.isin(masses, peaks)


def calibration_shift(peaks, protein_formulas):
    protein_peak = min([peak_isotope(formula) for formula in protein_formulas])
    return protein_peak - min(peaks, key=lambda x:abs(x-protein_peak))


def match_peaks(peak, binding_dict, bound_df, intensity, full=False, weight=10.):
    '''
    Match peaks to theoretical list
    '''
    start = time.time()
    # make new columns in binding_dict
    n_solns = len(binding_dict['Theoretical Peak Mass'])
    binding_dict['Experimental Peak'] = [peak]*n_solns
    binding_dict['Closeness of Fit (Loss)'] = [-1]*n_solns
    binding_dict['Closest Fit'] = [False]*n_solns

    # intensities f()
    # find_intensities = interpolate.interp1d(bound_df['m/z'].to_numpy(), bound_df['normalised_intensity'].to_numpy(), kind='linear')

    binding_site_record = calculate_score_no_interpolation(peak, binding_dict, bound_df, full, weight) #  allocate scores with no interpolation
    # binding_site_record['Intensity'] = find_intensities(binding_site_record['Theoretical Peak Mass'])
    binding_site_record['Intensity'] = intensity

    print('Elapsed (seconds):', str((time.time()-start)))
    if full:
        return pd.DataFrame(binding_site_record).sort_values(by=['Closeness of Fit (Loss)'])
    return binding_site_record


def plot_peaks(bound_df: pd.DataFrame, peaks: pd.DataFrame, keep: np.array, raw_MS=False):
    '''
    Plot MS and label peaks
    '''
    if raw_MS:
        plt.plot(bound_df['m/z'], bound_df['I'])
        plt.ylabel('Intensity')
    else:
        plt.plot(bound_df['m/z'], bound_df['normalised_intensity'])
        plt.plot(bound_df['m/z'][peaks], bound_df['normalised_intensity'][peaks], "x", label='Identified peaks')
        plt.plot(bound_df['m/z'][peaks[keep]], bound_df['normalised_intensity'][peaks[keep]], "ko", markersize=3, label='Filtered peaks')
        plt.ylabel('Relative abundance')
        plt.xlim([16500, 18500])

    plt.xlabel('m')
    plt.legend()
    plt.show()