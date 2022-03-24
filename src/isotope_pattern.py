from pyopenms import *
import numpy as np
from matplotlib import pyplot as plt
from icecream import ic
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw, dtw
import time

PROTON_MASS = 1.0078250319


def peak_isotope(formula: object) -> float:
    '''
    Use isotope pattern to find average mass
    '''
    iso_dict = find_nominal_masses(formula)
    best_abundance = 0.
    best_mass = 0.
    for _, value in iso_dict.items():
        abundance = sum([i for m,i in value])
        if abundance > best_abundance:
            best_mass = sum([m*i for m,i in value]) / abundance
            best_abundance = abundance

    return best_mass


def find_nominal_masses(formula_str: str):
    '''
    Return dictionary of isotopes at all possible nominal masses 
    '''
    seq_formula = EmpiricalFormula(formula_str)
    isotopes = seq_formula.getIsotopeDistribution( FineIsotopePatternGenerator(1e-3) )
    iso_dict = {}
    for iso in isotopes.getContainer():
        mz = iso.getMZ()
        group = int(mz)
        iso_dict[group] = iso_dict.get(group, []) + [(mz, iso.getIntensity())]
    return iso_dict


def find_isotope_pattern(formula_str: str):
    '''
    Return theoretical distribution of intensities
    '''

    iso_dict = find_nominal_masses(formula_str)
    isotopes = iso_dict.values()
    n = len(isotopes)
    masses = [0.]*n
    rel_abundance = [0.]*n

    for i, iso in enumerate(isotopes):
        for mass, abundance in iso:
            masses[i] += mass*abundance
            rel_abundance[i] += abundance
        masses[i] = masses[i] / rel_abundance[i]

    return masses, rel_abundance
        

def calculate_score_no_interpolation(peak_mass, binding_dict, bound_df, full=False):
    '''
    Return the compounds that best match the experimental distribution
    '''
    # find sites by matching
    binding_site_record = dict.fromkeys(binding_dict.keys())
    binding_dict['Theoretical Peak Mass'] = list(binding_dict['Theoretical Peak Mass'])
    
    exp_masses = bound_df['m/z'].to_numpy()
    exp_abundance = bound_df['I'].to_numpy()
    exp_abundance = exp_abundance / exp_abundance[exp_masses.searchsorted(peak_mass)]
    keep = np.where(np.logical_and(exp_abundance > 0.1, exp_abundance < 1.001))
    peak_compounds = binding_dict['Species']
    proton_offset = binding_dict['Proton Offset']

    best_idx = 0
    best_distance = np.inf

    for i, compound_list in enumerate(peak_compounds):
        
        distance, isotope_peak_mass = objective_func(''.join(compound_list), \
            exp_masses[keep], exp_abundance[keep], peak_mass, proton_offset[i])
        binding_dict['Closeness of Fit (Loss)'][i] = distance

        theoretical_peak_mass = isotope_peak_mass # - proton_offset[i] * PROTON_MASS
        binding_dict['Theoretical Peak Mass'][i] = theoretical_peak_mass
        binding_dict['ppm'][i] = abs(theoretical_peak_mass - peak_mass) / theoretical_peak_mass * 1000000 

        if distance < best_distance:
            best_distance = distance
            best_idx = i
    
    # plotWarpDTW(x, y, warp_path)
    binding_dict['Closest Fit'][best_idx]= True
    for key in binding_dict.keys():
        binding_site_record[key] = binding_dict[key][best_idx]
    
    if full:
        return binding_dict
    return binding_site_record


def objective_func(formula, exp_masses, exp_abundance, peak_mass, proton_offset, weight=100.):
    '''
    Returns the DTW loss
    '''
    if formula == '':
        return 99999

    masses, x = find_isotope_pattern(f'{formula}H-{proton_offset}')
    max_idx_x = maxInBitonic(x, 0, len(x) - 1)
    isotope_peak = masses[max_idx_x]
    difference = peak_mass - isotope_peak
    
    start = exp_masses.searchsorted(difference + masses[0])
    stop = exp_masses.searchsorted(difference + masses[-1], side='right')

    y = exp_abundance[start:stop]
    m = exp_masses[start:stop]

    distance, _ = fastdtw(list(zip(masses, np.array(x) / x[max_idx_x] * weight)), list(zip(m, y * weight)), dist=euclidean)
    plt.plot(m, y * weight)
    plt.plot(masses, np.array(x) * weight / x[max_idx_x])
    plt.clf()
    return distance / weight, isotope_peak


def maxInBitonic(arr, l, r) :
    '''
    Binary search in Bitonic array https://www.geeksforgeeks.org/find-the-maximum-element-in-an-array-which-is-first-increasing-and-then-decreasing/
    (Returns index of maximum)
    '''
    while (l <= r) :
        m = int(l + (r - l) / 2)
        if ((r == l + 1) and arr[l] >= arr[r]):
            return l
        if ((r == l + 1) and arr[l] < arr[r]):
            return r
        if (arr[m] > arr[m + 1] and arr[m] > arr[m - 1]):
            return m
        if (arr[m] > arr[m + 1] and arr[m] < arr[m - 1]) :
            r = m - 1
        else :
            l = m + 1
    return -1


def plotWarpDTW(x1, x2, warp_path):
    '''
    Plot curves https://towardsdatascience.com/an-illustrative-introduction-to-dynamic-time-warping-36aa98513b98
    '''
    fig, ax = plt.subplots(figsize=(16, 12))

    # Remove the border and axes ticks
    fig.patch.set_visible(False)
    ax.axis('off')

    for [map_x, map_y] in warp_path:
        ax.plot([map_x, map_y], [x1[map_x], x2[map_y]], '-k')

    ax.plot(x1, color='blue', marker='o', markersize=10, linewidth=5)
    ax.plot(x2, color='red', marker='o', markersize=10, linewidth=5)
    ax.tick_params(axis="both", which="major", labelsize=18)

    fig.show()


def plotIsotopeDistribution(isotopes, formula, save=False):
    '''
    Plot Isotope Distribution (https://pyopenms.readthedocs.io/en/latest/aasequences.html)
    '''
    xs = isotopes[0]
    ys = (np.array(isotopes[1]) - min(isotopes[1]))/(max(isotopes[1]) - min(isotopes[1]))
    
    plt.bar(xs, ys, width=0.05, snap=False)
    plt.xlabel("Atomic mass (u)")
    plt.ylabel("Relative abundance (%)")
    plt.ylim([0, 1.2])
    plt.title(f'{formula} Isotope Distribution')

    for x,y in zip(xs,ys):
        if y >= 0.01:
            label = "{:.4f}".format(y)
            plt.annotate(label, (x,y), textcoords="offset points", xytext=(0,10), ha='center')

    if save:
        plt.savefig(f'Data/{formula}_isotope_pattern.png')
    plt.show()


def missing_elements(fn=r'C:\Users\longd\Downloads\List of elements.csv'):
    '''
    Find missing elements
    '''
    import pandas as pd
    df = pd.read_csv(fn)
    element_symbols = df.Symbol.to_numpy()
    missing = []
    for s in element_symbols:
        try:
            EmpiricalFormula(s)
        except RuntimeError:
            missing.append(s)
    df[df.Symbol.isin(missing)][['Atomic number (Z)', 'Symbol','Name']].to_csv('missing_elements.csv', index_col=False)


if __name__ == "__main__":

    formula = 'C378H629N105O118S1PtNH3NH3Cl'
    isotopes = find_isotope_pattern(formula)
    plt.rcParams["figure.figsize"] = (20,10)
    plotIsotopeDistribution(isotopes, formula, True)