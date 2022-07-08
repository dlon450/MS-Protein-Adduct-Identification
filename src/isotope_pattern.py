# from multiprocessing.sharedctypes import Value
from pyopenms import FineIsotopePatternGenerator, EmpiricalFormula, CoarseIsotopePatternGenerator
import numpy as np
import similaritymeasures
from matplotlib import pyplot as plt
# from icecream import ic
# from scipy.spatial.distance import euclidean
# from fastdtw import fastdtw, dtw
# import time

PROTON_MASS = 1.0078250319


def peak_isotope(formula: object, accuracy=1e-3, peak=None) -> float:
    '''
    Use isotope pattern to find peak isotope
    '''
    masses, rel_abundance = find_isotope_pattern(formula, accuracy, peak)
    return masses[np.argmax(rel_abundance)]


def peak_isotope_old(formula: object) -> float:
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


def find_isotope_pattern(formula_str: str, accuracy=1e-3, peak=None, interval=5.):
    '''
    Return theoretical distribution of intensities
    '''
    seq_formula = EmpiricalFormula(formula_str)
    isotopes = seq_formula.getIsotopeDistribution( FineIsotopePatternGenerator(accuracy) )
    isotopes_container = isotopes.getContainer()
    masses_dict = {}
    abundances_dict = {}

    if peak is not None:
        lower, upper = binarySearchInterval(isotopes_container, np.floor(peak - interval), np.ceil(peak + interval))
        isotopes_container = isotopes_container[lower:upper]

    for iso in isotopes_container:
        mz = iso.getMZ()
        group = int(mz)
        abundance = iso.getIntensity()

        if group in masses_dict:
            masses_dict[group] += mz*abundance
            abundances_dict[group] += abundance
        else:
            masses_dict[group] = mz*abundance
            abundances_dict[group] = abundance

    rel_abundances = np.array(list(abundances_dict.values()))
    masses = np.array(list(masses_dict.values())) / rel_abundances
    return masses, rel_abundances


def find_isotope_pattern_coarse(formula_str: str, peak=None, interval=5.):
    '''
    Return coarse theoretical distribution of intensities
    '''
    seq_formula = EmpiricalFormula(formula_str)
    isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator() )
    isotopes_container = isotopes.getContainer()

    if peak is not None:
        lower, upper = binarySearchInterval(isotopes_container, np.floor(peak - interval), np.ceil(peak + interval))
        isotopes_container = isotopes_container[lower:upper]
    
    o = np.transpose([[iso.getMZ(), iso.getIntensity()] for iso in isotopes_container])
    return o[0], o[1]


def find_isotope_pattern_old(formula_str: str):
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
        

def calculate_score_no_interpolation(peak_mass, binding_dict, bound_df, full=False, weight=10.):
    '''
    Return the compounds that best match the experimental distribution
    '''
    # find sites by matching
    binding_site_record = dict.fromkeys(binding_dict.keys())
    binding_dict['Theoretical Peak Mass'] = list(binding_dict['Theoretical Peak Mass'])
    
    exp_masses = bound_df['m/z'].to_numpy()
    p = exp_masses.searchsorted(peak_mass)
    start = exp_masses.searchsorted(peak_mass - 6.)
    stop = exp_masses.searchsorted(peak_mass + 6.)
    exp_masses = exp_masses[start:stop]
    
    exp_abundance = bound_df['I'].to_numpy()
    exp_abundance = (exp_abundance / exp_abundance[p])[start:stop]
    keep = np.where(np.logical_and(exp_abundance > 0.15, exp_abundance < 1.001))
    # keep = np.where(exp_abundance > 0.)
    
    peak_compounds = binding_dict['Species']
    proton_offset = binding_dict['Proton Offset']

    best_idx = 0
    best_distance = np.inf

    for i, compound_list in enumerate(peak_compounds):
        
        # plt.plot(exp_masses, exp_abundance * weight, label='Experimental')
        distance, isotope_peak_mass = objective_func(''.join(compound_list), peak_mass, \
            exp_masses[keep], exp_abundance[keep], proton_offset[i], weight)
        binding_dict['Closeness of Fit (Loss)'][i] = distance

        theoretical_peak_mass = isotope_peak_mass
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


def objective_func(formula, peak_mass, m, y, proton_offset, weight=1.):
    '''
    Returns the DTW loss
    '''
    if formula == '':
        return 99999
    
    masses, x = find_isotope_pattern_coarse(f'{formula}H{-int(proton_offset)}', peak=peak_mass)
    max_idx_x = maxInBitonic(x, 0, len(x) - 1)
    isotope_peak = masses[max_idx_x]
    # masses = np.array(masses)
    # start = masses.searchsorted(isotope_peak - 6.)
    # stop = masses.searchsorted(isotope_peak + 6.)

    try:
        # distance, _ = fastdtw(list(zip(masses, np.array(x) / x[max_idx_x] * weight)), list(zip(m, y * weight)), dist=euclidean)
        # distance = similaritymeasures.area_between_two_curves(np.array(list(zip(m, y * weight))), np.array(list(zip(masses, np.array(x) / x[max_idx_x] * weight))))
        # distance = similaritymeasures.frechet_dist(list(zip(masses, np.array(x) / x[max_idx_x] * weight)), list(zip(m, y * weight)))
        # distance = similaritymeasures.pcm(np.array(list(zip(masses, np.array(x) / x[max_idx_x] * weight))), np.array(list(zip(m, y * weight))))
        distance = similaritymeasures.dtw(np.array(list(zip(m, y * weight))), np.array(list(zip(masses, np.array(x) / x[max_idx_x] * weight))), metric='euclidean')[0]
    except ValueError:
        return 99999., isotope_peak

    # plt.plot(m, y * weight, label='Experimental')
    # plt.plot(masses, np.array(x) * weight / x[max_idx_x], label='Theoretical')
    # plt.legend()
    # plt.xlim([peak_mass-10., peak_mass+10.])
    # plt.savefig(f'isotope_matches/{peak_mass}_{formula}.png')
    # plt.clf()

    return distance / weight, isotope_peak


def maxInBitonic(arr, l, r) :
    '''
    Binary search in Bitonic array https://www.geeksforgeeks.org/find-the-maximum-element-in-an-array-which-is-first-increasing-and-then-decreasing/
    (Returns index of maximum)
    '''
    while (l <= r) :
        try:
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
        except IndexError:
            return m
    return -1


def binarySearchInterval(isotopes, lbound, ubound):
    '''
    Find isotopes with masses in given interval
    '''
    low = 0
    high = len(isotopes) - 1

    lower_idx = binarySearch(isotopes, low, high, lbound)
    upper_idx = binarySearch(isotopes, lower_idx, high, ubound) - 1
    return lower_idx, upper_idx


def binarySearch(isotopes, low, high, X):
    '''
    Search for index position of mass X in isotopes
    '''
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        mass = isotopes[mid].getMZ()
        if mass < X:
            low = mid + 1
        elif mass > X:
            high = mid - 1
        else:
            return mid
    return high + 1


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


def find_species_additive_mass(species=["C769H1212N210O218S2", "Pt", "NH3", "H2O", "Cl"], counts=[1, 2, 4, 0, 2]):
    masses = [peak_isotope(s)*c for s, c in zip(species, counts)]
    return sum(masses)


if __name__ == "__main__":

    # formula = 'C378H629N105O118S1'                  # Bruker: 8564.630429
    # print(peak_isotope(formula, accuracy=1e-3))     # 8564.630447593447
    formula = 'C769H1212N210O218S2'
    # find_nominal_masses(formula)
    print(peak_isotope(formula))
    # isotopes = find_isotope_pattern(formula, accuracy=1e-3)
    # plt.rcParams["figure.figsize"] = (15,10)
    # plotIsotopeDistribution(isotopes, formula, save=False)