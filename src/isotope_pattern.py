from pyopenms import *
import numpy as np
from matplotlib import pyplot as plt
from icecream import ic
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw, dtw
import time


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


def find_isotope_pattern(formula_str: str, generator_num=16, paired=False, plot_distribution=False):
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

    exp_masses = bound_df['m/z'].to_numpy()
    exp_abundance = bound_df['I'].to_numpy()

    number_of_points = 16
    peak_compounds = binding_dict['Species']

    best_idx = 0
    best_distance = np.inf

    for i, compound_list in enumerate(peak_compounds):
        
        distance = objective_func(''.join(compound_list), exp_masses, exp_abundance, peak_mass, number_of_points)
        binding_dict['Closeness of Fit (Loss)'][i] = distance

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


def objective_func(formula, exp_masses, exp_abundance, peak_mass, number_of_points):
    '''
    Returns the DTW loss
    '''
    if formula == '':
        return 99999

    masses, x = find_isotope_pattern(formula, number_of_points)
    isotope_peak = masses[maxInBitonic(x, 0, number_of_points - 1)]
    difference = peak_mass - isotope_peak
    
    start = exp_masses.searchsorted(difference + masses[0])
    stop = exp_masses.searchsorted(difference + masses[-1], side='right')

    y = exp_abundance[start:stop]
    y = y / y.sum()
    # distance, _ = dtw(x, y, dist=None)

    m = exp_masses[start:stop]
    distance, _ = fastdtw(list(zip(masses, x)), list(zip(m, y)), dist=euclidean)

    return distance


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


def plotIsotopeDistribution(isotope_distribution, title="Isotope distribution"):
    '''
    Plot Isotope Distribution (https://pyopenms.readthedocs.io/en/latest/aasequences.html)
    '''
    plt.title(title)
    distribution = {"mass": [], "abundance": []}
    for iso in isotope_distribution.getContainer():
        distribution["mass"].append(iso.getMZ())
        distribution["abundance"].append(iso.getIntensity() * 100)

    bars = plt.bar(distribution["mass"], distribution["abundance"], width=0.05, snap=False) # snap ensures that all bars are rendered

    # plt.ylim([0, 110])
    # plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
    #                  math.ceil(distribution["mass"][-1]) + 2))
    plt.xlabel("Atomic mass (u)")
    plt.ylabel("Relative abundance (%)")
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

    # seq_formula = EmpiricalFormula(formula_str)
    # isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(generator_num) )
    # masses, rel_abundance = zip(*[(iso.getMZ(), iso.getIntensity()) for iso in isotopes.getContainer()])
    
    # if plot_distribution:
    #     plotIsotopeDistribution(isotopes)

    start = time.time()
    formula = 'C378H643Cl2N111O118Pt3S'
    formula = 'H'
    [print(peak_isotope(formula)) for _ in range(1)]
    print('Elapsed (seconds):', str((time.time()-start)))
    # print(peak_isotope('C378H629N105O118S1PtH2ONH3NH3') - 2*1.007825)
