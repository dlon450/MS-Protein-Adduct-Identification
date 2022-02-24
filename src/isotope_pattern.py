from pyopenms import *
import numpy as np
from matplotlib import pyplot as plt
from icecream import ic
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw, dtw


def peak_isotope(formula: object) -> float:
    '''
    Use isotope pattern to find average mass
    '''
    seq_formula = EmpiricalFormula(formula)
    isotopes = seq_formula.getIsotopeDistribution( FineIsotopePatternGenerator(1e-3) )
    iso_dict = {}
    for iso in isotopes.getContainer():
        mz = iso.getMZ()
        group = int(mz)
        iso_dict[group] = iso_dict.get(group, []) + [(mz, iso.getIntensity())]

    best_abundance = 0.
    best_mass = 0.
    for _, value in iso_dict.items():
        abundance = sum([i for m,i in value])
        if abundance > best_abundance:
            best_mass = sum([ m*i for m,i in value]) / abundance
            best_abundance = abundance

    return best_mass


def find_isotope_pattern(formula_str: str, generator_num=16, paired=False, plot_distribution=False):
    '''
    Return theoretical distribution of intensities
    '''
    seq_formula = EmpiricalFormula(formula_str)
    isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(generator_num) )
    masses, rel_abundance = zip(*[(iso.getMZ(), iso.getIntensity()) for iso in isotopes.getContainer()])
    
    if plot_distribution:
        plotIsotopeDistribution(isotopes)

    # isotope_peak = seq_formula.getMonoWeight()
    # average_mass = seq_formula.getAverageWeight()

    return masses, rel_abundance


def calculate_score(f, peak_mass, search_mask, binding_df):
    '''
    Find the distances between experimental and theoretical distributions using FastDTW https://github.com/slaypni/fastdtw
    '''    
    number_of_points = 20
    peak_compounds_df = binding_df[search_mask]
    idx = peak_compounds_df.index
    peak_compounds_df_list = peak_compounds_df.iloc[:,:2].to_dict('records')

    best_idx = 0
    best_distance = np.inf

    for i, compound_dict in enumerate(peak_compounds_df_list):
        masses, theo_rel_abundance = find_isotope_pattern(''.join(compound_dict['Compound']), number_of_points)
        interval = (masses[-1] - masses[0])/2
        exp_abundance = f(np.linspace(peak_mass-interval, peak_mass+interval, num=number_of_points))
        exp_rel_abundance = exp_abundance / exp_abundance.sum() * 100

        distance, _ = dtw(theo_rel_abundance, exp_rel_abundance, dist=None)
        binding_df.loc[idx[i], 'Loss'] = distance
        
        if distance < best_distance:
            best_distance = distance
            best_idx = idx[i]

    # plotWarpDTW(theo_rel_abundance, exp_rel_abundance, warp_path)
    binding_df.loc[best_idx, 'Best'] = True
        

def calculate_score_no_interpolation(peak_mass, binding_dict, bound_df, full=False):
    '''
    Return the compounds that best match the experimental distribution
    '''
    # find sites by matching
    binding_site_record = dict.fromkeys(binding_dict.keys())

    exp_masses = bound_df['m/z'].to_numpy()
    exp_abundance = bound_df['I'].to_numpy()

    number_of_points = 16
    peak_compounds = binding_dict['Compound']

    best_idx = 0
    best_distance = np.inf

    for i, compound_list in enumerate(peak_compounds):
        
        distance = objective_func(''.join(compound_list), exp_masses, exp_abundance, peak_mass, number_of_points)
        binding_dict['Loss'][i] = distance

        if distance < best_distance:
            best_distance = distance
            best_idx = i
    
    # plotWarpDTW(x, y, warp_path)
    binding_dict['Best'][best_idx]= True
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
    distance, _ = dtw(list(zip(masses, x)), list(zip(m, y)), dist=euclidean)

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


if __name__ == "__main__":

    formula = 'C378H629N105O118S1PtNH3NH3H2O'
    print(peak_isotope(formula))