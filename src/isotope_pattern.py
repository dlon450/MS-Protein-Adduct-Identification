from pyopenms import EmpiricalFormula, CoarseIsotopePatternGenerator
import numpy as np
from matplotlib import pyplot as plt
from icecream import ic
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw


def find_isotope_pattern(formula_str, generator_num=20, plot_distribution=False):
    '''
    Return theoretical distribution of intensities
    '''
    seq_formula = EmpiricalFormula(formula_str)
    isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(generator_num) )
    masses, rel_abundance = zip(*[(iso.getMZ(), iso.getIntensity()*100) for iso in isotopes.getContainer()])
    
    if plot_distribution:
        plotIsotopeDistribution(isotopes)

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

        distance, _ = fastdtw(theo_rel_abundance, exp_rel_abundance, dist=euclidean)
        binding_df.loc[idx[i], 'Loss'] = distance
        
        if distance < best_distance:
            best_distance = distance
            best_idx = idx[i]

    # plotWarpDTW(theo_rel_abundance, exp_rel_abundance, warp_path)
    binding_df.loc[best_idx, 'Best'] = True
        

def calculate_score_no_interpolation(f, peak_mass, search_mask, binding_df, bound_df):
    '''
    Return the compounds that best match the experimental distribution
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
        x = np.array(list(zip(masses, theo_rel_abundance)))
        y = bound_df[abs(bound_df['m/z'] - peak_mass) <= interval][['m/z', 'I']].to_numpy()
        y[:,1] = y[:,1] / y[:,1].sum() * 100
        
        distance, _ = fastdtw(x[:,1], y[:,1], dist=euclidean)
        binding_df.loc[idx[i], 'Loss'] = distance

        if distance < best_distance:
            best_distance = distance
            best_idx = idx[i]
    
    # plotWarpDTW(x[:,1], y[:,1], warp_path)
    binding_df.loc[best_idx, 'Best'] = True


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

    formula = 'C378H630N105O118S1H2OClK'
    masses, rel_abundance = find_isotope_pattern(formula, plot_distribution=True)
    for mass, rel_adun in zip(masses, rel_abundance):
        print("Isotope", mass, "has abundance", rel_adun, "%")