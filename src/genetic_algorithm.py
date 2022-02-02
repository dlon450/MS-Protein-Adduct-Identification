from isotope_pattern import objective_func
import numpy as np
from pyeasyga import pyeasyga

from file_io import read
from isotope_pattern import average_mass
from icecream import ic

PROTON_MASS = 1.007825


if __name__ == '__main__':

    base_path = "Data/"
    fn = "Ubiquitin_plusC_1in100_000001"
    unbound_file_path = base_path + "Deconvoluted Spectra/Ubi_1in100_broad band_000001.xlsx"
    compounds_file_path = base_path + "Compound Constraints/Compounds_CisOxTrans_nlp.xlsx"
    bound_file_path = base_path + "Deconvoluted Spectra/" + fn + ".xlsx"

    bound_df, unbound_df, compounds = read(bound_file_path, unbound_file_path, compounds_file_path)
    exp_masses = bound_df['m/z'].to_numpy()
    exp_abundance = bound_df['I'].to_numpy()
    formulas = compounds["Formula"].to_numpy()
    charges = compounds["Charge"].to_numpy()
    max_amount = compounds["Max"].to_numpy()
    min_amount = compounds["Min"].to_numpy()
    masses = np.vectorize(average_mass)(formulas) - np.dot(PROTON_MASS, charges)

    masses_dict = dict(zip(formulas, masses))
    max_amount_dict = dict(zip(formulas, max_amount))
    min_amount_dict = dict(zip(formulas, min_amount))
    peak_mass = 8775.

    data = []
    for c in formulas:
        if max_amount_dict[c] == 1:
            data.extend([(c*i, masses_dict[c]*i) for i in range(min_amount_dict[c],2) if i != 0])
        elif max_amount_dict[c] <= 3:
            data.extend([(c*i, masses_dict[c]*i) for i in range(min_amount_dict[c],3) if i != 0])
        elif  max_amount_dict[c] <= 6:
            data.extend([(c*i, masses_dict[c]*i) for i in range(min_amount_dict[c],4) if i != 0])
    ic(data)

    def fitness(individual, data):
        mass = 8565.768014818
        formula = ''
        for (selected, item) in zip(individual, data):
            if selected:
                formula += item[0]
                mass += item[1]
        difference = abs(mass - peak_mass)
        if abs(mass - peak_mass) <= 4:
            return objective_func(formula, exp_masses, exp_abundance, peak_mass)
        else:
            return 99999. + difference

    # initialise the GA
    ga = pyeasyga.GeneticAlgorithm(data[:-1],
                                population_size=200,
                                generations=100,
                                elitism=True,
                                maximise_fitness=False)
    ga.fitness_function = fitness
    ga.run()
    print(ga.best_individual())
    # (107.79452715029723, [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0])