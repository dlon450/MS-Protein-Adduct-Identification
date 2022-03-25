import pandas as pd 
import numpy as np
import config
from binding_site_search import search
from os import listdir
from os.path import isfile, join


def search_all(dirpath, compounds_file, adducts_file):
    '''
    Search multiple files and save results as .csv
    '''
    calibrated = '_softwareCal' if config.calibrate == 'on' else '_noSWCal'
    bound_files = [f for f in listdir(dirpath) if isfile(join(dirpath, f))]
    for bound_fn in bound_files:
        binding_sites = search(join(dirpath, bound_fn), compounds_file, adducts_file)
        output_fn = bound_fn[:bound_fn.rfind('.')] + f'{calibrated}.csv'
        binding_sites.to_csv(join(dirpath, output_fn), index=False)


def validate_ground_truth(ground_truth_file, binding_sites):
    '''
    Validate results against ground truth
    '''
    gt = pd.read_csv(ground_truth_file)['Compounds'].to_numpy()
    pred = binding_sites['Species'].to_numpy()
    best = binding_sites['Closest Fit'].to_numpy()

    result = np.array(['M']*len(gt))             # missing, default
    result[np.isin(gt, pred[best])] = 'R'        # returned, best
    result[np.isin(gt, pred[~best])] = 'F'       # feasible, not best

    return result


def generate_results(ground_truth, bound, compounds, adducts, results_file='Data/results_validation_hires.csv',\
    obj_func='DTWEucl', weights=[0., 1., 10., 100., 1000., 10000.]):
    '''
    Generate results for objective function and range of given weights 
    '''
    df = pd.read_csv(results_file)
    for w in weights:
        binding_sites = search(bound, compounds, adducts, weight=w)
        results = validate_ground_truth(ground_truth, binding_sites)
        df[f'{obj_func}_{w}'] = results
    df.to_csv(results_file, index=False)


if __name__ == "__main__":
    
    compounds = "Data/Compound Constraints/Compounds_CisOxTrans_latest.xlsx"
    adducts = "Data/Compound Constraints/Standard_Adducts.xlsx"
    bound = "Data/Deconvoluted Spectra/uc_hires_precal.xlsx"
    ground_truth = "Data/ground_truth.csv"

    # search_all('Data/Deconvoluted Spectra', compounds, adducts)
    generate_results(ground_truth, bound, compounds, adducts)