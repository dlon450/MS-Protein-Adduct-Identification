import pandas as pd 
import numpy as np
import config
import matplotlib.pyplot as plt
from binding_site_search import search
from peak_search import peak_find
from isotope_pattern import peak_isotope
from os import listdir
from os.path import isfile, join
from utils import *


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


def search_paired_files(dirpath, adducts_file):
    '''
    Search through paired bound and compound files
    '''
    all_files = [f for f in listdir(dirpath) if isfile(join(dirpath, f))]
    all_files.sort()
    mid_idx = len(all_files) // 2
    bound_files = all_files[:mid_idx]
    compound_files = all_files[mid_idx:]

    for bound_fn, compound_fn in zip(bound_files, compound_files):
        binding_sites = search(join(dirpath, bound_fn), join(dirpath, compound_fn), adducts_file)
        output_fn = 'output_' + bound_fn[:bound_fn.rfind('.')] + '.csv'
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


def generate_results(ground_truth, bound, compounds, adducts, results_file,\
        obj_func='CL', weights=[0., 1., 10., 100., 1000., 10000.]):
    '''
    Generate results for objective function and range of given weights 
    '''
    df = pd.read_csv(results_file)
    for w in weights:
        binding_sites = search(bound, compounds, adducts, weight=w)
        results = validate_ground_truth(ground_truth, binding_sites)
        df[f'{obj_func}_{w}'] = results
    df.to_csv(results_file, index=False)


def plot_resolutions_MS(x_range=[8400, 9400], lr_file='Data/Deconvoluted Spectra/uc_lowres_precal.xlsx', \
        mr_file='Data/Deconvoluted Spectra/uc_medres_precal.xlsx', \
        hr_file='Data/Deconvoluted Spectra/uc_hires_precal.xlsx', plot_peaks=True, save=True):

    lr = pd.read_excel(lr_file)
    mr = pd.read_excel(mr_file)
    hr = pd.read_excel(hr_file)

    subplot = 311
    if plot_peaks:
        titles = ['LR', 'MR', 'HR']
        for i, df in enumerate([lr, mr, hr]):
            df = normalise(df)
            plt.subplot(subplot+i), plt.plot(df['m/z'], df['normalised_intensity']), plt.ylabel('Relative Abundance'), \
                plt.xlim(x_range), plt.title(titles[i])
            _, peaks, keep = peak_find(df, 0.01)
            plt.subplot(subplot+i), plt.plot(df['m/z'][peaks[keep]], df['normalised_intensity'][peaks[keep]], "kx")
        plt.xlabel('Atomic mass')
        if save:
            plt.savefig('Data/Deconvoluted Spectra/peak_identification_res.png')
    else:
        plt.subplot(311), plt.plot(lr['m/z'], lr['I']), plt.ylabel('Intensity'), \
            plt.xlim(x_range), plt.title('Low Resolution')
        plt.subplot(312), plt.plot(mr['m/z'], mr['I']), plt.ylabel('Intensity'), \
            plt.xlim(x_range), plt.title('Medium Resolution')
        plt.subplot(313), plt.plot(hr['m/z'], hr['I']), plt.ylabel('Intensity'), plt.xlabel('m'), \
            plt.xlim(x_range), plt.title('High Resolution')
        if save:
            plt.savefig('Data/Deconvoluted Spectra/MS_resolution.png')

    # plt.bar(hr['m/z'], hr['I'], width=0.05), plt.ylabel('Intensity'), plt.xlabel('m'), plt.xlim([8558.6, 8577.6])
    plt.show()


def accuracy_ppm_comparison(formulas=['C378H629N105O118S1'], n=8):
    '''
    Compare ppm at different accuracies for peak in isotope pattern
    '''
    total_ppms = np.zeros(n)
    for formula in formulas:
        baseline = peak_isotope(formula, accuracy=0.1 ** n)
        accuracies = 0.1 ** np.arange(1, n)
        for i, acc in enumerate(accuracies):
            current = peak_isotope(formula, acc)
            ppm = abs(current - baseline) / baseline * 1000000
            total_ppms[i] += ppm
    return total_ppms/n


if __name__ == "__main__":
    
    compounds = "Data/Compound Constraints/Compounds_CisOxTrans_latest.xlsx"
    adducts = "Data/Compound Constraints/Standard_Adducts.xlsx"
    bound = "Data/Deconvoluted Spectra/uc_hires_precal.xlsx"
    results = 'Data/results_validation_hires.csv'
    ground_truth = "Data/ground_truth.csv"

    # search_all('Data/Deconvoluted Spectra', compounds, adducts)
    # search_paired_files('Data/Input Data/Cisplatin-20220421T034342Z-001/Cisplatin', adducts)
    # generate_results(ground_truth, bound, compounds, adducts, results_file=results)
    # plt.rcParams["figure.figsize"] = (25,12)
    # plot_resolutions_MS(x_range=[8840, 8860], plot_peaks=False, save=False)
    print(accuracy_ppm_comparison(['C378H629N105O118S1', 'C560H874Fe1N148O156S4', 'C613H951O185N193S10', 'C769H1212N210O218S2']))
