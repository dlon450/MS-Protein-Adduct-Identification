# from lib2to3.pgen2.parse import ParseError
import pandas as pd 
import numpy as np
import config
import matplotlib.pyplot as plt
from binding_site_search import search
from peak_search import peak_find
from isotope_pattern import peak_isotope
from os import listdir
from os.path import isfile, join, dirname
from utils import normalise
import sys
# import csv


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


def search_paired_files(dirpath, adducts_file, analysis=True):
    '''
    Search through paired bound and compound files
    '''
    all_files = [f for f in listdir(dirpath) if isfile(join(dirpath, f))]
    all_files.sort()
    mid_idx = len(all_files) // 2
    bound_files = all_files[:mid_idx]
    compound_files = all_files[mid_idx:]

    if not analysis:
        for bound_fn, compound_fn in zip(bound_files, compound_files):
            sys.stdout = open(join(dirpath, 'peaks', bound_fn[:bound_fn.rfind('.')] + '.txt'),'wt')
            binding_sites = search(join(dirpath, bound_fn), join(dirpath, compound_fn), adducts_file)
            output_fn = 'output_' + bound_fn[:bound_fn.rfind('.')] + '.csv'
            binding_sites.to_csv(join(dirpath, output_fn), index=False)
    else:
        ground_truth_files = listdir(join(dirpath, 'ground_truth'))
        ground_truth_files.sort()
        tolerances = [3.1]
        peak_heights = [0.01]
        # first = True
        for bound_fn, compound_fn, ground_truth in zip(bound_files, compound_files, ground_truth_files):
            # if first:
            #     first = False
            #     continue
            for ph in peak_heights:
                for tol in tolerances:
                    generate_results(join(dirpath, 'ground_truth', ground_truth), join(dirpath, bound_fn), join(dirpath, compound_fn), \
                        adducts_file, join(dirpath, 'ground_truth', ground_truth), column_name='cDTW_CB', tolerance=tol, peak_height=ph)
          

def validate_ground_truth(ground_truth_file, binding_sites, excel_file=False):
    '''
    Validate results against ground truth
    '''
    if not excel_file:
        gt = pd.read_csv(ground_truth_file)['Compounds'].to_numpy()
    else:
        gt = pd.read_excel(ground_truth_file)['Compounds'].to_numpy()
    pred = binding_sites['Species'].to_numpy()
    best = binding_sites['Closest Fit'].to_numpy()

    result = np.array(['M']*len(gt))             # missing, default
    result[np.isin(gt, pred[best])] = 'R'        # returned, best
    result[np.isin(gt, pred[~best])] = 'F'       # feasible, not best

    return result


def generate_results(ground_truth, bound, compounds, adducts, results_file,\
        column_name='CL', weights=[0., 1., 10., 100., 1000., 10000.], tolerance=config.tolerance, peak_height=config.peak_height):
    '''
    Generate results for objective function and range of given weights 
    '''
    extension = results_file[results_file.rfind('.'):]
    if extension == '.csv':
        df = pd.read_csv(results_file)
        excel = False
    elif extension == '.xlsx':
        df = pd.read_excel(results_file)
        excel = True
    else:
        raise ValueError('Incorrect file type.')

    for w in weights:
        binding_sites = search(bound, compounds, adducts, weight=w, tolerance=tolerance, peak_height=peak_height)
        results = validate_ground_truth(ground_truth, binding_sites, excel)
        df[f'{column_name}_{w}'] = results
    
    if excel:
        df.to_excel(results_file, index=False)
    else:
        df.to_csv(results_file, index=False)


def plot_resolutions_MS(x_range=[8400, 9400], lr_file='Data/Deconvoluted Spectra/uc_lowres_precal.xlsx', \
        mr_file='Data/Deconvoluted Spectra/uc_medres_precal.xlsx', \
        hr_file='Data/Deconvoluted Spectra/uc_hires_precal.xlsx', plot_peaks=True, save=True):

    lr = pd.read_excel(lr_file)
    mr = pd.read_excel(mr_file)
    hr = pd.read_excel(hr_file)

    subplot = 311
    if plot_peaks:
        titles = ['Low Resolution', 'Medium Resolution', 'High Resolution']
        for i, df in enumerate([lr, mr, hr]):
            df = normalise(df)
            plt.subplot(subplot+i), plt.plot(df['m/z'], df['normalised_intensity']), plt.ylabel('Relative Abundance'), \
                plt.xlim(x_range), plt.title(titles[i])
            _, peaks, keep, _ = peak_find(df, 0.01)
            plt.subplot(subplot+i), plt.plot(df['m/z'][peaks[keep]], df['normalised_intensity'][peaks[keep]], "kx")
        plt.xlabel('Atomic mass')
        if save:
            plt.savefig('Data/Deconvoluted Spectra/peak_identification_res.pdf')
    else:
        plt.subplot(311), plt.plot(lr['m/z'], lr['I']), plt.ylabel('Intensity'), \
            plt.xlim(x_range), plt.title('Low Resolution')
        plt.subplot(312), plt.plot(mr['m/z'], mr['I']), plt.ylabel('Intensity'), \
            plt.xlim(x_range), plt.title('Medium Resolution')
        plt.subplot(313), plt.plot(hr['m/z'], hr['I']), plt.ylabel('Intensity'), plt.xlabel('m'), \
            plt.xlim(x_range), plt.title('High Resolution')
        if save:
            plt.savefig('Data/Deconvoluted Spectra/MS_resolution.pdf')

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


def find_peaks_with_solutions(dirpath, adducts_file):
    '''
    Search through paired bound and compound files
    '''
    all_files = [f for f in listdir(dirpath) if isfile(join(dirpath, f))]
    all_files.sort()
    mid_idx = len(all_files) // 2
    bound_files = all_files[:mid_idx]
    compound_files = all_files[mid_idx:]
    ground_truth_files = listdir(join(dirname(dirpath), 'ground_truth'))
    ground_truth_files.sort()

    tolerances = [1.1, 2.1, 3.1, 4.1, 5.1, 6.1]
    peak_heights = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1]
    fp = {'Dataset': []}
    fp.update({f'{ph*100.}%_{tol}': [] for ph in peak_heights for tol in tolerances})

    for bound_fn, compound_fn, ground_truth in zip(bound_files, compound_files, ground_truth_files):

        extension = ground_truth[ground_truth.rfind('.'):]
        if extension == '.csv':
            gt = pd.read_csv(join(dirname(dirpath), 'ground_truth', ground_truth))
        elif extension == '.xlsx':
            gt = pd.read_excel(join(dirname(dirpath), 'ground_truth', ground_truth))
        else:
            raise ValueError('Incorrect file type.')
        
        df = {}
        df['Compounds'] = gt['Compounds'].to_numpy()
        if 'Experimental Mass' in gt.columns:
            df['Experimental Mass'] = gt['Experimental Mass'].to_numpy()
            f, _, _ = search(join(dirpath, bound_fn), join(dirpath, compound_fn), adducts_file, \
                    tolerance=3.1, peak_height=0.01, weight=1., return_peaks=True)
            df['Intensities (%)'] = f(df['Experimental Mass']).round(5) * 100.

        fp['Dataset'].append(f'{ground_truth[13:]}')
        for ph in peak_heights:
            for tol in tolerances:
                f, all_peaks, peaks = search(join(dirpath, bound_fn), join(dirpath, compound_fn), adducts_file, \
                    tolerance=tol, peak_height=ph, weight=1., return_peaks=True)
                df[f'{ph*100.}%_{tol}'] = peaks

                # fp[f'{ph*100.}%_{tol}'].append(len(peaks)) # total peaks

                if 'Experimental Mass' in df:
                    count = 0
                    for gp in df['Experimental Mass']:
                        if not np.isclose(peaks, gp, atol=5.).any(): # not this will be missing
                            count += 1
                    # fp[f'{ph*100.}%_{tol}'].append(len(peaks) - count)
                    fp[f'{ph*100.}%_{tol}'].append(count)

                else:
                    fp[f'{ph*100.}%_{tol}'].append(0)
                
        # print(df)
        # pd.DataFrame.from_dict(df, orient='index').transpose().to_csv(join(dirpath, 'peaks', '_' + bound_fn[:bound_fn.rfind('.')] + '.csv'), index=False)
    return fp

def insert_intensities(dirpath):
    '''
    Insert intensities column into results
    '''
    peak_solns_files = listdir(join(dirpath, 'peaks'))
    peak_solns_files.sort()

    results_files = listdir(join(dirpath, 'ground_truth'))
    results_files.sort()

    missing = []
    for r, p in zip(results_files, peak_solns_files):
        df = pd.read_excel(join(dirpath, 'ground_truth', r))
        try:
            if 'Intensities (%)' not in df.columns:
                df['Intensities (%)'] = pd.read_csv(join(dirpath, 'peaks', p))['Intensities (%)'].dropna()
                cols = df.columns.tolist()
                cols = [cols[0]] + [cols[-1]] + cols[1:-1]
                df[cols].to_excel(join(dirpath, 'ground_truth', r), index=False)
            else:
                df['Intensities (%)'] = pd.read_csv(join(dirpath, 'peaks', p))['Intensities (%)'].dropna()
                df.to_excel(join(dirpath, 'ground_truth', r), index=False)
        except KeyError:
            missing.append(p)
            continue
    
    if missing != []:
        print(missing)


def find_peak_distances_distribution(dirpath):
    # minimum is 15.98 --> use 15
    ground_truth_files = listdir(join(dirname(dirpath), 'ground_truth'))
    distances = []

    for ground_truth in ground_truth_files:
        extension = ground_truth[ground_truth.rfind('.'):]
        if extension == '.csv':
            gt = pd.read_csv(join(dirname(dirpath), 'ground_truth', ground_truth))
        elif extension == '.xlsx':
            gt = pd.read_excel(join(dirname(dirpath), 'ground_truth', ground_truth))
        else:
            raise ValueError('Incorrect file type.')
        
        if 'Experimental Mass' in gt.columns:
            p = gt['Experimental Mass'].to_numpy()
            distances.append([t - s for s, t in zip(
                p, p[1:])])

    return distances


if __name__ == "__main__":
    
    compounds = "Data/Compound Constraints/Compounds_CisOxTrans_latest.xlsx"
    adducts = "Data/Compound Constraints/Standard_Adducts.xlsx"
    bound = "Data/Deconvoluted Spectra/uc_hires_precal.xlsx"
    results = 'Data/results_validation_hires.csv'
    ground_truth = "Data/ground_truth.csv"

    folders = [
        r'Data\Input Data\Cisplatin-20220421T034342Z-001\Cisplatin',
        r'Data\Input Data\Gold Complexes-20220421T034327Z-001\Gold Complexes',
        r'Data\Input Data\Oxaliplatin-20220421T034341Z-001\Oxaliplatin',
        r'Data\Input Data\RAPTA-C-20220421T034334Z-001\RAPTA-C',
        r'Data\Input Data\Ru Complexes-20220421T034332Z-001\Ru Complexes',
        r'Data\Input Data\ox_rapc_hires\inner'
    ]
    # folders = [r'Data\Input Data\ox_rapc_hires\inner']

    # search_all('Data/Deconvoluted Spectra', compounds, adducts)
    # search_paired_files(r'Data\Input Data\RAPTA-C-20220421T034334Z-001\RAPTA-C', adducts, analysis=True)
    # false_positives = {}
    # for folder in folders:
    #     # search_paired_files(folder, adducts, analysis=True)

    #     fp = find_peaks_with_solutions(folder, adducts)
    #     for col, val in fp.items():
    #         if col not in false_positives:
    #             false_positives[col] = val
    #         else:
    #             false_positives[col].extend(val)

    #     # insert_intensities(folder)
    # pd.DataFrame.from_dict(false_positives, orient='index').transpose().to_csv(r'Data\Input Data\missing_peaks.csv', index=False)
    plt.rcParams["figure.figsize"] = (15,12)
    plot_resolutions_MS()

    # distances = [find_peak_distances_distribution(folder) for folder in folders]
    # print(distances)

    # generate_results(ground_truth, bound, compounds, adducts, results_file=results)
    # plt.rcParams["figure.figsize"] = (25,12)
    # plot_resolutions_MS(x_range=[8840, 8860], plot_peaks=False, save=False)
    # print(accuracy_ppm_comparison(['C378H629N105O118S1', 'C560H874Fe1N148O156S4', 'C613H951O185N193S10', 'C769H1212N210O218S2']))
