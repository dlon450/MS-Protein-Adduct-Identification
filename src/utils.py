from sklearn.preprocessing import MinMaxScaler 
import pandas as pd 
import numpy as np

from os import listdir
import config


def read(bound_file_path, unbound_file_path, compounds_file_path):
    '''
    Read in excel files as pd.DataFrame objects
    '''
    bound_df = pd.read_excel(bound_file_path)
    unbound_df = pd.read_excel(unbound_file_path)
    compounds_df = pd.read_excel(compounds_file_path)

    return bound_df, unbound_df, compounds_df


def normalise(spectrum):
    sc = MinMaxScaler()
    spectrum["normalised_intensity"] = sc.fit_transform(spectrum["I"].values.reshape(-1,1))

    return spectrum


def normalise_spectrums(unbound_df, bound_df):
    '''
    Normalise unbound/bound spectrums
    '''
    bound_df = normalise(bound_df)
    unbound_df = normalise(unbound_df)

    return unbound_df, bound_df