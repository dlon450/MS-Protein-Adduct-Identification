from sklearn.preprocessing import MinMaxScaler 
import pandas as pd 
import numpy as np

def normalise(spectrum):
    sc = MinMaxScaler()
    spectrum["normalised_intensity"] = sc.fit_transform(spectrum["I"].values.reshape(-1,1))

    return spectrum

def normalise_spectrums(unbound_df, bound_df):
    bound_df = normalise(bound_df)
    unbound_df = normalise(unbound_df)

    return unbound_df, bound_df