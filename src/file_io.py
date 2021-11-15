import pandas as pd
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