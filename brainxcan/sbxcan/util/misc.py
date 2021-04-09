import pandas as pd
import os.path
import re

import numpy as np
import scipy.stats

def file_exists(fn):
    return os.path.isfile(fn)

def read_table(fn, indiv_col):
    _, fn_ext = os.path.splitext(fn)
    compress_args = {}
    if fn_ext == '.gz':
        fn_new = re.sub('.gz$', '', fn)
        compress_args = {'compression': 'gzip'}
        _, fn_ext = os.path.splitext(fn_new)
    if fn_ext == '.parquet':
        df = pd.read_parquet(fn)
    elif fn_ext == '.csv':
        df = pd.read_csv(fn, **compress_args)
    elif fn_ext == '.txt' or fn_ext == '.tsv':
        df = pd.read_csv(fn, sep='\s+', **compress_args)
    for i in range(df.shape[1]):
        if df.columns[i] == indiv_col:
            break
    col_list = df.columns.to_list()
    col_list.pop(i)
    col_list = [ indiv_col ] + col_list
    df = df.reindex(columns=col_list)
    df.rename(columns={indiv_col: 'indiv'}, inplace=True)
    df.indiv = df.indiv.astype(str)
    return df

def z2p(zscore):
    '''
    Input 1d np.array zscore and return the corresponding two-sided p-value.
    '''
    return scipy.stats.norm.sf(np.abs(zscore)) * 2
    