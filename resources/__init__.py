
import itertools
import os
import pandas as pd

from workflow.modules import *

def read_afdb(single_pdb_only=True, v=True):
    fp_ = 'resources/afdb/UP000005640_9606_HUMAN.txt'
    df_ = pd.read_csv(fp_, names=('pdb_',))
    if v: print(uf(len(df_)), 'raw-structures')
    df_.insert(loc=0, column='uniprot_', value=df_['pdb_'].str.split('-').map(lambda l: l[1]))
    if single_pdb_only:
        #df_ = df_.query('pdb_.str.endswith("-F1-model_v1.pdb.gz")', engine='python').copy()
        df_ = df_[ ~df_['uniprot_'].duplicated(keep=False) ].reset_index(drop=True)
        if v: print(uf(len(df_)), 'single-pdb structures')
    return df_

def read_swiss(v=False):
    fp_ = 'resources/SWISS_2021_11_30/SWISS-MODEL_Repository/INDEX'
    df_ = pd.read_csv(fp_, comment='#', sep='\t')
    if v: print(f'{len(df)}\traw records')
    return df_
