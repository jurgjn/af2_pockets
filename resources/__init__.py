
import itertools
import os
import pandas as pd

from workflow.modules import *

@functools.lru_cache()
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

@functools.lru_cache()
def read_swiss(v=False):
    fp_ = 'resources/SWISS_2021_11_30/SWISS-MODEL_Repository/INDEX'
    df_ = pd.read_csv(fp_, comment='#', sep='\t')
    if v: print(f'{len(df)}\traw records')
    return df_

@functools.lru_cache()
def family_index():
    fp_ = 'resources/Clark2020/LBSp_dataset/family_index.csv'
    df_ = pd.read_csv(fp_, skiprows=9).iloc[: , 1:]
    df_.columns = ['PDBid', 'Family', 'Type', 'In_CryptoSite', 'Resolution']
    return df_

@functools.lru_cache()
def ubs_index():
    fp_ = 'resources/Clark2020/LBSp_dataset/UBS_index.csv'
    df_ = pd.read_csv(fp_, skiprows=4).iloc[: , 1:]
    df_.columns = ['Fam_number', 'Res', 'Res_number', 'Chain']
    return df_

@functools.lru_cache()
def perf_struct():
    fp_ = 'resources/Clark2020/41598_2020_72906_MOESM2_ESM.txt'
    return pd.read_csv(fp_, sep='\t', index_col=0)

@functools.lru_cache()
def perf_family():
    fp_ = 'resources/Clark2020/41598_2020_72906_MOESM3_ESM.txt'
    return pd.read_csv(fp_, sep='\t', index_col=0)

@functools.lru_cache()
def pdb_chain_uniprot():
    fp_ = 'resources/sifts/pdb_chain_uniprot.tsv'
    return pd.read_csv(fp_, sep='\t', comment='#')
