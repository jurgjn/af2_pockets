
rule lbsp:
    output:
        pdb = pfile(struct_id='{}', step='lbsp', suffix='.pdb'),
    shell: """
        gunzip -c resources/Clark2020/LBSp_dataset/data_files/{wildcards.struct_id}.pdb1.gz > {output.pdb} || gunzip -c resources/Clark2020/LBSp_dataset/CryptoSite_overlap_data_files/{wildcards.struct_id}.pdb1.gz > {output.pdb}
    """

'''
@functools.lru_cache()
def read_sifts_unique(v=True):
    """ Select families with a unique SIFTS mapping 
    """
    df_family_index = resources.family_index()
    df_pdb_chain_uniprot = resources.pdb_chain_uniprot()
    df_sifts_unique = df_family_index.merge(df_pdb_chain_uniprot, left_on='PDBid', right_on='PDB',).groupby(['Family']).agg(
        SP_PRIMARY_SET = ('SP_PRIMARY', lambda x: set(x).pop()),
        SP_PRIMARY_LEN = ('SP_PRIMARY', lambda x: len(set(x))),
    )
    df_sifts_unique = df_sifts_unique.query('SP_PRIMARY_LEN == 1').reset_index()
    if v: print(f'{len(df_sifts_unique)} protein families with unique SIFTS mappings selected')
    # Manual blacklist for a variety of different reasons, e.g. protein size or other difficulties modeling; one structure mapping to multiple LBSP "families"
    blacklist_ = ['C1KBQ3', 'Q97UY8', 'P24941', 'P16088', 'P68400', 'Q72498', 'P27708', 'P26660', 'P08659', 'Q16539', 'P23467', 'P23470', 'Q9Y6E0', 'Q8CHT0', 'P04191', 'O92972', 'P0C6X7', 'Q92793', 'O54438']\
               + ['P02754', 'P02766', 'P00183', 'P00760']
    return df_sifts_unique.query('~(SP_PRIMARY_SET in @blacklist_)')

checkpoint lbsp_af2_check_resid:
    """
    snakemake lbsp_af2_check_resid --cores $LSB_DJOB_NUMPROC --use-conda --dry-run
    """
    input:
        pdb = [ pfile(struct_id=struct_id, step='lbsp_af2', suffix='.pdb', base='results') for struct_id in read_sifts_unique(v=True)['SP_PRIMARY_SET'] ],
    output:
        tsv = 'results/lbsp_af2_check_resid.tsv',
    run:
        import Bio, Bio.PDB, resources
        df_ = read_sifts_unique(v=True).merge(resources.ubs_index(), left_on='Family', right_on='Fam_number')

        def res_(sp_primary_set, pos_):
            fp_ = pfile(struct_id=sp_primary_set, step='lbsp_af2', suffix='.pdb', base='results')
            structure = Bio.PDB.PDBParser().get_structure(sp_primary_set, fp_)
            ppb = Bio.PDB.PPBuilder()
            for pp in ppb.build_peptides(structure):
                return (pp[pos_ - 1].get_resname())

        df_['AF_res'] = [* map(res_, df_['SP_PRIMARY_SET'], df_['Res_number']) ]
        df_['res_ok'] = df_['Res'] == df_['AF_res']

        print(df_)
        print(df_['res_ok'].value_counts())

        df_['ubs_ok'] = df_.groupby('Family')['res_ok'].transform('all')
        print(df_[['SP_PRIMARY_SET', 'ubs_ok']].drop_duplicates()['ubs_ok'].value_counts())

        df_.query('ubs_ok == True').reset_index(drop=True).to_csv(output.tsv, sep='\t', index=False)

def fig4ABC_source_data_input(wildcards):
    # Final list of structures limited to the ones passing sanity checks (e.g. UBS residues)
    fp_ = checkpoints.lbsp_af2_check_resid.get(**wildcards).output[0]
    df_ = pd.read_csv(fp_, sep='\t')
    tsv_lbsp = [ pfile(struct_id=struct_id, step='lbsp.obabel_hxr.autosite.summary', suffix='.tsv', base='results') for struct_id in resources.family_index()['PDBid'] ]
    tsv_lbsp_af2 = [ pfile(struct_id=struct_id, step='lbsp_af2.obabel_hxr.autosite.summary', suffix='.tsv', base='results') for struct_id in set(df_['SP_PRIMARY_SET'].tolist()) ]
    return tsv_lbsp + tsv_lbsp_af2

rule fig4ABC_source_data:
    """
    ./smk_codon fig4ABC_source_data --dry-run
    snakemake fig4ABC_source_data --cores $LSB_DJOB_NUMPROC --use-conda --dry-run
    """
    input:
        fig4ABC_source_data_input
    output:
        tsv = 'fig4ABC.tsv'
    run:
        l_df_ = [* map(lambda fp: pd.read_csv(fp, sep='\t', dtype={'struct_id': 'str'}), input) ]
        df_ = pd.concat(l_df_, axis=0)
        df_.to_csv(output.tsv, index=False, sep='\t')
'''
