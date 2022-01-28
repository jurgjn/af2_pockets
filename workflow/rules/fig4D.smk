
@functools.lru_cache()
def read_af2_not_swiss():
    df_af2_ = resources.read_afdb().rename({'uniprot_': 'UniProtKB_ac'}, axis=1)
    in_swiss_ = set(resources.read_swiss()['UniProtKB_ac'])
    print(uf(len(in_swiss_)), 'UniProtKB_ac in SWISS-MODEL index')
    df_af2_['in_swiss'] = [* map(lambda uniprot_: uniprot_ in in_swiss_, df_af2_['UniProtKB_ac'])]
    df_af2_['in_swiss'].value_counts()
    df_af2_ = df_af2_.query('in_swiss == False')
    assert any(df_af2_['UniProtKB_ac'].duplicated()) == False
    print(uf(len(df_af2_)), 'single-pdb AF2 models not in SWISS-MODEL (checked unique UniProtKB_ac)')
    return df_af2_
    #return df_af2_.head(10)
    #blacklist_ = ['Q5VWN6', 'P46821'] # Two >2000 aa structures with slow pocket detection
    #return df_af2_.query('~(UniProtKB_ac in @blacklist_)')#.head(10)

rule af2:
    # Small number of structures have `v2` structures, e.g. Q9BXP8, Q13219
    output:
        pdb = pfile(struct_id='{}', step='af2', suffix='.pdb'),
    shell: """
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v1.pdb ||\
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v2.pdb
    """

rule af2_not_swiss:
    """
    snakemake af2_not_swiss --cores $LSB_DJOB_NUMPROC --use-conda --dry-run
    """
    input:
        pdb = [ pfile(struct_id=struct_id, step='af2', suffix='.pdb', base='results') for struct_id in read_af2_not_swiss()['UniProtKB_ac'] ],
    output:
        tsv = 'results/af2_not_swiss.tsv',
    run:
        def get_af2_stats(fp_):
            resseq_pLDDT = get_resid_pLDDT(fp_)
            return len(resseq_pLDDT), np.mean(list(resseq_pLDDT.values()))

        col_ = ['UniProtKB_ac', 'n_resid', 'mean_pLDDT']
        df_ = read_af2_not_swiss()
        df_[['n_resid', 'mean_pLDDT']] = [* map(get_af2_stats, input.pdb) ]
        df_[col_].to_csv(output.tsv, sep='\t', index=False, header=True, float_format='%.2f')

rule uniprot_txt:
    output:
        txt = pfile(struct_id='{}', step='uniprot_txt', suffix='.txt'),
    shell: """
        wget -O {output.txt} https://www.uniprot.org/uniprot/{wildcards.struct_id}.txt
    """

def has_catalytic_activity(fp):
    try:
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec170
        with open(fp) as fh:
            rec = Bio.SwissProt.read(fh)
            #print(dir(rec))
            for _ in rec.comments:
                if _.startswith('CATALYTIC ACTIVITY:'):
                    #print(_)
                    return _
    except ValueError as e:
        print(fp, e)
    return False

def has_EC(fp):
    try:
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec170
        with open(fp) as fh:
            rec = Bio.SwissProt.read(fh)
            #print(dir(rec))
            if '; EC=' in rec.description:
                return rec.description
    except ValueError as e:
        print(fp, e)
    return False

rule fig4D_source_data:
    """
    snakemake fig4D_source_data --cores $LSB_DJOB_NUMPROC --use-conda --dry-run
    """
    input:
        af2_not_swiss = 'results/af2_not_swiss.tsv',
        txt = [ pfile(struct_id=struct_id, step='uniprot_txt', suffix='.txt', base='results') for struct_id in read_af2_not_swiss()['UniProtKB_ac'] ],
        tsv = [ pfile(struct_id=struct_id, step='af2.obabel_hxr.autosite.summary', suffix='.tsv', base='results') for struct_id in read_af2_not_swiss()['UniProtKB_ac'] ],
    output:
        tsv = 'fig4D.tsv'
    run:
        df_ = pd.read_csv(input.af2_not_swiss, sep='\t')#.head(10)
        df_['has_catalytic_activity'] = [ *map(has_catalytic_activity, input.txt) ]
        df_['has_EC'] = [ *map(has_EC, input.txt) ]
        df_['has_activity'] = (df_['has_catalytic_activity'] != False) | (df_['has_EC'] != False)

        def score_(fp_, mean_pLDDT_thresh_=0):
            df_ = pd.read_csv(fp_, sep='\t').query('mean_pLDDT >= @mean_pLDDT_thresh_')
            if len(df_) > 0:
                return df_.head(1)['score'].tolist()[0]
            else:
                return 0

        df_['autosite_score'] = [ *map(score_, input.tsv) ]
        df_['autosite_score_50'] = [ *map(lambda fp_: score_(fp_, mean_pLDDT_thresh_=50), input.tsv) ]
        df_['autosite_score_70'] = [ *map(lambda fp_: score_(fp_, mean_pLDDT_thresh_=70), input.tsv) ]
        df_['autosite_score_90'] = [ *map(lambda fp_: score_(fp_, mean_pLDDT_thresh_=90), input.tsv) ]
        df_['autosite_score_95'] = [ *map(lambda fp_: score_(fp_, mean_pLDDT_thresh_=95), input.tsv) ]

        df_.to_csv(output.tsv, sep='\t', index=False, float_format='%.2f')
