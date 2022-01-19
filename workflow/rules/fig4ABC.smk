


rule lbsp:
    output:
        pdb = pfile(struct_id='{}', step='lbsp', suffix='.pdb'),
    shell: """
        gunzip -c resources/Clark2020/LBSp_dataset/data_files/{wildcards.struct_id}.pdb1.gz > {output.pdb} || gunzip -c resources/Clark2020/LBSp_dataset/CryptoSite_overlap_data_files/{wildcards.struct_id}.pdb1.gz > {output.pdb}
    """

rule fig1A_source_data:
    """
    snakemake fig1A_source_data --cores $LSB_DJOB_NUMPROC --use-conda --dry-run
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
