
rule afswiss_proteom_reduced_log10evalue:
    """
    time snakemake afswiss_proteom_ava_accessions --cores 1 --dry-run
    time snakemake --profile euler afswiss_proteom_ava_accessions --cores 1 --dry-run
    """
    input:
        m8_gz = 'resources/foldseek_all_vs_all/afswiss_proteom_ava.m8.gz',
        #output:
        #tsv_gz = 'results_ava/afswiss_proteom_ava.reduced_log10_evalue.m8.gz',
    resources:
        runtime = '01:00', # Runtime in hrs
        memory = '256000', # RAM in MB
    params:
        #nrows = 10000,
        nrows = None,
    run:
        names_ = 'query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits'
        df_raw = pd.read_csv(input.m8_gz, sep='\t', names=names_, nrows=params.nrows)
        df_raw['query'] = df_raw['query'].str.replace(r'^AF-', '').str.replace(r'-F1-model_v2.pdb.gz$', '')
        df_raw['target'] = df_raw['target'].str.replace(r'^AF-', '').str.replace(r'-F1-model_v2.pdb.gz$', '')
        df_raw['log10_evalue'] = (-np.log10(df_raw['evalue'])).replace([np.inf, -np.inf], np.nan).round(2)
        print(df_raw)
        #se_ids = pd.concat([df_raw['query'], df_raw['target']], axis=0).drop_duplicates(keep='first').sort_values()
        #se_ids.to_csv(output.tsv_gz, index=False, header=False, sep='\t')


