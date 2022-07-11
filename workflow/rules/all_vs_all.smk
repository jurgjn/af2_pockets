
rule afswiss_proteom_log10evalue:
    """ <2hr on a 
    time snakemake afswiss_proteom_log10evalue --cores 1 --dry-run
    time snakemake afswiss_proteom_log10evalue --profile euler --dry-run
    """
    input:
        m8_gz = 'resources/foldseek_all_vs_all/afswiss_proteom_ava.m8.gz',
    output:
        tsv_gz = 'results_ava/afswiss_proteom_ava.log10_evalue.tsv.gz',
        parquet = 'results_ava/afswiss_proteom_ava.log10_evalue.parquet',
    resources:
        runtime = '06:00', # Runtime in hrs
        memory = '500000', # RAM in MB
    params:
        #nrows = 1000000,
        nrows = None
    run:
        names_ = 'query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits'
        df_raw = pd.read_csv(input.m8_gz, sep='\t', names=names_, nrows=params.nrows)
        df_raw['query'] = df_raw['query'].str.replace(r'^AF-', '').str.replace(r'-F1-model_v2.pdb.gz$', '')
        df_raw['target'] = df_raw['target'].str.replace(r'^AF-', '').str.replace(r'-F1-model_v2.pdb.gz$', '')
        df_raw['log10_evalue'] = -np.log10(df_raw['evalue']).round(2)

        # log10 of evalue seems to overflow for particularly small values; set to maximum non-infinite value instead as a conservative lower bound
        max_ = df_raw['log10_evalue'].replace(np.inf, np.nan).dropna().max()
        df_raw['log10_evalue'] = df_raw['log10_evalue'].replace(np.inf, max_)

        df_raw = df_raw.sort_values(['log10_evalue'], ascending=False)
        
        df_raw.to_csv(output.tsv_gz, index=False, header=True, sep='\t')
        df_raw.to_parquet(output.parquet)

rule afswiss_proteom_cc:
    """
    time snakemake afswiss_proteom_cc --cores 1 --dry-run
    time snakemake afswiss_proteom_cc --profile euler --dry-run
    """
    input:
        tsv_gz = 'results_ava/afswiss_proteom_ava.log10_evalue.tsv.gz',
        #pkl = 'results_ava/afswiss_proteom_ava.log10_evalue.pkl',
    output:
        tsv_gz = 'results_ava/afswiss_proteom_ava.log10_evalue.cc.tsv.gz',
    resources:
        runtime = '12:00', # Runtime in hrs
        memory = '500000', # RAM in MB
    params:
        #nrows = 1000000,
        nrows = None
    run:
        print('Reading data')
        df_ = pd.read_csv(input.tsv_gz, sep='\t', nrows=params.nrows)
        #df_ = pd.read_pickle(input.pkl)

        print('Creating nodes')
        nodes_sorted = sorted(set(df_['query']) | set(df_['target']))
        G = nk.Graph(n=len(nodes_sorted), weighted=True, directed=False)
        G_accession = G.attachNodeAttribute('accession', str)
        for i, n in enumerate(nodes_sorted):
            G_accession[i] = n

        print('Adding edges')
        node_to_id = {n: i for i, n in enumerate(nodes_sorted)}
        for i, r in itertools.islice(df_.iterrows(), None):
            #print(i, r['query'], r['target'], r['log10_evalue'], se_[r['query']], se_[r['target']])
            G.addEdge(u=node_to_id[r['query']], v=node_to_id[r['target']], w=r['log10_evalue'])

        print('Finding connected components')
        cc = nk.components.ConnectedComponents(G)
        cc.run()

        print(nk.overview(G))
        def df_cc_():
            for c_i, c in itertools.islice(enumerate(cc.getComponents()), None):
                for node_i in c:
                    #print(node_i, G_accession[node_i], c_i)
                    yield(G_accession[node_i], c_i)

        df_cc = pd.DataFrame.from_records(df_cc_(), columns=['accession', 'connected_component_id'])
        df_cc.to_csv(output.tsv_gz, sep='\t')

rule afswiss_proteom_communities:
    """
    time snakemake afswiss_proteom_communities --cores 1 --dry-run
    time snakemake afswiss_proteom_communities --profile euler --dry-run
    """
    input:
        tsv_gz = 'results_ava/afswiss_proteom_ava.log10_evalue.tsv.gz',
        #pkl = 'results_ava/afswiss_proteom_ava.log10_evalue.pkl',
    output:
        tsv_gz = 'results_ava/afswiss_proteom_ava.log10_evalue.communities.tsv.gz',
    resources:
        cores = 10,
        runtime = '12:00', # Runtime in hrs
        memory = '50000', # RAM in MB
    params:
        #nrows = 10000,
        nrows = None
    run:
        nk.setNumberOfThreads(resources.cores)

        print('Reading data')
        df_ = pd.read_csv(input.tsv_gz, sep='\t', nrows=params.nrows)
        #df_ = pd.read_pickle(input.pkl)

        print('Creating nodes')
        nodes_sorted = sorted(set(df_['query']) | set(df_['target']))
        G = nk.Graph(n=len(nodes_sorted), weighted=True, directed=False)
        G_accession = G.attachNodeAttribute('accession', str)
        for i, n in enumerate(nodes_sorted):
            G_accession[i] = n

        print('Adding edges')
        node_to_id = {n: i for i, n in enumerate(nodes_sorted)}
        for i, r in itertools.islice(df_.iterrows(), None):
            #print(i, r['query'], r['target'], r['log10_evalue'], se_[r['query']], se_[r['target']])
            G.addEdge(u=node_to_id[r['query']], v=node_to_id[r['target']], w=r['log10_evalue'])

        # https://networkit.github.io/dev-docs/notebooks/Community.html
        print('Communities (PLM)')
        plmCommunities = nk.community.detectCommunities(G, algo=nk.community.PLM(G))
        print(plmCommunities)
        print("{0} elements assigned to {1} subsets".format(plmCommunities.numberOfElements(), plmCommunities.numberOfSubsets()))
        print("the biggest subset has size {0}".format(max(plmCommunities.subsetSizes())))

        print('Communities (PLP)')
        plpCommunities = nk.community.detectCommunities(G, algo=nk.community.PLP(G))
        print(plpCommunities)
        print("{0} elements assigned to {1} subsets".format(plpCommunities.numberOfElements(), plpCommunities.numberOfSubsets()))
        print("the biggest subset has size {0}".format(max(plpCommunities.subsetSizes())))

        # Output
        df_communities = pd.DataFrame.from_dict({
            'accession': nodes_sorted,
            'plm_community_id': plmCommunities.getVector(),
            'plp_commynity_id': plpCommunities.getVector(),
        })
        df_communities.to_csv(output.tsv_gz, sep='\t')

