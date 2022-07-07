
import gzip

# Replaced with: ln -sf resources/afdb/pdb results_af2/af2
#rule af2:
#    input:
#        pdb = 'resources/afdb/pdb/{struct_pref}/{struct_id}.pdb.gz',
#    output:
#        pdb = pfile(struct_id='{}', step='af2', suffix='.pdb.gz'),
#    shell: """
#        ln {input.pdb} {output.pdb}
#    """

configfile: 'config/config.yaml'

def af2_resid_pLDDT(fp_):
    resseq_pLDDT = collections.OrderedDict()
    parser = Bio.PDB.PDBParser(QUIET=True)
    fh_ = gzip.open(fp_, 'rt') # https://github.com/biopython/biopython/issues/3498#issuecomment-796563490
    structure = parser.get_structure(fp_, fh_)
    for chains in structure:
        for chain in chains:
            for residue in chain:
                resname = residue.get_resname()
                hetflag, resseq, icode = residue.get_id()
                for atom in residue:
                    resseq_pLDDT[resseq] = atom.bfactor
    return resseq_pLDDT

def af2_n_resid_mean_pLDDT(fp_):
    resseq_pLDDT = af2_resid_pLDDT(fp_)
    return len(resseq_pLDDT), np.mean(list(resseq_pLDDT.values()))

rule af2_bulk_stats:
    input:
        txt = 'results_af2/af2/{af2_bulk_id}.txt',
    output:
        tsv = 'results_af2/af2_bulk_stats/{af2_bulk_id}.tsv',
    resources:
        runtime = '02:00', # Runtime in hrs
        memory = '10000', # RAM in MB
    run:
        df_ = pd.read_csv(input.txt, names=['pdb_gz'])
        df_['pdb_gz'] = 'results_af2/af2' + df_['pdb_gz'].str.removeprefix('pdb')
        df_['accession'] = df_['pdb_gz'].str.removesuffix('.pdb.gz').map(os.path.basename)
        df_ = df_.set_index('accession', drop=True)
        df_[['af2_n_resid', 'af2_mean_pLDDT']] = [* map(af2_n_resid_mean_pLDDT, df_['pdb_gz']) ]
        df_.af2_n_resid = df_.af2_n_resid.astype(int)
        df_.to_csv(output.tsv, sep='\t', index=True, header=True, float_format='%.2f')

rule af2_bulk:
    """
    time snakemake af2_bulk --cores 1 --dry-run
    time snakemake af2_bulk --profile euler --cores 1 --dry-run
    """
    input:
        [f'results_af2/af2_bulk_stats/{af2_bulk_id}.tsv' for af2_bulk_id in config['af2_bulk_id'] ],

'''
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
'''
