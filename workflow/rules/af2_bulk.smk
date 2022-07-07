
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

localrules: af2_wget # compute nodes don't have Internet access

# Replaced with: ln -sf resources/afdb/pdb results_af2/af2
#rule af2:
#    input:
#        pdb = 'resources/afdb/pdb/{struct_pref}/{struct_id}.pdb.gz',
#    output:
#        pdb = pfile(struct_id='{}', step='af2', suffix='.pdb.gz'),
#    shell: """
#        ln {input.pdb} {output.pdb}
#    """

rule af2_wget: # https://www.alphafold.ebi.ac.uk/download
    output:
        tar = 'resources/afdb/{af2_bulk_id}.tar',
    shell: """
        wget -O {output.tar} https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/{wildcards.af2_bulk_id}.tar
    """

rule af2_pdb:
    """Extract AF2 bulk archive into prefix-derived subdirectories to keep individual directory sizes manageable

    Testing:
        tar -tf UP000000625_83333_ECOLI_v2.tar --exclude='*.cif.gz' --xform='s|\(^AF-\)\(..\)\(..\)\(..\)\(.*\)\(-F1-model_v2.pdb.gz$\)|pdb/\2/\3/\4/\2\3\4\5.pdb.gz|' --verbose --show-transformed-names

    Refs:
    https://unix.stackexchange.com/questions/198151/tar-extract-into-directory-with-same-base-name
    https://www.gnu.org/software/sed/manual/html_node/Back_002dreferences-and-Subexpressions.html
    """
    input:
        tar = 'resources/afdb/{af2_bulk_id}.tar',
    output:
        txt = 'resources/afdb/pdb/{af2_bulk_id}.txt',
        txt_ln = 'results_af2/af2/{af2_bulk_id}.txt',
    resources:
        runtime = '04:00', # Runtime in hrs
        memory = '10000', # RAM in MB
    shell: """
        tar -xf {input.tar} --exclude='*.cif.gz' --xform='s|\\(^AF-\\)\\(..\\)\\(..\\)\\(..\\)\\(.*\\)\\(-model_v2.pdb.gz$\\)|resources/afdb/pdb/\\2/\\3/\\4/\\2\\3\\4\\5.pdb.gz|'
        tar -tf {input.tar} --exclude='*.cif.gz' --xform='s|\\(^AF-\\)\\(..\\)\\(..\\)\\(..\\)\\(.*\\)\\(-model_v2.pdb.gz$\\)|resources/afdb/pdb/\\2/\\3/\\4/\\2\\3\\4\\5.pdb.gz|' --show-transformed-names > {output.txt}
    """

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
        df_['pdb_gz'] = 'results_af2/af2' + df_['pdb_gz'].str.removeprefix('resources/afdb/pdb')
        df_['accession'] = df_['pdb_gz'].str.removesuffix('.pdb.gz').map(lambda s: os.path.basename(s.rsplit('-', 1)[0]))
        df_['af2_model_id'] = df_['pdb_gz'].str.removesuffix('.pdb.gz').map(os.path.basename)
        df_ = df_.set_index('af2_model_id', drop=True)
        df_[['af2_n_resid', 'af2_mean_pLDDT']] = [* map(af2_n_resid_mean_pLDDT, df_['pdb_gz']) ]
        df_.af2_n_resid = df_.af2_n_resid.astype(int)
        df_.to_csv(output.tsv, sep='\t', index=True, header=True, float_format='%.2f')

rule af2_bulk_prank:
    """
    software/p2rank_2.4/prank predict -c alphafold -o results_af2/af2_bulk_p2rank UP000005640_9606_HUMAN_v2.txt
    """
    input:
        txt = 'resources/afdb/pdb/{af2_bulk_id}.txt',
    output:
        txt = temp('{af2_bulk_id}.txt'),
        dir = directory('results_af2/af2_bulk_prank/{af2_bulk_id}'),
    resources:
        runtime = '12:00', # Runtime in hrs
        memory = '10000', # RAM in MB
    threads: 10
    shell: """
        cp {input.txt} {output.txt}
        software/p2rank_2.4/prank predict -threads 10 -c alphafold -o {output.dir} {output.txt}
    """

rule af2_bulk:
    """
    time snakemake af2_bulk --cores 1 --until af2_wget --dry-run
    time snakemake af2_bulk --cores 1 --until af2_pdb --dry-run
    time snakemake af2_bulk --profile euler --cores 1 --dry-run
    """
    input:
        #[f'resources/afdb/{af2_bulk_id}.tar' for af2_bulk_id in config['af2_bulk_id'] ],
        #[f'resources/afdb/pdb/{af2_bulk_id}.txt' for af2_bulk_id in config['af2_bulk_id'] ],
        #[f'results_af2/af2_bulk_stats/{af2_bulk_id}.tsv' for af2_bulk_id in config['af2_bulk_id'] ],
        [f'results_af2/af2_bulk_prank/{af2_bulk_id}/' for af2_bulk_id in config['af2_bulk_id'] ],
