"""
Extract subsets from UniProt .xml bulk downloads (https://www.uniprot.org/help/downloads)

Refs:
https://riptutorial.com/python/example/25995/opening-and-reading-large-xml-files-using-iterparse--incremental-parsing-
https://stackoverflow.com/questions/10074200/should-memory-usage-increase-when-using-elementtree-iterparse-when-clearing

Downloading:
$ mamba create -n aria2c -c bioconda aria2
$ conda activate aria2c
$ aria2c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/RELEASE.metalink
"""

def extract_n_entries(fp, stop=1000):
    events_ = ('start', 'end')
    path = []
    n_entries = 0
    with gzip.open(fp, 'r') as fh:
        for event, elem in itertools.islice(ET.iterparse(fh, events=events_), stop):
            _, _, elem_tag = elem.tag.rpartition('}')
            if event == 'start':
                path.append(elem_tag)
            elif event == 'end':
                if path == ['uniprot', 'entry']:
                    n_entries += 1
                    elem.clear()
                path.pop()
    return n_entries

rule extract_n_entries:
    """
    time snakemake extract_n_entries --cores 1 --dry-run
    time snakemake extract_n_entries --profile euler --cores 1 --dry-run
    """
    input:
        sprot_xml_gz = 'current_release/uniprot_sprot.xml.gz', # 0.8GB
        trembl_xml_gz = 'current_release/uniprot_trembl.xml.gz', # 190GB
    output:
        txt = 'results/extract_n_entries.txt',
    resources:
        runtime = '36:00', # Runtime in hrs
        memory = '4096', # RAM in MB
    run:
        with open(output.txt, 'w') as fh:
            stop = 100000 #None
            print(f'{extract_n_entries(fp=input.sprot_xml_gz, stop=stop):,} records in {input.sprot_xml_gz} (stop={stop})', file=fh) # 5min
            stop = 1000000
            print(f'{extract_n_entries(fp=input.trembl_xml_gz, stop=stop):,} records in {input.trembl_xml_gz} (stop={stop})', file=fh)
        # Estimated runtime for uniprot_trembl.xml.gz:
        # 5 / .8 * 190 / 60
        # 19.791666666666668

def extract_all_vs_all(fp, accessions, stop=1000):
    events_ = ('start', 'end')
    path = []
    path_elem = []
    attr = {}
    with gzip.open(fp, 'r') as fh:
        for event, elem in itertools.islice(ET.iterparse(fh, events=events_), stop):
            _, _, elem_tag = elem.tag.rpartition('}')
            if event == 'start':
                path.append(elem_tag)
                path_elem.append(elem)
            elif event == 'end':
                if path == ['uniprot', 'entry', 'accession']:
                    attr['accession'] = set([elem.text]) if not('accession' in attr.keys()) else attr['accession'] | set([elem.text])
                elif path == ['uniprot', 'entry', 'gene', 'name'] and elem.attrib['type'] == 'primary':
                    attr['gene_name_primary'] = elem.text
                elif path == ['uniprot', 'entry', 'organism', 'name'] and (elem.attrib['type'] == 'scientific'):
                    attr['organism_name_scientific'] = elem.text
                elif path == ['uniprot', 'entry', 'organism', 'dbReference'] and (elem.attrib['type'] == 'NCBI Taxonomy'):
                    attr['organism_dbReference_NCBI_Taxonomy'] = elem.attrib['id']
                elif path == ['uniprot', 'entry', 'dbReference'] and (elem.attrib['type'] == 'EC'):
                    attr['dbReference_EC_Id'] = elem.attrib['id']
                elif path == ['uniprot', 'entry', 'dbReference'] and (elem.attrib['type'] == 'PDB'):
                    attr['dbReference_PDB_Id'] = elem.attrib['id'] if not('dbReference_PDB' in attr.keys()) else ';'.join([attr['dbReference_PDB'], elem.attrib['id']])
                elif path == ['uniprot', 'entry', 'dbReference'] and (elem.attrib['type'] == 'Pfam'):
                    attr['dbReference_Pfam_Id'] = elem.attrib['id']
                elif path == ['uniprot', 'entry', 'dbReference', 'property'] and (path_elem[2].attrib['type'] == 'Pfam') and (elem.attrib['type'] == 'entry name'):
                    attr['dbReference_Pfam_entry_name'] = elem.attrib['value']
                elif path == ['uniprot', 'entry', 'dbReference'] and (elem.attrib['type'] == 'InterPro'):
                    attr['dbReference_InterPro_Id'] = elem.attrib['id']
                elif path == ['uniprot', 'entry', 'dbReference', 'property'] and (path_elem[2].attrib['type'] == 'InterPro') and (elem.attrib['type'] == 'entry name'):
                    attr['dbReference_InterPro_entry_name'] = elem.attrib['value']
                elif path == ['uniprot', 'entry', 'dbReference'] and (elem.attrib['type'] == 'PROSITE'):
                    attr['dbReference_PROSITE_Id'] = elem.attrib['id']
                elif path == ['uniprot', 'entry', 'dbReference', 'property'] and (path_elem[2].attrib['type'] == 'PROSITE') and (elem.attrib['type'] == 'entry name'):
                    attr['dbReference_PROSITE_entry_name'] = elem.attrib['value']
                elif path == ['uniprot', 'entry', 'dbReference'] and (elem.attrib['type'] == 'SUPFAM'):
                    attr['dbReference_SUPFAM_Id'] = elem.attrib['id']
                elif path == ['uniprot', 'entry', 'dbReference', 'property'] and (path_elem[2].attrib['type'] == 'SUPFAM') and (elem.attrib['type'] == 'entry name'):
                    attr['dbReference_SUPFAM_entry_name'] = elem.attrib['value']
                elif path == ['uniprot', 'entry']:
                    if len(attr['accession'] & accessions) > 0:
                        attr['accession'] = ';'.join(sorted(attr['accession'] & accessions))
                        yield(attr)
                    attr = {}
                path.pop()
                path_elem.pop()
                elem.clear()

rule extract_all_vs_all:
    """
    time snakemake extract_all_vs_all --cores 1 --dry-run
    time snakemake --profile euler extract_all_vs_all --cores 1 --dry-run
    """
    input:
        sprot_xml_gz = 'current_release/uniprot_sprot.xml.gz',
        trembl_xml_gz = 'current_release/uniprot_trembl.xml.gz',
        tsv_gz = '../../22.07.01_foldseek_all_vs_all/afswiss_proteom_ava.accessions.m8.gz',
    output:
        tsv_gz = '../../22.07.01_foldseek_all_vs_all/afswiss_proteom_ava.accessions_with_metadata-YY.MM.DD.tsv.gz',
    resources:
        runtime = '48:00', # Runtime in hrs
        memory = '100000', # RAM in MB
    run:
        #stop = 1000000
        stop = None
        # Note: some entries have *multiple* accession numbers!
        #accessions = set(['P0C9F2', 'P0C9F6', 'A0A816AIK7', 'A0A8J4G4B6', 'Q1HR36', 'P0C9F0', 'P31937', 'Q9UDN3', 'P63103',
        #                  'P29213', 'Q3ZCF9', # P29213 is an alias of Q3ZCF9; currently not captured
        #                  'P0C9H3', #dbReference_Pfam_Id, dbReference_Pfam_entry name
        #])
        accessions = set(pd.read_csv(input.tsv_gz, names=['acc'])['acc'])

        cols_ = [
            'accession',
            'gene_name_primary',
            'organism_name_scientific',
            'organism_dbReference_NCBI_Taxonomy',
            'dbReference_EC_Id',
            'dbReference_PDB_Id',
            'dbReference_Pfam_Id',
            'dbReference_Pfam_entry_name',
            'dbReference_InterPro_Id',
            'dbReference_InterPro_entry_name',
            'dbReference_PROSITE_Id',
            'dbReference_PROSITE_entry_name',
            'dbReference_SUPFAM_Id',
            'dbReference_SUPFAM_entry_name',
        ]
        df_tsv = pd.DataFrame.from_records(itertools.chain(
                extract_all_vs_all(fp=input.sprot_xml_gz, accessions=accessions, stop=stop),
                extract_all_vs_all(fp=input.trembl_xml_gz, accessions=accessions, stop=stop),
        ), index='accession', columns=cols_)

        #print(df_tsv.loc[['Q1HR36', 'P0C9H3']].head(5).transpose())
        df_tsv.to_csv(output.tsv_gz, index=True, header=True, sep='\t')
