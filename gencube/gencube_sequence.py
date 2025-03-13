from tabulate import tabulate   # to print table form output

# Custom functions
from .utils import (
    check_argument,
    check_now,
    search_assembly,
    json_to_dataframe,
    check_access_database,
    check_access_full_sequence,
    print_warning,
    mkdir_raw_output,
    save_metadata,
    download_sequence,
    )
from .constants import (
    LS_GENCUBE_SEQUENCE_LABEL,
    )

## gencube sequence ------------------------------
def sequence (
    keywords,
    level,
    refseq,
    ucsc,
    latest,
    metadata,
    download,
    recursive,
    ):
    # Check invalid arguments
    ls_level = ['complete', 'chromosome', 'scaffold', 'contig']
    ls_invalid_level = check_argument (level, ls_level, '-v/--level')

    ls_invalid_d = []
    if download:
        ls_download = ['refseq_rna', 'refseq_rna_genomic', 'refseq_cds_genomic', 'refseq_pseudo', 'ensembl_cdna', 'refseq_pep', 'refseq_pep_cds', 'ensembl_pep']
        ls_invalid_d = check_argument (download, ls_download, '-d/--download')

    if len(ls_invalid_d) > 0 or len(ls_invalid_level) > 0:
        print('')
        return

    # Check the current time
    now = check_now ()

    # Search keywords in NCBI Assembly DB
    ls_search = search_assembly(keywords)

    if ls_search:
        # Convert JSON to DataFrame
        df_genome = json_to_dataframe (ls_search, level, refseq, ucsc, latest)

        # Check access to NCBI RefSeq & Ensembl Rapid Release
        df_genome_plus, dic_ensembl_meta = check_access_database (df_genome, mode = 'sequence')

        print(tabulate(df_genome_plus[LS_GENCUBE_SEQUENCE_LABEL], headers='keys', tablefmt='grid'))
        print('')
        # Print warning message
        print_warning(df_genome_plus, 100)

        # Check full accessibility
        df_full_sequence = check_access_full_sequence (df_genome_plus, dic_ensembl_meta)
        # Print warning message
        print_warning(df_genome_plus, 100)

        # Make output folders
        if metadata or download:
            mkdir_raw_output('sequence')

        # Save metadata
        if metadata:
            # Drop unnecessary columns
            df_genome_plus = df_genome_plus.drop(['Ensembl'], axis=1)
            # Rename columns
            df_genome_plus[['RefSeq_db', 'Ensembl_db']] = df_full_sequence[['RefSeq', 'Ensembl']]
            # Save metadata
            save_metadata(df_genome_plus, 'sequence', keywords, level, now)

        # Save sequence files
        if download:
            download_sequence(df_full_sequence, df_genome_plus, dic_ensembl_meta, download, recursive)
