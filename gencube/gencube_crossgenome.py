from tabulate import tabulate   # to print table form output

# Custom functions
from .utils import (
    check_argument,
    check_now,
    search_assembly,
    json_to_dataframe,
    check_access_database,
    check_access_full_crossgenome,
    print_warning,
    mkdir_raw_output,
    save_metadata,
    download_crossgenome,
    )
from .constants import (
    LS_GENCUBE_CROSSGENOME_LABEL,
    )

## gencube crossgenome ------------------------------
def crossgenome (
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
        ls_download = ['ensembl_homology', 'toga_homology', 'toga_align_codon', 'toga_align_protein', 'toga_inact_mut']
        ls_invalid_d = check_argument (download, ls_download, '-d/--download')

    if len(ls_invalid_d) > 0 or len(ls_invalid_level) > 0:
        print('')
        return

    # Check the current time
    now = check_now ()

    # Search for keywords in the NCBI Assembly DB
    ls_search = search_assembly(keywords)

    if ls_search:
        # Convert JSON to DataFrame
        df_genome = json_to_dataframe (ls_search, level, refseq, ucsc, latest)

        # Check access to Ensembl Rapid Release & Zoonomia TOGA
        df_genome_plus, dic_ensembl_meta, df_zoonomia = check_access_database (df_genome, mode = 'crossgenome')

        print(tabulate(df_genome_plus[LS_GENCUBE_CROSSGENOME_LABEL], headers='keys', tablefmt='grid'))
        print('')
        # Print warning message
        print_warning(df_genome_plus, 100)

        # Check full accessibility
        df_full_crossgenome = check_access_full_crossgenome (df_genome_plus, dic_ensembl_meta, df_zoonomia)
        # Print warning message
        print_warning(df_genome_plus, 100)

        # Make output folders
        if metadata or download:
            mkdir_raw_output('crossgenome')

        # Save metadata
        if metadata:
            # Drop unnecessary columns
            df_genome_plus = df_genome_plus.drop(['Ensembl', 'Zoonomia'], axis=1)
            # Rename columns
            df_genome_plus[['Ensembl_db', 'Zoonomia_db']] = df_full_crossgenome[['Ensembl', 'Zoonomia']]
            # Save metadata
            save_metadata(df_genome_plus, 'crossgenome', keywords, level, now)

        # Save comparative genomics files
        if download:
            download_crossgenome(df_full_crossgenome, df_genome_plus, dic_ensembl_meta, df_zoonomia, download, recursive)
