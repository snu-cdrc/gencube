from tabulate import tabulate   # to print table form output

# Custom functions
from .utils import (
    check_argument,
    check_now,
    search_assembly,
    json_to_dataframe,
    check_access_database,
    print_warning,
    mkdir_raw_output,
    save_metadata,
    download_genome,
    convert_chr_label_genome,
    )
from .constants import (
    LS_GENCUBE_GENOME_LABEL,
    )

## gencube genome ------------------------------
def genome (
    keywords,
    level,
    refseq,
    ucsc,
    latest,
    metadata,
    download,
    database,
    chr_style,
    masking,
    compresslevel,
    recursive,
    ):
    # Check invalid arguments
    ls_level = ['complete', 'chromosome', 'scaffold', 'contig']
    ls_invalid_level = check_argument (level, ls_level, '-v/--level')

    ls_invalid_db = []
    if download:
        ls_download = ['refseq', 'genbank', 'genark', 'ensembl']
        ls_invalid_db = check_argument (database, ls_download, '-db/--database')

    if len(ls_invalid_db) > 0 or len(ls_invalid_level) > 0:
        print('')
        return

    # Check the current time
    now = check_now()

    # Search for keywords in the NCBI Assembly DB
    ls_search = search_assembly(keywords)

    if ls_search:
        # Convert JSON to DataFrame
        df_genome = json_to_dataframe(ls_search, level, refseq, ucsc, latest)

        # Check access to NCBI (GenBank & RefSeq), GenArk & Ensembl Rapid Release
        df_genome_plus, dic_genark_meta, dic_ensembl_meta = check_access_database(df_genome, mode='genome')

        print(tabulate(df_genome_plus[LS_GENCUBE_GENOME_LABEL], headers='keys', tablefmt='grid'))
        print('')

        # Print warning message
        print_warning(df_genome_plus, 100)

        # Make output folders
        if metadata or download:
            mkdir_raw_output('genome')

        # Save metadata
        if metadata:
            save_metadata(df_genome_plus, 'genome', keywords, level, now)

        # Download genome files
        if download:

            dic_download = download_genome(df_genome_plus, database, dic_genark_meta, dic_ensembl_meta, recursive)
            # Change chromosome label style
            convert_chr_label_genome(df_genome_plus, dic_download, chr_style, masking, compresslevel, recursive)


