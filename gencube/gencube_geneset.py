from tabulate import tabulate   # to print table form output

# Custom functions
from .utils import (
    check_argument,
    check_now,
    search_assembly,
    json_to_dataframe,
    check_access_database,
    check_access_full_geneset,
    print_warning,
    mkdir_raw_output,
    save_metadata,
    download_geneset,
    convert_chr_label_geneset,
    #save_pickle,
    )
from .constants import (
    LS_GENCUBE_GENSET_LABEL,
    )

## gencube geneset ------------------------------
def geneset (
    keywords,
    level,
    refseq,
    ucsc,
    latest,
    metadata,
    download,
    chr_style,
    recursive,
    ):
    # Check invalid arguments
    ls_level = ['complete', 'chromosome', 'scaffold', 'contig']
    ls_invalid_level = check_argument (level, ls_level, '-v/--level')

    ls_invalid_d = []
    if download:
        ls_download = ['refseq_gtf', 'refseq_gff', 'gnomon', 'cross', 'same', 'augustus', 'xenoref', 'genark_ref', 'ensembl_gtf', 'ensembl_gff', 'toga_gtf', 'toga_bed', 'toga_pseudo']
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

        # Check access to NCBI RefSeq, GenArk, Ensembl Rapid Release & Zoonomia TOGA
        df_genome_plus, dic_genark_meta, dic_ensembl_meta, df_zoonomia = check_access_database (df_genome, mode = 'geneset')

        print(tabulate(df_genome_plus[LS_GENCUBE_GENSET_LABEL], headers='keys', tablefmt='grid'))
        print('')
        # Print warning message
        print_warning(df_genome_plus, 100)

        # Check full accessibility
        df_full_geneset = check_access_full_geneset (df_genome_plus, dic_genark_meta, dic_ensembl_meta, df_zoonomia)
        # Print warning message
        print_warning(df_genome_plus, 100)

        # Make output folders
        if metadata or download:
            mkdir_raw_output('geneset')

        # Save metadata
        if metadata:
            # Drop unnecessary columns
            df_genome_plus = df_genome_plus.drop(['GenArk', 'Ensembl', 'Zoonomia'], axis=1)
            # Rename columns
            df_genome_plus[['RefSeq_db', 'GenArk_db', 'Ensembl_db', 'Zoonomia_db']] = df_full_geneset[['RefSeq', 'GenArk', 'Ensembl', 'Zoonomia']]
            # Save metadata
            save_metadata(df_genome_plus, 'geneset', keywords, level, now)

        # Save geneset files
        if download:
            dic_download = download_geneset(df_full_geneset, df_genome_plus, dic_ensembl_meta, dic_genark_meta, df_zoonomia, download, recursive)
            # Change chromosome label style
            convert_chr_label_geneset (df_genome_plus, dic_download, chr_style, recursive)
