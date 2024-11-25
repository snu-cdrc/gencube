from tabulate import tabulate   # to print table form output

# Custom functions
from .utils import (
    check_argument,
    check_now,
    search_assembly, 
    json_to_dataframe,
    check_access_database,
    check_access_full_annotation,
    mkdir_raw_output,
    download_annotation,
    convert_chr_label_annotation,
    )
from .constants import (
    LS_GENCUBE_ANNOTATION_LABEL,
    )

## gencube annotation ------------------------------
def annotation (
    keywords,
    level,
    refseq,
    ucsc,
    latest,
    download,
    chr_style,
    recursive,
    ):
    # Check invalid arguments
    ls_level = ['complete', 'chromosome', 'scaffold', 'contig']
    ls_invalid_level = check_argument (level, ls_level, '-v/--level')
    
    ls_invalid_d = []
    if download:
        ls_download = ['gap', 'sr', 'td', 'wm', 'rmsk', 'cpg', 'gc']
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
        
        # Check_access to NCBI (Genbank & RefSeq) & Ensembl Rapid Release
        df_genome_plus, dic_genark_meta = check_access_database (df_genome, mode = 'annotation')
        
        print(tabulate(df_genome_plus[LS_GENCUBE_ANNOTATION_LABEL], headers='keys', tablefmt='grid'))
        print('')
        
        # Check full accessibility
        df_full_annotation = check_access_full_annotation (df_genome_plus, dic_genark_meta)
    
        # Save annotation files
        if download:
            mkdir_raw_output('annotation') # Make output folders
            dic_download = download_annotation(df_full_annotation, df_genome_plus, dic_genark_meta, download, recursive)
            # Change chromosome label style
            convert_chr_label_annotation (df_genome_plus, dic_download, chr_style, recursive)



