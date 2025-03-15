# Custom functions
from .utils import (
    check_argument,
    check_now,
    make_query,
    search_sra,
    fetch_meta,
    convert_format,
    mkdir_raw_output,
    save_seq_metadata,
    get_fastq_dataframe,
    )

## gencube seqmeta ---------------------------------
def seqmeta(
    keywords,
    organism,
    strategy,
    source,
    platform,
    selection,
    filter,
    properties,
    layout,
    access,
    bioproject,
    biosample,
    accession,
    title,
    author,
    mbases,
    publication,
    modification,
    readlength,
    textword,
    exclude,
    detail,
    metadata,
    url,
    LS_STRATEGY,
    LS_SOURCE,
    LS_PLATFORM,
    LS_SELECTION,
    LS_FILTER,
    LS_PROPERTIES,
    ):
    # Check invalid arguments
    ls_invalid_strategy = check_argument (strategy, LS_STRATEGY, 'strategy')
    ls_invalid_source = check_argument (source, LS_SOURCE, 'source')
    ls_invalid_platform = check_argument (platform, LS_PLATFORM, 'platform')
    ls_invalid_selection = check_argument (selection, LS_SELECTION, 'selection')
    ls_invalid_filter = check_argument (filter, LS_FILTER, 'filter')
    ls_invalid_properties = check_argument (properties, LS_PROPERTIES, 'properties')

    if len(ls_invalid_strategy) > 0 or len(ls_invalid_source) > 0 or \
        len(ls_invalid_platform) > 0 or len(ls_invalid_selection) > 0 or \
        len(ls_invalid_filter) > 0 or len(ls_invalid_properties) > 0:
        print('')
        return

    # Check the current time
    now = check_now()

    # Input query
    query, dic_query_pars, dic_query_kwds, dic_query_excs = make_query(
        organism, strategy, source, platform, selection,
        filter, layout, access, bioproject, biosample, accession,
        title, author, publication, modification,
        properties, readlength, mbases, textword, keywords, exclude)

    print('# Search experimental sequencing data in NCBI SRA database')
    print(f'  Search query: {query} \n')

    # Check the number of searched results in each step
    if detail:
        print('# Check the number of searched result in each step')

        # Options
        if dic_query_pars:
            print('- Options')
            for key in list(dic_query_pars.keys()):
                out_count = search_sra(dic_query_pars[key])["Count"]
                if 'Intersection' not in key:
                    if dic_query_pars[key] == ' ':
                        print(f'  {key}')
                    else:
                        print(f'    {key}: {out_count}')
                else:
                    print(f'  Intersection (options): {out_count}')
            print('')

        # Keywords
        if keywords:
            print('- Keywords with options')
            for key in list(dic_query_kwds.keys()):
                if 'space' in key:
                    continue
                out_count = search_sra(dic_query_kwds[key])["Count"]
                print(f'  {key}: {out_count}')
            print('')

        # Keywords for exclusion
        if exclude:
            print('- Keywords for exclusion with options')
            for key in list(dic_query_excs.keys()):
                out_count = search_sra(dic_query_excs[key])["Count"]
                print(f'  {key}: {out_count}')
            print('')

    # Search query in the SRA database
    record = search_sra(query)
    search_ids = record['IdList']

    print('# Searched result')
    n = int(record["Count"])
    if n < 2:
        print(f'  Total {record["Count"]} experiment-level ID is searched.\n')
    else:
        print(f'  Total {record["Count"]} experiment-level IDs are searched.\n')

    # Fetch metadata, re-format, and save study-level and experiment-level tables
    if metadata:
        # Integrated metadata
        mkdir_raw_output('seqmeta') # Make output folders
        out_fetch = fetch_meta(search_ids)
        #out_fetch.to_csv('seqmeta_full.txt', sep='\t', header=True, index=False)
        #print('done')
        df_study, df_experiment = convert_format(out_fetch, query)
        out_name_url = save_seq_metadata(df_study, df_experiment, organism, now)

        if url:
            # Ouput file for download link
            get_fastq_dataframe(df_experiment, out_name_url)


    else:
        print('!! If you want to save the metadata of the searched datasets, please use the -m or --metadata option. \n')


