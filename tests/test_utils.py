import pytest
from unittest.mock import patch, mock_open, Mock, MagicMock, call
import pandas as pd
import pickle
import os
from datetime import datetime
import requests
from bs4 import BeautifulSoup
from io import StringIO
import hashlib
from Bio import Entrez
from pathlib import Path
import urllib3                # For advanced HTTP client functionalities


from gencube.utils import (
    get_entrez_info,
    save_pickle,
    load_pickle,
    check_now,
    add_string,
    check_argument,
    check_url,
    list_ftp_directory,
    list_http_folders,
    list_http_files,
    download_csv,
    mkdir_raw_output,
    save_metadata,
    calculate_md5,
    download_genome_url,
    download_url,
    get_md5,
    delete_file,
    make_chr_dataframe,
    search_assembly,
    json_to_dataframe,
    check_access_database,
    
    
    DOWNLOAD_FOLDER_NAME,
    OUT_FOLDER_NAME,
    script_dir,
    LS_ASSEMBLY_REPORT_LABEL,
    

)

# Get email and api key information for E-utility
get_entrez_info()

## -----------------------------------------------------------
## gencube all subcommands
## -----------------------------------------------------------
@patch('builtins.input', side_effect=['test@example.com', 'test_api_key'])
@patch('builtins.open', new_callable=mock_open)
@patch('pathlib.Path.exists', return_value=False)
@patch('pathlib.Path.home', return_value=Path('/mock/home'))
def test_get_entrez_info_new_file(mock_home, mock_exists, mock_open, mock_input):
    get_entrez_info()

    mock_open.assert_called_once_with(Path('/mock/home/.gencube_entrez_info'), 'w')
    handle = mock_open()
    handle.write.assert_any_call("# This information is used in E-utilities\n")
    handle.write.assert_any_call("email = test@example.com\n")
    handle.write.assert_any_call("api_key = test_api_key\n")

    assert Entrez.email == 'test@example.com'
    assert Entrez.api_key == 'test_api_key'

@patch('builtins.open', new_callable=mock_open, read_data="# This information is used in E-utilities.\nemail = test@example.com\napi_key = test_api_key\n")
@patch('pathlib.Path.exists', return_value=True)
@patch('pathlib.Path.home', return_value=Path('/mock/home'))
def test_get_entrez_info_existing_file(mock_home, mock_exists, mock_open):
    get_entrez_info()

    mock_open.assert_called_once_with(Path('/mock/home/.gencube_entrez_info'), 'r')

    assert Entrez.email == 'test@example.com'
    assert Entrez.api_key == 'test_api_key'
    
## -----------------------------------------------------------
## gencube genome, geneset, sequence, annotation, crossgenome
## -----------------------------------------------------------
@patch('builtins.open', new_callable=mock_open)
@patch('pandas.DataFrame.to_pickle')
def test_save_pickle_dataframe(mock_to_pickle, mock_open):
    df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
    save_pickle(df, 'test_df')

    mock_to_pickle.assert_called_once_with('tests/data/test_df.pkl')
    mock_open.assert_not_called()  # open is not called for DataFrame

@patch('builtins.open', new_callable=mock_open)
@patch('pickle.dump')
def test_save_pickle_other(mock_pickle_dump, mock_open):
    data = {'key': 'value'}
    save_pickle(data, 'test_dict')

    mock_open.assert_called_once_with('tests/data/test_dict.pkl', 'wb')
    mock_pickle_dump.assert_called_once()

@patch('builtins.open', new_callable=mock_open)
@patch('pandas.read_pickle')
def test_load_pickle_dataframe(mock_read_pickle, mock_open):
    mock_read_pickle.return_value = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
    result = load_pickle('test_df', type='dataframe')

    mock_read_pickle.assert_called_once_with('tests/data/test_df.pkl')
    mock_open.assert_not_called()  # open is not called for DataFrame
    assert isinstance(result, pd.DataFrame)

@patch('builtins.open', new_callable=mock_open)
@patch('pickle.load')
def test_load_pickle_other(mock_pickle_load, mock_open):
    mock_pickle_load.return_value = {'key': 'value'}
    result = load_pickle('test_dict')

    mock_open.assert_called_once_with('tests/data/test_dict.pkl', 'rb')
    mock_pickle_load.assert_called_once()
    assert isinstance(result, dict)
    assert result == {'key': 'value'}
    
@patch('gencube.utils.datetime')
def test_check_now(mock_datetime):
    # Mock datetime to return a fixed date and time
    mock_datetime.now.return_value = datetime(2024, 7, 15, 12, 34, 56)
    mock_datetime.strftime = datetime.strftime  # Use the actual strftime method

    result = check_now()
    expected_result = "240715_123456"
    assert result == expected_result

@pytest.mark.parametrize("pre, string, expected", [
    ('Hello', 'world', 'Hello, world'),
    ('', 'world', 'world'),
    (None, 'world', 'world'),
    ('', '', ''),
    ('Hello', '', 'Hello, ')
])
def test_add_string(pre, string, expected):
    result = add_string(pre, string)
    assert result == expected

@pytest.mark.parametrize("input, ls, str, expected_result, expected_print", [
    ('a,b,x', ['a', 'b', 'c'], 'test', ['x'], 'Invalid argument test: x'),
    ('a,b,x,y', ['a', 'b', 'c'], 'test', ['x', 'y'], 'Invalid argument test: x, y'),
    ('', ['a', 'b', 'c'], 'test', [], None),
    (None, ['a', 'b', 'c'], 'test', [], None),
])
@patch('builtins.print')
def test_check_argument(mock_print, input, ls, str, expected_result, expected_print):
    result = check_argument(input, ls, str)
    assert result == expected_result
    if expected_print:
        mock_print.assert_called_once_with(expected_print)
    else:
        mock_print.assert_not_called()

@pytest.mark.parametrize("verify, show_output, file_name, mock_status_code, mock_headers, expected_result, expected_print", [
    (True, True, '', 200, {}, True, None),
    (True, True, '', 301, {'Location': 'http://example.com/new_location'}, True, "  file: !! Redirected to http://example.com/new_location"),
    (True, True, '', 403, {}, False, "  file: !! Access denied. File cannot be downloaded"),
    (True, True, '', 404, {}, False, "  file: !! File not found"),
    (True, True, '', 429, {}, False, "  file: !! Too many requests. Please try again later"),
    (True, True, '', 500, {}, False, "  file: !! Server error. Please try again later"),
    (True, True, '', 418, {}, False, "  file: !! Unexpected status code: 418"),
    (True, False, '', 404, {}, False, None),
])
@patch('requests.head')
@patch('builtins.print')
def test_check_url(mock_print, mock_requests_head, verify, show_output, file_name, mock_status_code, mock_headers, expected_result, expected_print):
    mock_response = Mock()
    mock_response.status_code = mock_status_code
    mock_response.headers = mock_headers
    mock_requests_head.return_value = mock_response

    url = "http://example.com/file"
    result = check_url(url, verify=verify, show_output=show_output, file_name=file_name)
    
    assert result == expected_result
    
    if expected_print:
        mock_print.assert_called_once_with(expected_print)
    else:
        mock_print.assert_not_called()

@patch('ftplib.FTP')
def test_list_ftp_directory(mock_ftp_cls):
    # Mock FTP instance and its methods
    mock_ftp = MagicMock()
    mock_ftp_cls.return_value = mock_ftp

    # Define the behavior of the mock methods
    mock_ftp.nlst.return_value = ['file1.txt', 'file2.txt', 'file3.txt']

    host = 'ftp.example.com'
    directory = '/example_directory'
    expected_files = ['file1.txt', 'file2.txt', 'file3.txt']

    result = list_ftp_directory(host, directory)

    # Assertions
    mock_ftp_cls.assert_called_once_with(host)
    mock_ftp.login.assert_called_once()
    mock_ftp.cwd.assert_called_once_with(directory)
    mock_ftp.nlst.assert_called_once()
    mock_ftp.quit.assert_called_once()

    assert result == expected_files

@patch('requests.get')
@patch('builtins.print')
def test_list_http_folders_success(mock_print, mock_requests_get):
    # Mock the response from requests.get
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.content = '''
    <html>
        <body>
            <a href="folder1/">folder1/</a>
            <a href="folder2/">folder2/</a>
            <a href="file.txt">file.txt</a>
        </body>
    </html>
    '''
    mock_requests_get.return_value = mock_response

    url = "http://example.com"
    expected_folders = ['folder1', 'folder2']

    result = list_http_folders(url)

    assert result == expected_folders
    mock_requests_get.assert_called_once_with(url, verify=False)
    mock_print.assert_not_called()

@patch('requests.get')
@patch('builtins.print')
def test_list_http_folders_request_exception(mock_print, mock_requests_get):
    # Mock requests.get to raise a RequestException
    mock_requests_get.side_effect = requests.exceptions.RequestException("Test exception")

    url = "http://example.com"
    expected_folders = []

    result = list_http_folders(url)

    assert result == expected_folders
    mock_requests_get.assert_called_once_with(url, verify=False)
    mock_print.assert_called_once_with("Error accessing URL: Test exception")

@patch('requests.get')
@patch('builtins.print')
def test_list_http_files_success(mock_print, mock_requests_get):
    # Mock the response from requests.get
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.content = '''
    <html>
        <body>
            <a href="file1.txt">file1.txt</a>
            <a href="file2.txt">file2.txt</a>
            <a href="folder/">folder/</a>
        </body>
    </html>
    '''
    mock_requests_get.return_value = mock_response

    url = "http://example.com"
    expected_files = ['file1.txt', 'file2.txt']

    result = list_http_files(url)

    assert result == expected_files
    mock_requests_get.assert_called_once_with(url, verify=False)
    mock_print.assert_not_called()

@patch('requests.get')
@patch('builtins.print')
def test_list_http_files_request_exception(mock_print, mock_requests_get):
    # Mock requests.get to raise a RequestException
    mock_requests_get.side_effect = requests.exceptions.RequestException("Test exception")

    url = "http://example.com"
    expected_files = []

    result = list_http_files(url)

    assert result == expected_files
    mock_requests_get.assert_called_once_with(url, verify=False)
    mock_print.assert_called_once_with("Error accessing URL: Test exception")

@patch('requests.get')
def test_download_csv_success(mock_requests_get):
    # Mock the response from requests.get
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.text = "col1,col2,col3\nval1,val2,val3"
    mock_requests_get.return_value = mock_response

    url = "http://example.com/file.csv"
    expected_content = "col1,col2,col3\nval1,val2,val3"

    result = download_csv(url)

    assert isinstance(result, StringIO)
    assert result.getvalue() == expected_content
    mock_requests_get.assert_called_once_with(url, verify=True)

@patch('requests.get')
def test_download_csv_request_exception(mock_requests_get):
    # Mock requests.get to raise a RequestException
    mock_requests_get.side_effect = requests.exceptions.RequestException("Test exception")

    url = "http://example.com/file.csv"

    with pytest.raises(requests.exceptions.RequestException):
        download_csv(url)

    mock_requests_get.assert_called_once_with(url, verify=True)

@patch('os.path.exists')
@patch('os.mkdir')
def test_mkdir_raw_output(mock_mkdir, mock_exists):
    # Setup the mock return values
    mock_exists.side_effect = lambda x: x in (DOWNLOAD_FOLDER_NAME, OUT_FOLDER_NAME)
    
    # Call the function
    mkdir_raw_output()
    
    # Check that os.path.exists was called correctly
    mock_exists.assert_any_call(DOWNLOAD_FOLDER_NAME)
    mock_exists.assert_any_call(OUT_FOLDER_NAME)
    
    # Check that os.mkdir was not called because folders already exist
    mock_mkdir.assert_not_called()
    
    # Change the mock return values to simulate non-existing directories
    mock_exists.side_effect = lambda x: False
    
    # Call the function again
    mkdir_raw_output()
    
    # Check that os.path.exists was called correctly
    mock_exists.assert_any_call(DOWNLOAD_FOLDER_NAME)
    mock_exists.assert_any_call(OUT_FOLDER_NAME)
    
    # Check that os.mkdir was called twice (once for each directory)
    mock_mkdir.assert_any_call(DOWNLOAD_FOLDER_NAME)
    mock_mkdir.assert_any_call(OUT_FOLDER_NAME)
    assert mock_mkdir.call_count == 2
"""
@patch('pandas.DataFrame.to_csv')
@patch('builtins.print')
def test_save_metadata(mock_print, mock_to_csv):
    # Create a mock dataframe
    data = {'A': [1, 2], 'B': [3, 4], 'NCBI': [5, 6]}
    df = pd.DataFrame(data)
    
    function = 'test_function'
    keywords = ['keyword1', 'keyword2']
    level = 'test_level'
    now = '20220101_123456'
    
    save_metadata(df, function, keywords, level, now)
    
    # Expected metadata dataframe (without 'NCBI' column)
    expected_df_meta = df.drop(columns=['NCBI'])
    
    # Expected file name
    expected_out_name = f"Meta_{function}_2keywords_{level.replace(',', '-')}_{now}.txt"
    expected_file_path = f'{OUT_FOLDER_NAME}/{expected_out_name}'
    
    # Assert to_csv was called correctly
    mock_to_csv.assert_called_once_with(expected_file_path, sep='\t', index=False)
    
    # Assert print statements
    expected_print_calls = [
        '# Metadata are saved:',
        f'  {expected_out_name}\n'
    ]
    mock_print.assert_has_calls([patch('builtins.print', call) for call in expected_print_calls], any_order=False)
"""
@patch('builtins.open', new_callable=mock_open, read_data=b'This is a test file.')
def test_calculate_md5(mock_file):
    # Expected MD5 hash for the given test file content
    expected_md5 = hashlib.md5(b'This is a test file.').hexdigest()

    filename = 'testfile.txt'
    result = calculate_md5(filename)

    assert result == expected_md5
    mock_file.assert_called_once_with(filename, 'rb')

@patch('requests.get')
def test_get_md5(mock_get):
    # Define the mock response data
    md5_data = "d41d8cd98f00b204e9800998ecf8427e  file1.txt\n" \
               "0cc175b9c0f1b6a831c399e269772661  file2.txt\n"

    # Create a mock response object with the specified text
    mock_response = Mock()
    mock_response.text = md5_data
    mock_response.status_code = 200
    mock_get.return_value = mock_response

    url_md5sum = 'http://example.com/md5sum'

    # Call the function
    df = get_md5(url_md5sum)

    # Create the expected DataFrame
    expected_df = pd.DataFrame({
        'MD5': ['d41d8cd98f00b204e9800998ecf8427e', '0cc175b9c0f1b6a831c399e269772661'],
        'Filename': ['file1.txt', 'file2.txt']
    })

    # Assert that the returned DataFrame is equal to the expected DataFrame
    pd.testing.assert_frame_equal(df, expected_df)

@pytest.mark.parametrize("exists, remove_side_effect, expected_print_calls, expected_remove_calls", [
    (True, None, [], [call('testfile.txt')]),  # File exists, no exception
    (False, None, [], []),  # File does not exist
    (True, Exception("Test exception"), [call("Error deleting testfile.txt: Test exception")], [call('testfile.txt')])  # Exception raised
])
@patch('os.path.exists')
@patch('os.remove')
@patch('builtins.print')
def test_delete_file(mock_print, mock_remove, mock_exists, exists, remove_side_effect, expected_print_calls, expected_remove_calls):
    # Set the return value of os.path.exists
    mock_exists.return_value = exists
    
    # Set the side effect of os.remove
    if remove_side_effect:
        mock_remove.side_effect = remove_side_effect
    else:
        mock_remove.side_effect = None

    file_path = 'testfile.txt'
    delete_file(file_path)

    # Assert os.path.exists was called with the correct file path
    mock_exists.assert_called_once_with(file_path)
    # Assert os.remove was called with the correct file path if expected
    if expected_remove_calls:
        mock_remove.assert_has_calls(expected_remove_calls, any_order=False)
    else:
        mock_remove.assert_not_called()
    # Assert print was called with the correct message if expected
    if expected_print_calls:
        mock_print.assert_has_calls(expected_print_calls, any_order=False)
    else:
        mock_print.assert_not_called()

report_content = """
# Assembly name:  Dog10K_Boxer_Tasha
# Organism name:  Canis lupus familiaris (dog)
# Infraspecific name:  breed=boxer
# Isolate:  Tasha
# Sex:  female
# Taxid:          9615
# BioSample:      SAMN02953603
# BioProject:     PRJNA13179
# Submitter:      Dog Genome Sequencing Consortium
# Date:           2020-10-06
# Synonyms:       canFam6	
# Assembly type:  haploid
# Release type:   major
# Assembly level: Chromosome
# Genome representation: full
# WGS project:    AAEX04
# Assembly method: Canu v. 1.8
# Expected final version: yes
# Genome coverage: 100.0x
# Sequencing technology: PacBio Sequel
# GenBank assembly accession: GCA_000002285.4
# RefSeq assembly accession: GCF_000002285.5
# RefSeq assembly and GenBank assemblies identical: no
#
## Assembly-Units:
## GenBank Unit Accession	RefSeq Unit Accession	Assembly-Unit name
## GCA_000000145.4	GCF_000000145.4	Primary Assembly
## GCA_013123455.1	GCF_000184225.1	non-nuclear
#
# Ordered by chromosome/plasmid; the chromosomes/plasmids are followed by
# unlocalized scaffolds.
# Unplaced scaffolds are listed at the end.
# RefSeq is equal or derived from GenBank object.
#
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000001.4	=	NC_006583.4	Primary Assembly	122014068	chr1
2	assembled-molecule	2	Chromosome	CM000002.4	=	NC_006584.4	Primary Assembly	82037489	chr2
chrUn6	unplaced-scaffold	na	na	AAEX04000039.1	=	NW_023329752.1	Primary Assembly	551	chrUn_NW_023329752v1
MT	assembled-molecule	MT	Mitochondrion	CM023446.1	<>	na	non-nuclear	16735	na
MT	assembled-molecule	MT	Mitochondrion	na	<>	NC_002008.4	non-nuclear	16727	chrM
"""
expected_output = pd.DataFrame({
    'GenBank': ['CM000001.4', 'CM000002.4', 'AAEX04000039.1', 'CM023446.1'],
    'RefSeq': ['NC_006583.4', 'NC_006584.4', 'NW_023329752.1', 'NC_002008.4'],
    'UCSC': ['chr1', 'chr2', 'chrUn_NW_023329752v1', 'chrM'],
    'Ensembl': ['1', '2', 'AAEX04000039.1', 'MT'],
    'Gencode': ['chr1', 'chr2', 'AAEX04000039.1', 'chrM']
})
@patch('os.path.exists', return_value=True)
@patch('builtins.open', new_callable=mock_open, read_data=report_content)
def test_make_chr_dataframe(mock_file, mock_exists):
    organism = 'dog'
    assembly_id = 'GCA_000002285.4'
    
    # Call the function
    result = make_chr_dataframe(organism, assembly_id)

    # Assert that the result matches the expected output
    pd.testing.assert_frame_equal(result.reset_index(drop=True), expected_output)


def test_search_assembly():
    result = search_assembly(['GCF_000002285.5'])
    result_null = search_assembly(['GCF_XXXXXXXX.X'])
    
    expected_result = load_pickle('utils_search_assembly', 'dataframe')
    
    assert result == expected_result
    assert isinstance(result_null, type(None))
    
def test_json_to_dataframe():
    input = load_pickle('utils_search_assembly', 'dataframe')
    level = 'complete,chromosome'
    refseq = False
    ucsc = False
    latest = False
    
    result = json_to_dataframe (input, level, refseq, ucsc, latest)
    expected_result = load_pickle('utils_json_to_dataframe', 'dataframe')
    
    print(result)
    print(expected_result)
    
    assert result.equals(expected_result)
    
def test_check_access_database():
    # Count return values
    def count_return_values (return_value):
        return len(return_value) if isinstance(return_value, tuple) else 1
        
    intput = load_pickle('utils_json_to_dataframe', 'dataframe')
    
    tp_values = check_access_database(intput, mode = 'genome')
    assert count_return_values(tp_values) == 3
    assert isinstance(tp_values[0], pd.DataFrame)
    assert isinstance(tp_values[1], dict)
    assert isinstance(tp_values[2], dict)
    
    tp_values = check_access_database(intput, mode = 'geneset')
    assert count_return_values(tp_values) == 4
    assert isinstance(tp_values[0], pd.DataFrame)
    assert isinstance(tp_values[1], dict)
    assert isinstance(tp_values[2], dict)
    assert isinstance(tp_values[3], pd.DataFrame)
    
    tp_values = check_access_database(intput, mode = 'sequence')
    assert count_return_values(tp_values) == 2
    assert isinstance(tp_values[0], pd.DataFrame)
    assert isinstance(tp_values[1], dict)
    
    tp_values = check_access_database(intput, mode = 'annotation')
    assert count_return_values(tp_values) == 2
    assert isinstance(tp_values[0], pd.DataFrame)
    assert isinstance(tp_values[1], dict)
    
    tp_values = check_access_database(intput, mode = 'crossgenome')
    assert count_return_values(tp_values) == 3
    assert isinstance(tp_values[0], pd.DataFrame)
    assert isinstance(tp_values[1], dict)
    assert isinstance(tp_values[2], pd.DataFrame)
