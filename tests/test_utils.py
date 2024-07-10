import pytest
from Bio import Entrez        # For accessing NCBI Entrez databases
import pandas as pd           # For data manipulation and analysis
import numpy as np            # For numerical operations
import requests               # For making HTTP requests
import ftplib                 # For handling FTP connectionssave_metadata
import urllib3                # For advanced HTTP client functionalities
import xmltodict              # For converting XML data to Python dictionaries
import json                   # For parsing JSON data
from bs4 import BeautifulSoup # For parsing HTML and XML documents
from io import StringIO       # For in-memory file-like objects
import os                     # For interacting with the operating system
import re                     # For regular expression operations
import gzip                   # For reading and writing gzip-compressed files
import hashlib                # For generating MD5 checksums
from tqdm import tqdm         # For displaying progress bars
from tabulate import tabulate # For printing data in a tabular format
from datetime import datetime # For manipulating date and time objects
import time                   # For handling time-related tasks
import pickle

import requests_mock

from unittest.mock import patch, MagicMock
from gencube.utils import check_url  # 적절한 모듈 이름으로 대체하세요
# To suppress the SSL warnings when using verify=False in the requests library
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning) 

from gencube.utils import (
    load_pickle,
    check_now,
    search_assembly,

    check_access_database,
    json_to_dataframe,
    download_genome,
    convert_chr_label_genome,
    check_access_full_geneset,
    download_geneset,
    convert_chr_label_geneset,
    check_access_full_sequence,
    download_sequence,
    check_access_full_annotation,
    download_annotation,
    check_access_full_crossgenome,
    download_crossgenome,
    make_chr_dataframe,
    make_query,
    add_operator,
    search_sra,
    fetch_meta,
    convert_format,
    save_seq_metadata,
    save_metadata,
    
    
    check_url,
    list_ftp_directory,
    list_http_folders,
    download_csv,
    download_genome_url,
    download_url,
    get_md5,
    calculate_md5,
    delete_file,
    
    
    DOWNLOAD_FOLDER_NAME,
)

## -----------------------------------------------------------
## gencube genome, geneset, sequence, annotation, crossgenome
## -----------------------------------------------------------
def test_check_now():
    now = check_now()
    assert isinstance(now, str)
    assert len(now) > 0

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

## ---------------------------------------------
## gencube genome
## ---------------------------------------------
def test_download_genome():
    print('')

def test_convert_chr_label_genome():
    print('')

## ---------------------------------------------
## gencube geneset
## ---------------------------------------------
def test_check_access_full_geneset():
    
    
    check_access_full_geneset (df_genome_plus, dic_genark_meta, dic_ensembl_meta, df_zoonomia)
    

def test_download_geneset():
    print('')

def test_convert_chr_label_geneset():
    print('')

## ---------------------------------------------
## gencube sequence
## ---------------------------------------------
def test_check_access_full_sequence():
    print('')

def test_download_sequence():
    print('')

## ---------------------------------------------
## gencube annotation
## ---------------------------------------------
def test_check_access_full_annotation():
    print('')

def test_download_annotation():
    print('')

## ---------------------------------------------
## gencube crossgenome
## ---------------------------------------------
def test_check_access_full_crossgenome():
    print('')

def test_download_crossgenome():
    print('')

## -----------------------------------------------------------
## gencube geneset, sequence, annotation, crossgenome
## -----------------------------------------------------------
def test_add_string():
    print('')

## -----------------------------------------------------------
## gencube genome, geneset, sequence, annotation, crossgenome
## -----------------------------------------------------------
# check_url(url, status=False, verify=True, show_output=True)
@pytest.mark.parametrize("status_code, expected_output, expected_status", [
    (200, True, 200),
    (301, True, 301),
    (302, True, 302),
    (307, True, 307),
    (403, False, 403),
    (404, False, 404),
    (429, False, 429),
    (500, False, 500),
    (400, False, 400),
])
def test_check_url(status_code, expected_output, expected_status):
    mock_response = MagicMock()
    mock_response.status_code = status_code
    mock_response.headers = {'Location': 'http://redirected.url'} if status_code in (301, 302, 307) else {}

    with patch('requests.head', return_value=mock_response):
        assert check_url('http://example.com', show_output=False) == expected_output
        assert check_url('http://example.com', status=True, show_output=False) == expected_status

def test_check_url_exception():
    with patch('requests.head', side_effect=requests.exceptions.RequestException("Error")):
        assert not check_url('http://example.com', show_output=False)


# list_ftp_directory(host, directory)        
def test_list_ftp_directory():
    # Mock the FTP connection
    mock_ftp = MagicMock()
    
    # Set up the mock to return a list of files when nlst is called
    mock_ftp.nlst.return_value = ['file1.txt', 'file2.txt', 'dir1', 'dir2']
    
    # Patch the ftplib.FTP object to return our mock instead of making a real connection
    with patch('ftplib.FTP', return_value=mock_ftp):
        # Call the function with test parameters
        result = list_ftp_directory('test_host', 'test_directory')
        
        # Verify the correct directory was accessed
        mock_ftp.cwd.assert_called_with('test_directory')
        
        # Verify the FTP connection was closed
        mock_ftp.quit.assert_called_once()
        
        # Verify the result matches our expectation
        assert result == ['file1.txt', 'file2.txt', 'dir1', 'dir2']

def test_list_ftp_directory_exception():
    # Mock the FTP connection to raise an error on login
    mock_ftp = MagicMock()
    mock_ftp.login.side_effect = ftplib.error_perm("530 Login incorrect.")
    
    with patch('ftplib.FTP', return_value=mock_ftp):
        with pytest.raises(ftplib.error_perm):
            list_ftp_directory('test_host', 'test_directory')


# list_http_folders(url)
def test_list_http_folders():
    # HTML content to be returned by the mock request
    html_content = """
    <html>
    <body>
        <a href="folder1/">folder1/</a>
        <a href="folder2/">folder2/</a>
        <a href="file1.txt">file1.txt</a>
        <a href="/folder3/">/folder3/</a>
    </body>
    </html>
    """
    
    # Mock the response object
    mock_response = MagicMock()
    mock_response.content = html_content
    mock_response.status_code = 200
    
    with patch('requests.get', return_value=mock_response):
        # Call the function with the test URL
        folders = list_http_folders('http://example.com')
        
        # Check the returned folder list
        assert folders == ['folder1', 'folder2']

def test_list_http_folders_request_exception():
    with patch('requests.get', side_effect=requests.exceptions.RequestException("Error")):
        # Call the function with the test URL
        folders = list_http_folders('http://example.com')
        
        # Check that the function returns an empty list in case of an exception
        assert folders == []


# def download_csv(url, verify=True):
def test_download_csv():
    # CSV content to be returned by the mock request
    csv_content = "col1,col2,col3\nval1,val2,val3"

    # Mock the response object
    mock_response = MagicMock()
    mock_response.text = csv_content
    mock_response.status_code = 200

    with patch('requests.get', return_value=mock_response):
        # Call the function with the test URL
        result = download_csv('http://example.com/test.csv')
        
        # Check that the result is a StringIO object
        assert isinstance(result, StringIO)
        
        # Check that the content of the StringIO object matches the expected CSV content
        assert result.getvalue() == csv_content

def test_download_csv_request_exception():
    with patch('requests.get', side_effect=requests.exceptions.RequestException("Error")):
        with pytest.raises(requests.exceptions.RequestException):
            download_csv('http://example.com/test.csv')


# save_metadata (df, function, keywords, level, now)
def test_save_metadata():
    # Create a sample DataFrame
    df = pd.DataFrame({
        'A': [1, 2, 3],
        'B': [4, 5, 6],
        'NCBI': [7, 8, 9]
    })
    
    # Expected DataFrame after dropping 'NCBI' column
    expected_df = df.drop(columns=['NCBI'])
    
    # Mock parameters
    function = 'test_function'
    keywords = ['keyword1', 'keyword2']
    level = 'level1,level2'
    now = '20220101'
    
    # Expected output file name
    expected_out_name = "Meta_test_function_2keywords_level1-level2_20220101.txt"
    
    # Patch the to_csv method of DataFrame
    with patch.object(pd.DataFrame, 'to_csv') as mock_to_csv:
        save_metadata(df, function, keywords, level, now)
        
        # Check that to_csv was called once with correct parameters
        mock_to_csv.assert_called_once_with(expected_out_name, sep='\t', index=False)
        
        # Check the DataFrame passed to to_csv
        args, kwargs = mock_to_csv.call_args
        pd.testing.assert_frame_equal(args[0], expected_df)
        
def test_save_metadata_single_keyword():
    # Create a sample DataFrame
    df = pd.DataFrame({
        'A': [1, 2, 3],
        'B': [4, 5, 6],
        'NCBI': [7, 8, 9]
    })
    
    # Expected DataFrame after dropping 'NCBI' column
    expected_df = df.drop(columns=['NCBI'])
    
    # Mock parameters
    function = 'test_function'
    keywords = ['keyword1']
    level = 'level1,level2'
    now = '20220101'
    
    # Expected output file name
    expected_out_name = "Meta_test_function_keyword1_level1-level2_20220101.txt"
    
    # Patch the to_csv method of DataFrame
    with patch.object(pd.DataFrame, 'to_csv') as mock_to_csv:
        save_metadata(df, function, keywords, level, now)
        
        # Check that to_csv was called once with correct parameters
        mock_to_csv.assert_called_once_with(expected_out_name, sep='\t', index=False)
        
        # Check the DataFrame passed to to_csv
        args, kwargs = mock_to_csv.call_args
        pd.testing.assert_frame_equal(args[0], expected_df)


# calculate_md5(filename)
@pytest.fixture
def setup_test_file(tmpdir):
    """Setup a temporary file for testing."""
    test_file = tmpdir.join("testfile.txt")
    test_content = b"Test content for MD5 checksum"
    with open(test_file, 'wb') as f:
        f.write(test_content)
    return str(test_file), test_content

def test_calculate_md5(setup_test_file):
    test_file, test_content = setup_test_file
    
    # Calculate the expected MD5 hash
    expected_md5 = hashlib.md5(test_content).hexdigest()
    
    # Call the function to test
    calculated_md5 = calculate_md5(test_file)
    
    # Assert that the calculated MD5 matches the expected MD5
    assert calculated_md5 == expected_md5





# download_genome_url(url, local_filename=None, url_md5sum=None, verify=True)
# Sample data
MD5SUM_DATA = """
0d6915df43fc5b9cf988e529fa85800d  ./GCF_007004135.1_ASM700413v1_feature_table.txt.gz
a8361829aa947b850c8ce5c074c1d90b  ./GCF_007004135.1_ASM700413v1_genomic.fna.gz
"""
GENOME_DATA = b"FAKE GENOME DATA"

@pytest.fixture
def setup_environment(tmpdir):
    """Setup a temporary download folder for testing."""
    download_folder = tmpdir.mkdir("downloads")
    return str(download_folder)

def test_download_genome_url(setup_environment, requests_mock):
    url = "https://fakeurl.com/GCF_007004135.1_ASM700413v1_genomic.fna.gz"
    url_md5sum = "https://fakeurl.com/md5sum.txt"
    
    # Mock the MD5SUM URL
    requests_mock.get(url_md5sum, text=MD5SUM_DATA)
    # Mock the file download URL
    requests_mock.get(url, content=GENOME_DATA)
    
    local_filename = os.path.join(setup_environment, "GCF_007004135.1_ASM700413v1_genomic.fna.gz")
    
    with patch('gencube.utils.DOWNLOAD_FOLDER_NAME', setup_environment):
        # Call the function to test
        download_genome_url(url, url_md5sum=url_md5sum)
    
        # Check if the file was downloaded
        assert os.path.exists(local_filename)
        
        # Check the file content
        with open(local_filename, 'rb') as f:
            file_content = f.read()
            assert file_content == GENOME_DATA
        
        # Calculate the MD5 sum of the downloaded file
        md5_download = calculate_md5(local_filename)
        print(f"Expected MD5: a8361829aa947b850c8ce5c074c1d90b")
        print(f"Calculated MD5: {md5_download}")
        
        # Check the MD5 sum of the downloaded file
        assert md5_download == "a8361829aa947b850c8ce5c074c1d90b"

def test_md5sum_mismatch(setup_environment, requests_mock):
    url = "https://fakeurl.com/GCF_007004135.1_ASM700413v1_genomic.fna.gz"
    url_md5sum = "https://fakeurl.com/md5sum.txt"
    
    # Mock the MD5SUM URL with incorrect MD5 sum
    requests_mock.get(url_md5sum, text=MD5SUM_DATA.replace("a8361829aa947b850c8ce5c074c1d90b", "incorrect_md5sum"))
    # Mock the file download URL
    requests_mock.get(url, content=GENOME_DATA)
    
    local_filename = os.path.join(setup_environment, "GCF_007004135.1_ASM700413v1_genomic.fna.gz")
    
    with patch('gencube.utils.DOWNLOAD_FOLDER_NAME', setup_environment):
        # Call the function to test
        download_genome_url(url, url_md5sum=url_md5sum)
        
        # Ensure the file does not exist after failed download
        assert not os.path.exists(local_filename)

def test_download_without_md5sum(setup_environment, requests_mock):
    url = "https://fakeurl.com/GCF_007004135.1_ASM700413v1_genomic.fna.gz"
    
    # Mock the file download URL
    requests_mock.get(url, content=GENOME_DATA)
    
    local_filename = os.path.join(setup_environment, "GCF_007004135.1_ASM700413v1_genomic.fna.gz")
    
    with patch('gencube.utils.DOWNLOAD_FOLDER_NAME', setup_environment):
        # Call the function to test without md5sum URL
        download_genome_url(url)
        
        # Check if the file was downloaded
        assert os.path.exists(local_filename)
        
        # Check the file content
        with open(local_filename, 'rb') as f:
            assert f.read() == GENOME_DATA




# get_md5 (url_md5sum, verify=True
MD5SUM_DATA = """
0d6915df43fc5b9cf988e529fa85800d  ./GCF_007004135.1_ASM700413v1_feature_table.txt.gz
a8361829aa947b850c8ce5c074c1d90b  ./GCF_007004135.1_ASM700413v1_genomic.fna.gz
"""

def test_get_md5(requests_mock):
    url_md5sum = "https://fakeurl.com/md5sum.txt"
    
    # Mock the MD5SUM URL
    requests_mock.get(url_md5sum, text=MD5SUM_DATA)
    
    # Call the function to test
    df = get_md5(url_md5sum)
    
    # Check the DataFrame content
    expected_data = [
        ["0d6915df43fc5b9cf988e529fa85800d", "./GCF_007004135.1_ASM700413v1_feature_table.txt.gz"],
        ["a8361829aa947b850c8ce5c074c1d90b", "./GCF_007004135.1_ASM700413v1_genomic.fna.gz"]
    ]
    expected_df = pd.DataFrame(expected_data, columns=['MD5', 'Filename'])
    
    pd.testing.assert_frame_equal(df, expected_df)

# delete_file(file_path)
@pytest.fixture
def create_test_file(tmpdir):
    """Create a temporary file for testing."""
    test_file = tmpdir.join("testfile.txt")
    with open(test_file, 'w') as f:
        f.write("Test content")
    return str(test_file)

def test_delete_file(create_test_file):
    test_file = create_test_file
    
    # Ensure the file exists before deletion
    assert os.path.exists(test_file)
    
    # Call the function to test
    delete_file(test_file)
    
    # Check that the file has been deleted
    assert not os.path.exists(test_file)

def test_delete_nonexistent_file(tmpdir):
    test_file = tmpdir.join("nonexistent.txt")
    
    # Ensure the file does not exist
    assert not os.path.exists(test_file)
    
    # Call the function to test
    delete_file(str(test_file))
    
    # Check that the function completes without error and the file still does not exist
    assert not os.path.exists(test_file)
    
    
# make_chr_dataframe (organism, assembly_id)
# Sample data and expected output
LS_ASSEMBLY_REPORT_LABEL = [
    'Sequence-Name', 'Sequence-Role', 'Assigned-Molecule', 'Assigned-Molecule-Location/Type', 
    'GenBank-Accn', 'Relationship', 'RefSeq-Accn', 'Assembly-Unit', 'Sequence-Length', 'UCSC-style-name'
]

SAMPLE_REPORT_CONTENT = """
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
## GenBank Unit Accession    RefSeq Unit Accession    Assembly-Unit name
## GCA_000000145.4   GCF_000000145.4   Primary Assembly
## GCA_013123455.1   GCF_000184225.1   non-nuclear
#
# Ordered by chromosome/plasmid; the chromosomes/plasmids are followed by
# unlocalized scaffolds.
# Unplaced scaffolds are listed at the end.
# RefSeq is equal or derived from GenBank object.
#
# Sequence-Name    Sequence-Role    Assigned-Molecule    Assigned-Molecule-Location/Type    GenBank-Accn    Relationship    RefSeq-Accn    Assembly-Unit    Sequence-Length    UCSC-style-name
1    assembled-molecule    1    Chromosome    CM000001.4    =    NC_006583.4    Primary Assembly    122014068    chr1
2    assembled-molecule    2    Chromosome    CM000002.4    =    NC_006584.4    Primary Assembly    82037489    chr2
MT    assembled-molecule    MT    Mitochondrion    CM023446.1    <>    na    non-nuclear    16735    na
MT    assembled-molecule    MT    Mitochondrion    na    <>    NC_002008.4    non-nuclear    16727    chrM
"""

mock_report_content = """\
1\tassembled-molecule\t1\tChromosome\tCM000001.4\t=\tNC_006583.4\tPrimary Assembly\t122014068\tchr1
MT\tassembled-molecule\tMT\tMitochondrion\tCM023446.1\t<>\tna\tnon-nuclear\t16735\tna
MT\tassembled-molecule\tMT\tMitochondrion\tna\t<>\tNC_002008.4\tnon-nuclear\t16727\tchrM
"""

@pytest.fixture
def setup_environment(tmp_path):
    download_folder = tmp_path / "downloads"
    output_folder = tmp_path / "outputs"
    download_folder.mkdir()
    output_folder.mkdir()
    
    in_report_name = download_folder / "test_organism-test_assembly_assembly_report.txt"
    with open(in_report_name, "w") as f:
        f.write(mock_report_content)
    
    return download_folder, output_folder, in_report_name

def test_make_chr_dataframe(setup_environment):
    download_folder, output_folder, in_report_name = setup_environment
    organism = "test_organism"
    assembly_id = "test_assembly"
    out_report_name = os.path.join(output_folder, f"{organism}-{assembly_id}_chr_name.txt")
    
    with patch('gencube.utils.DOWNLOAD_FOLDER_NAME', str(download_folder)), patch('gencube.utils.OUT_FOLDER_NAME', str(output_folder)):
        df_result = make_chr_dataframe(organism, assembly_id)
        
        # 예상되는 결과
        expected_data = {
            "genbank": ["CM000001.4", "CM023446.1"],
            "refseq": ["NC_006583.4", "NC_002008.4"],
            "ucsc": ["chr1", "chrM"],
            "ensembl": ["1", "MT"],
            "gencode": ["chr1", "chrM"]
        }
        expected_df = pd.DataFrame(expected_data)
        
        pd.testing.assert_frame_equal(df_result.reset_index(drop=True), expected_df.reset_index(drop=True))
        
        # 생성된 파일 확인
        assert os.path.exists(out_report_name)
        df_output = pd.read_csv(out_report_name, sep='\t')
        pd.testing.assert_frame_equal(df_output.reset_index(drop=True), expected_df.reset_index(drop=True))


## ---------------------------------------------
## gencube seqmeta
## ---------------------------------------------
# make_query(organism, strategy, source, layout, keywords, operators, now):

# Test cases
now = '240601-041301'
@pytest.mark.parametrize("organism, strategy, source, layout, keywords, operators, now, expected_query, expected_out_name", [
    # Test case 1: Single organism
    ("human", "", "", "", "", "", now, '("homo sapiens"[Organism])', "human_240601-041301"),
    # Test case 2: Multiple organisms
    ("human,mouse", "", "", "", "", "", now, '("homo sapiens"[Organism] OR "mus musculus"[Organism])', "2organisms_240601-041301"),
    # Test case 3: Strategy included
    ("human", "rna", "", "", "", "", now, '(("homo sapiens"[Organism]) AND ("rna seq"[Strategy]))', "human_rna_240601-041301"),
    # Test case 4: Source included
    ("human", "", "genomic", "", "", "", now, '(("homo sapiens"[Organism]) AND ("genomic"[Source]))', "human_genomic_240601-041301"),
    # Test case 5: Layout included
    ("human", "", "", "paired", "", "", now, '(("homo sapiens"[Organism]) AND "paired"[Layout])', "human_paired_240601-041301"),
    # Test case 6: Keywords and operators included
    ("human", "", "", "", "cancer", "and", now, '(("homo sapiens"[Organism]) AND "cancer")', "human_and-cancer_240601-041301"),
    
    # Test case 7: Multiple keywords and operators
    ("human", "", "", "", "cancer,tumor", "and,or", now, '((("homo sapiens"[Organism]) AND "cancer") OR "tumor")', "human_and-cancer_or-tumor_240601-041301"),
    # Test case 8: Complex case with all fields
    ("human,mouse", "rna", "genomic", "paired", "cancer,tumor", "and,or", now, 
     '(((((("homo sapiens"[Organism] OR "mus musculus"[Organism]) AND ("rna seq"[Strategy])) AND ("genomic"[Source])) AND "paired"[Layout]) AND "cancer") OR "tumor")', 
     "2organisms_rna_genomic_paired_and-cancer_or-tumor_240601-041301"),
])
def test_make_query(organism, strategy, source, layout, keywords, operators, now, expected_query, expected_out_name):
    query, out_name = make_query(organism, strategy, source, layout, keywords, operators, now)
    assert query == expected_query
    assert out_name == expected_out_name
    

# search_sra(query)
# The function to be tested
def search_sra(query):
    from Bio import Entrez
    print('# Search experimental sequencing data in NCBI SRA database')
    print(f'  Search query: \n  {query}\n')
    
    handle = Entrez.esearch(db="sra",  # Database to search
                            term=query,  # Search term
                            retmax=10000  # Number of results to return
                            )
    record = Entrez.read(handle)
    
    n = int(record["Count"])
    if n < 2:
        print(f'  Total {record["Count"]} sample-level ID is searched.\n')
    else:
        print(f'  Total {record["Count"]} sample-level IDs are searched.\n')
    
    return record

# Mock data
mock_query = "(homo sapiens[Organism])"
mock_record_single = {"Count": "1", "IdList": ["SRR123456"]}
mock_record_multiple = {"Count": "10", "IdList": ["SRR123456", "SRR123457", "SRR123458", "SRR123459", "SRR123460",
                                                 "SRR123461", "SRR123462", "SRR123463", "SRR123464", "SRR123465"]}

@patch('Bio.Entrez.esearch')
@patch('Bio.Entrez.read')
def test_search_sra_single_result(mock_entrez_read, mock_entrez_esearch, capsys):
    # Setup the mock to return the mock_record_single
    mock_entrez_esearch.return_value = MagicMock()
    mock_entrez_read.return_value = mock_record_single

    # Call the function
    result = search_sra(mock_query)

    # Capture the output
    captured = capsys.readouterr()

    # Assert the expected results
    assert result == mock_record_single
    assert "Total 1 sample-level ID is searched." in captured.out

@patch('Bio.Entrez.esearch')
@patch('Bio.Entrez.read')
def test_search_sra_multiple_results(mock_entrez_read, mock_entrez_esearch, capsys):
    # Setup the mock to return the mock_record_multiple
    mock_entrez_esearch.return_value = MagicMock()
    mock_entrez_read.return_value = mock_record_multiple

    # Call the function
    result = search_sra(mock_query)

    # Capture the output
    captured = capsys.readouterr()

    # Assert the expected results
    assert result == mock_record_multiple
    assert "Total 10 sample-level IDs are searched." in captured.out


# fetch_meta(ls_id)


def test_convert_format():
    print('')

def test_save_seq_metadata():
    print('')




if __name__ == "__main__":
    pytest.main()