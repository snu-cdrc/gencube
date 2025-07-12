# For accessing NCBI Entrez databases
from Bio import Entrez
# For convertion of data format or data manipulation
import pandas as pd           # For data manipulation and analysis
import numpy as np            # For numerical operations
import xmltodict              # For converting XML data to Python dictionaries
import json                   # For parsing JSON data
import ast
# For access to or manipulation of web data
import requests               # For making HTTP requests
import urllib3                # For advanced HTTP client functionalities
import ftplib                 # For handling FTP connectionssave_metadata
from bs4 import BeautifulSoup # For parsing HTML and XML documents
from io import StringIO       # For in-memory file-like objects

from pathlib import Path      # For handling filesystem paths in a platform-independent way
import os                     # For interacting with the operating system
import subprocess             # For running and interacting with system commands and processes

import re                     # For regular expression operations
import gzip                   # For reading and writing gzip-compressed files

import hashlib                # For generating MD5 checksums
# For displayed output
from tqdm import tqdm         # For displaying progress bars
from tabulate import tabulate # For printing data in a tabular format
# For management of time information
from datetime import datetime # For manipulating date and time objects
import time                   # For handling time-related tasks
# For operating system and architecture check
import platform
# For multi-threading
from concurrent.futures import ThreadPoolExecutor, as_completed

# To suppress the SSL warnings when using verify=False in the requests library
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# To suppress the pandas warnings
import warnings
from pandas.errors import PerformanceWarning
warnings.simplefilter("ignore", PerformanceWarning)

# Constant variables
from .constants import (
    NCBI_FTP_HOST,
    NCBI_FTP_URL,
    #ENSEMBL_FTP_HOST,
    #ENSEMBL_RAPID_FTP_URL,
    #ENSEMBL_RM_FTP_URL,
    ENSEMBL_BETA_FTP_URL,
    ENSEMBL_TREE_JSON,
    GENARK_URL,
    ZOONOMIA_URL,
    ZOONOMIA_META,
    UCSC_KENT_URL,
    LS_NCBI_ASSEMBLY_META_KEY,
    LS_NCBI_ASSEMBLY_META_LABEL,
    LS_SRA_META_STUDY_KEY,
    LS_SRA_META_EXPERIMENT_KEY,
    LS_SRA_META_STUDY_LABEL,
    LS_SRA_META_EXPERIMENT_LABEL,
    LS_ASSEMBLY_REPORT_LABEL,
    DIC_ZOONOMIA,
    )

# Global variables
DOWNLOAD_FOLDER_NAME = 'gencube_raw_download'
OUT_FOLDER_NAME = 'gencube_output'
script_dir = os.path.dirname(os.path.abspath(__file__))
## -----------------------------------------------------------
## gencube all subcommands
## -----------------------------------------------------------
# Get email and api key information for E-utilities
def get_entrez_info():
    global api_key

    home_dir = Path.home()
    config_file = home_dir / '.gencube_entrez_info'

    if config_file.exists():
        # Read the file
        with open(config_file, 'r') as file:
            lines = file.readlines()
            email = lines[1].split('=')[1].strip()
            api_key = lines[2].split('=')[1].strip()
    else:
        # Prompt the user for email with validation
        print(
            "\n"
            '--------------------------------------------------------------------------------------\n'
            "All gencube subcommands use NCBI's Entrez Utilities (E-Utilities), requiring an email.\n"
            "Without an NCBI API key, you can make 3 requests per second; with an NCBI API key,\n"
            "this limit increases to 10 requests per second.\n"
            "If you submit your NCBI API key, you can perform tasks at more than three times \n"
            "the speed when using the seqmeta subcommand, especially when fetching metadata.\n"
            "If possible, it is recommended to submit your API key.\n"
            "\n"
            "The submitted information is stored in the file at the path below for reuse.\n"
            f"Path: {config_file}\n"
            "\n"
            "If you want to resubmit this information, run '$ gencube info'.\n"
            '--------------------------------------------------------------------------------------\n'
        )
        # Prompt the user for email
        email = ''
        while not re.match(r"[^@]+@[^@]+\.[^@]+", email):
            email = input("Email address: ")
            print('')
            if not re.match(r"[^@]+@[^@]+\.[^@]+", email):
                print('Invalid email format. Please try again\n')

        # Prompt the user for API key with skip option
        api_key = ''
        api_skip = False
        while not re.match(r'^[a-zA-Z0-9]{36}$', api_key) and not api_skip:
            api_key = input("NCBI API key (type 'no' to skip): ").strip()
            if api_key.lower() in ['no', 'n']:
                api_key = ''
                api_skip = True
            elif not re.match(r'^[a-zA-Z0-9]{36}$', api_key):
                print(
                    'The key should be a 36-character mix of letters and numbers.\n'
                    'Invalid API key format. Please try again\n'
                    )

        # Save the information to the file
        with open(config_file, 'w') as out_f:
            out_f.write("# This information is used in E-utilities\n")
            out_f.write(f"email = {email}\n")
            out_f.write(f"api_key = {api_key}\n")

        # Inform the user where the information is saved
        print(f'\n!! The information is saved to "{config_file}"\n')

        return True

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

## -----------------------------------------------------------
## gencube genome, geneset, sequence, annotation, crossgenome
## -----------------------------------------------------------
# Check starting time
def check_now ():
    now = datetime.now()
    formatted_now = now.strftime("%y%m%d_%H%M%S")
    return formatted_now

# Merge strings with comma
def add_string (pre, string):
    if pre:
        pre += f', {string}'
    else:
        pre = string

    return pre

# Check right and wrong inputs
def check_argument (input, ls, str):

    if input:
        ls_input_raw = input.split(',')
        ls_input = []
        ls_wrong = []

        # Check right inputs
        for input_tmp in ls:

            if input_tmp in ls_input_raw:
                ls_input.append(input_tmp)

        # Check wrong inputs
        for type_raw in ls_input_raw:
            if type_raw not in ls_input:
                ls_wrong.append(type_raw)

        if len(ls_wrong) > 0:
            print_arg = ''
            for wrong in ls_wrong:
                print_arg = add_string(print_arg, wrong)

            print(f'Invalid argument {str}: {print_arg}')

        return ls_wrong
    else:
        tmp = []
        return tmp

# Check url
def check_url(url, verify=True, show_output=True, file_name=''):
    try:
        if not file_name:
            file_name = url.split('/')[-1]
        # Use the HEAD request to fetch only the headers of the resource and allow redirects.
        response = requests.head(url, allow_redirects=True, verify=verify)
        # Check for successful response and file accessibility.
        if response.status_code == 200:
            return True
        # Check if the response is a redirect, and the file might still be accessible.
        elif response.status_code in (301, 302, 307):
            if show_output:
                print(f"  {file_name}: !! Redirected to {response.headers['Location']}")
            return True
        elif response.status_code == 403:
            if show_output:
                print(f"  {file_name}: !! Access denied. File cannot be downloaded")
            return False
        elif response.status_code == 404:
            if show_output:
                print(f"  {file_name}: !! File not found")
            return False
        elif response.status_code == 429:
            if show_output:
                print(f"  {file_name}: !! Too many requests. Please try again later")
            return False
        elif response.status_code >= 500:
            if show_output:
                print(f"  {file_name}: !! Server error. Please try again later")
            return False
        else:
            # Handle any other unexpected status codes.
            if show_output:
                print(f"  {file_name}: !! Unexpected status code: {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        # Handle any exceptions that occur during the request.
        if show_output:
            print(f"  {file_name}: !! Error checking URL: {e}")
        return False

# Fetch information of folders located in ftp server
def list_ftp_directory(host, directory):
    ftp = ftplib.FTP(host)
    ftp.login()  # Login (anonymous access)
    ftp.cwd(directory)  # Change to the specified directory
    files = ftp.nlst()  # Retrieve directory listing
    ftp.quit()  # Close FTP connection
    return files

# Fetch information of folders located in ftp server
def list_http_folders(url):
    try:
        # Send a GET request to the URL
        response = requests.get(url, verify=False)
        response.raise_for_status()  # Check for request errors

        # Parse the HTML content
        soup = BeautifulSoup(response.content, 'html.parser')

        # Find all links that are directories
        folders = []
        for link in soup.find_all('a'):
            href = link.get('href')
            if href and href.endswith('/') and not href.startswith('/'):
                folders.append(href.replace('/', ''))

        return folders

    except requests.exceptions.RequestException as e:
        print(f"Error accessing URL: {e}")
        return []

# Fetch information of files located in an HTTP server directory
def list_http_files(url):
    try:
        # Send a GET request to the URL
        response = requests.get(url, verify=False)
        response.raise_for_status()  # Check for request errors

        # Parse the HTML content
        soup = BeautifulSoup(response.content, 'html.parser')

        # Find all links that are files
        files = []
        for link in soup.find_all('a'):
            href = link.get('href')
            if href and not href.endswith('/') and not href.startswith('/'):
                files.append(href)

        return files

    except requests.exceptions.RequestException as e:
        print(f"Error accessing URL: {e}")
        return []

# Fetch information using url
def download_csv(url, verify=True):
    # Download the file with SSL verification disabled
    response = requests.get(url, verify=verify)
    response.raise_for_status()  # Check for request errors
    return StringIO(response.text)

# Print warning message
def print_warning (df, n):
    if df.shape[0] > n:
        print(f'Warning: The number of genomes searched is very large (> {n}), resulting in extensive output on the terminal.')
        print('If reviewing all results in the terminal is difficult, consider using the "--metadata" option to save metadata.\n')

# Make download and output folder
def mkdir_raw_output (folder_name):
    # Make download folder
    if not os.path.exists(DOWNLOAD_FOLDER_NAME):
        os.mkdir(DOWNLOAD_FOLDER_NAME)
    if not os.path.exists(OUT_FOLDER_NAME):
        os.mkdir(OUT_FOLDER_NAME)
    global out_subfolder
    out_subfolder = f'{OUT_FOLDER_NAME}/{folder_name}'
    if not os.path.exists(out_subfolder):
        os.mkdir(out_subfolder)

# Save metadata
def save_metadata (df_meta, function, keywords, level, now):
    # Remove some columns not required in metadata
    if function == 'genome':
        df_meta = df_meta.drop(columns=['NCBI'])

    # Save the metadata
    if len(keywords) != 1:
        keyword = f'{len(keywords)}keywords'
    else:
        keyword = keywords[0]

    out_name = f"Meta_{function}_{keyword.replace(' ', '-')}_{level.replace(',', '-')}_{now}.txt"
    df_meta.to_csv(f'{out_subfolder}/{out_name}', sep='\t', index=False)

    print('# Metadata are saved:')
    print(f'  {out_name}\n')

# Calculate and return the MD5 hash of a file
def calculate_md5(filename):
    return hashlib.md5(open(filename,'rb').read()).hexdigest()

# Download genome file and check md5sum
def download_genome_url(url, local_filename=None, url_md5sum=None, verify=True, recursive=None):
    # File name
    raw_filename = url.split('/')[-1]
    if local_filename is None:
        local_filename = raw_filename
    path_local_file = os.path.join(DOWNLOAD_FOLDER_NAME, local_filename)

    if url_md5sum:
        response = requests.get(url_md5sum, verify=verify)
        # This will raise an error for bad responses
        response.raise_for_status()
        # Split the content by new lines and parse it
        data = response.text.split('\n')
        data = [line.split() for line in data if line]  # Split each line into parts and remove empty lines
        # Create a DataFrame
        df = pd.DataFrame(data, columns=['MD5', 'Filename'])

        # Get md5sum of file
        for idx in df.index:
            md5sum_filename = df.loc[idx, 'Filename'].split('/')[-1]
            if raw_filename == md5sum_filename:
                md5sum_original = df.loc[idx, 'MD5']

    count = 0
    recursive_download = False
    while True:
        if url_md5sum:
            if os.path.exists(path_local_file):
                # Get md5sum from download file
                md5sum_download = calculate_md5(path_local_file)
                # If md5sums are same between original and download file
                if md5sum_download == md5sum_original:
                    if count == 0:
                        if '_assembly_report.txt' in path_local_file:
                            print('  Assembly Report was already downloaded')
                        elif 'fa.gz' in path_local_file:
                            print(f'  {local_filename} was already downloaded')
                    break
                # If md5sums are different between original and download file
                else:
                    print(
                        '  md5sum value is not same with original file\n'
                        f'  {md5sum_original}: original file \n'
                        f'  {md5sum_download}: download file \n'
                    )
                    if count < 1:
                        print('  Re-try downloading the file')
                        os.remove(path_local_file)
                    else:
                        print('  Try download file again')
                        os.remove(path_local_file)
                        break
        else:
            # If there is no url_md5sum (can't be checked)
            if 'fa.gz' in local_filename and os.path.exists(path_local_file):
                if recursive is False:
                    if count == 0:
                        print(f'  {local_filename} was already downloaded')
                    break
                else:
                    if recursive_download:
                        break
                    else:
                        os.remove(path_local_file)

        # Check the size of the local file if it already exists.
        existing_size = os.path.getsize(path_local_file) if os.path.exists(path_local_file) else 0

        # Add the Range header to the HTTP request to download the part of the file that is not yet downloaded.
        headers = {"Range": f"bytes={existing_size}-"}
        response = requests.get(url, headers=headers, stream=True, verify=verify)
        total_size = int(response.headers.get('content-length', 0)) + existing_size

        # Open the file in append mode and proceed with the download.
        with open(path_local_file, 'ab') as f, tqdm(
            desc='  ' + local_filename,
            initial=existing_size,
            total=total_size,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                bar.update(len(chunk))

        count += 1
        recursive_download = True

    if not url_md5sum:
        print("  !! MD5 check failed. If issues, delete and retry, or use --recursive")

# Download file and check md5sum
def download_url(url, local_filename=None, url_md5sum=None, verify=True, recursive=False, warning=True, out_path=False):
    # File name
    raw_filename = url.split('/')[-1]
    if local_filename is None:
        local_filename = raw_filename
    if out_path:
        path_local_file = os.path.join(out_path, local_filename)
    else:
        path_local_file = os.path.join(DOWNLOAD_FOLDER_NAME, local_filename)

    md5sum = False
    if url_md5sum:
        response = requests.get(url_md5sum, verify=verify)
        # This will raise an error for bad responses
        response.raise_for_status()
        # Split the content by new lines and parse it
        data = response.text.split('\n')
        data = [line.split() for line in data if line]  # Split each line into parts and remove empty lines
        # Create a DataFrame
        df = pd.DataFrame(data, columns=['MD5', 'Filename'])

        # Get md5sum of file
        for idx in df.index:
            md5sum_filename = df.loc[idx, 'Filename'].split('/')[-1]
            if raw_filename == md5sum_filename:
                md5sum_original = df.loc[idx, 'MD5']
                md5sum = True

    count = 0
    recursive_download = False
    while True:
        if url_md5sum and md5sum:
            if os.path.exists(path_local_file):
                # Get md5sum from download file
                md5sum_download = calculate_md5(path_local_file)
                # If md5sums are same between original and download file
                if md5sum_download == md5sum_original:
                    if count == 0:
                        print(f'  {local_filename} was already downloaded')
                    break
                # If md5sums are different between original and download file
                else:
                    print(
                        '  md5sum value is not same with original file\n'
                        f'  {md5sum_original}: original file \n'
                        f'  {md5sum_download}: download file \n'
                    )
                    if count < 1:
                        print('  Re-try downloading the file')
                        os.remove(path_local_file)
                    else:
                        print('  Try download file again')
                        os.remove(path_local_file)
                        break
        else:
            # If there is no url_md5sum (can't be checked)
            if os.path.exists(path_local_file):
                if not recursive:
                    if count == 0:
                        print(f'  {local_filename} was already downloaded')
                    break
                else:
                    if recursive_download:
                        break
                    else:
                        os.remove(path_local_file)

        # Check the size of the local file if it already exists.
        existing_size = os.path.getsize(path_local_file) if os.path.exists(path_local_file) else 0

        # Add the Range header to the HTTP request to download the part of the file that is not yet downloaded.
        headers = {"Range": f"bytes={existing_size}-"}
        response = requests.get(url, headers=headers, stream=True, verify=verify)
        total_size = int(response.headers.get('content-length', 0)) + existing_size

        # Open the file in append mode and proceed with the download.
        with open(path_local_file, 'ab') as f, tqdm(
            desc='  ' + local_filename,
            initial=existing_size,
            total=total_size,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                bar.update(len(chunk))

        count += 1
        recursive_download = True

    if not url_md5sum or not md5sum:
        if warning:
            print("  !! MD5 check failed. If issues, delete and retry, or use --recursive")

# Fetch md5sum information
def get_md5 (url_md5sum, verify=True):
    response = requests.get(url_md5sum, verify=verify)
    # This will raise an error for bad responses
    response.raise_for_status()
    # Split the content by new lines and parse it
    data = response.text.split('\n')
    data = [line.split() for line in data if line]  # Split each line into parts and remove empty lines
    # Create a DataFrame
    df = pd.DataFrame(data, columns=['MD5', 'Filename'])

    return df

# Remove file
def delete_file(file_path):
    try:
        if os.path.exists(file_path):
            os.remove(file_path)
    except Exception as e:
        print(f"Error deleting {file_path}: {e}")

# Make integrative dataframe of chromosome label of databases
def make_chr_dataframe (organism, assembly_id):
    # read report file
    ls_report = []
    in_report_name = f'{DOWNLOAD_FOLDER_NAME}/{organism}-{assembly_id}_assembly_report.txt'
    out_report_name = f'{OUT_FOLDER_NAME}/{organism}-{assembly_id}_chr_name.txt'

    if os.path.exists(in_report_name):
        with open (in_report_name, 'r') as file:
            for line in file:
                tmp = line.split()
                if len(tmp) > 0 and line[0] != '#':
                    if line[0] == '<':
                        break
                    ls_report.append(line.strip().split('\t'))

    else:
        print(f'  {organism}-{assembly_id}_assembly_report.txt file is not found \n')
        return ''

    df_report_raw = pd.DataFrame(ls_report, columns=LS_ASSEMBLY_REPORT_LABEL)
    df_report = df_report_raw[['GenBank-Accn', 'RefSeq-Accn', 'Assigned-Molecule',  'UCSC-style-name']]
    df_report_edit = df_report.copy()

    # Create ensembl-style and gencode-style label
    # - Ensembl
    df_report_edit.loc[:, 'Ensembl'] = df_report.apply(
        lambda row: row['Assigned-Molecule'] if row['GenBank-Accn'].startswith('CM') or row['RefSeq-Accn'].startswith('NC') else row['GenBank-Accn'],
        axis=1)
    # - Gencode
    df_report_edit.loc[:, 'Gencode'] = df_report.apply(
        lambda row: 'chr' + row['Assigned-Molecule'] if row['GenBank-Accn'].startswith('CM') or row['RefSeq-Accn'].startswith('NC') else row['GenBank-Accn'],
        axis=1)
    df_report_edit['Gencode'] = df_report_edit['Gencode'].replace('chrMT', 'chrM')

    df_report_edit.drop(['Assigned-Molecule'], axis=1, inplace=True)
    # Change column names
    df_report_edit.columns = ['GenBank', 'RefSeq', 'UCSC', 'Ensembl', 'Gencode']

    # Merge rows including chrM information
    if len(df_report_edit[df_report_edit['Ensembl'] == 'MT'].index) == 2:
        # Find the first valid 'GenBank' value that is not 'na' among rows where 'Ensembl' is 'MT'
        genbank_mt = df_report_edit.loc[(df_report_edit['Ensembl'] == 'MT') & (df_report_edit['GenBank'] != 'na'), 'GenBank'].iloc[0]
        # Replace all 'na' values in 'GenBank' for rows where 'Ensembl' is 'MT'
        df_report_edit.loc[(df_report_edit['Ensembl'] == 'MT') & (df_report_edit['GenBank'] == 'na'), 'GenBank'] = genbank_mt
        # After replacement, remove the rows that originally had a 'GenBank' value (that was not 'na')
        df_report_edit = df_report_edit.drop(df_report_edit[(df_report_edit['Ensembl'] == 'MT') & (df_report_edit['RefSeq'] == 'na')].index)
    # UCSC chrM name
    df_report_edit.loc[df_report_edit['Ensembl'] == 'MT', 'UCSC'] = 'chrM'

    # Save dataframe
    #df_report_edit_out = df_report_edit.copy()
    #df_report_edit_out.columns = ['GenBank', 'RefSeq', 'UCSC', 'Ensembl', 'Gencode']
    df_report_edit.to_csv(out_report_name, sep='\t', index=False)

    # Remove report file
    #delete_file(in_report_name)

    return df_report_edit

# Search genome using users keywords
def search_assembly (keywords):
    print('# Search assemblies in NCBI database')
    print(f'  Keyword: {keywords}\n')

    ls_search = []
    for keyword in keywords:

        # Search keyword in NCBI Assembly database
        try:
            handle = Entrez.esearch(db="assembly", term=keyword, retmax="10000")
            record = Entrez.read(handle)
            handle.close()

            assembly_id = record["IdList"]
            # Check if the id list is not empty
            if assembly_id:
                # Fetch detailed information using searched ids in NCBI Assembly database
                handle = Entrez.efetch(db="assembly", id=assembly_id, rettype="docsum")
                record = Entrez.read(handle)
                handle.close()

                ls_search += record['DocumentSummarySet']['DocumentSummary']
            else:
                print(f"  No results found for keyword '{keyword}'\n")
        except Exception as e:
            print(f"  !! An error occurred while searching for keyword '{keyword}': {e} \n")
            continue

    if ls_search:
        return ls_search
    else:
        print("  !! No results found for all keyword")

# Convert the data format
def json_to_dataframe (search, level, refseq, ucsc, latest):
    # Inner function
    def get_values_from_keys(dic):
        ls_value = []
        for key in LS_NCBI_ASSEMBLY_META_KEY:
            if key == 'Synonym_GCA':
                ls_value.append(dic['Synonym']['Genbank'])
            elif key == 'Synonym_GCF':
                ls_value.append(dic['Synonym']['RefSeq'])
            elif key == 'GB_BioProjects':
                if len(dic['GB_BioProjects']) > 0:
                    ls_value.append(dic['GB_BioProjects'][0]['BioprojectAccn'])
                else:
                    ls_value.append('')
            elif key in dic:
                    ls_value.append(dic[key])

        return pd.DataFrame([ls_value], columns=LS_NCBI_ASSEMBLY_META_LABEL)

    n = len(search)
    if n < 2:
        print(f'  Total {len(search)} genomes is searched\n')
    else:
        print(f'  Total {len(search)} genomes are searched\n')
    print('# Filter genomes based on the following criteria')
    for i in range(len(search)):
        df_tmp = get_values_from_keys(search[i])

        if i == 0:
            df_assembly = df_tmp
        else:
            df_assembly = pd.concat([df_assembly, df_tmp])

    df_assembly = df_assembly.reset_index(drop=True)
    # Remove
    df_assembly['Release'] = df_assembly['Release'].apply(lambda x: x.split(' ')[0])
    df_assembly['Update'] = df_assembly['Update'].apply(lambda x: x.split(' ')[0])
    df_assembly['Release (RefSeq)'] = df_assembly['Release (RefSeq)'].apply(lambda x: x.split(' ')[0])
    df_assembly['SeqReleaseDate'] = df_assembly['SeqReleaseDate'].apply(lambda x: x.split(' ')[0])
    df_assembly['SubmissionDate'] = df_assembly['SubmissionDate'].apply(lambda x: x.split(' ')[0])
    df_assembly['LastUpdateDate'] = df_assembly['LastUpdateDate'].apply(lambda x: x.split(' ')[0])

    df_assembly['Level'] = df_assembly['Level'].replace('Complete Genome', 'Complete')
    # Merge RefSeq and Genbank accession in a column. (Priority: RefSeq > Genbank)
    df_assembly['NCBI'] = df_assembly.apply(
        lambda x: x['RefSeq'] if x['RefSeq'] != '' else x['GenBank'],
        axis=1)

    # Check Level (Assembly status)
    ls_level = level.split(',')
    ls_level = [x.capitalize() for x in ls_level] # capitalize the first letter
    df_assembly_out = df_assembly[df_assembly['Level'].isin(ls_level)]

    # Extract latest genomes
    if latest:
        ls_index = []

        for idx in df_assembly_out.index:
            latest_acc = df_assembly_out.loc[idx, 'Latest accession']
            if latest_acc == '':
                ls_index.append(idx)
            elif  latest_acc[:3] == 'GCF':
                ls_index.append(df_assembly_out[df_assembly_out['RefSeq'] == latest_acc].index[0])
            elif  latest_acc[:3] == 'GCA':
                ls_index.append(df_assembly_out[df_assembly_out['GenBank'] == latest_acc].index[0])
            else: # 지울 부분!
                print('Latest accession is other type!')

        ls_index = list(set(ls_index))
        df_assembly_out = df_assembly_out.loc[ls_index]

    # Extract only genoms that RefSeq accession is issued.
    if refseq:
        df_assembly_out = df_assembly_out[df_assembly_out['RefSeq'] != '']
    # Extract only genoms that UCSC name is issued.
    if ucsc:
        df_assembly_out = df_assembly_out[df_assembly_out['UCSC'] != '']
    # Sort by Update
    df_assembly_out = df_assembly_out.sort_values(by='Update', ascending=False).reset_index(drop=True)
    # Print searched result
    print(f'  Level:   {ls_level}')
    print(f'  RefSeq:  {refseq}')
    print(f'  UCSC:    {ucsc}')
    print(f'  Latest:  {latest}\n')

    return df_assembly_out

# Check the available genomes in GenArk and Ensembl (Rapid Release) server
def check_access_database (df, mode):
    if mode == 'genome':
        print('# Check accessibility to GenArk, Ensembl Beta')
    elif mode == 'geneset':
        print('# Check accessibility to GenArk, Ensembl Beta and Zoonomia server')
    elif mode == 'sequence':
        print('# Check accessibility to Ensembl Beta')
    elif mode == 'annotation':
        print('# Check accessibility to GenArk')
    elif mode == 'crossgenome':
        print('# Check accessibility to Ensembl Beta and Zoonomia server')

    ls_genark_mode = ['genome', 'geneset', 'annotation']
    ls_ensembl_mode = ['genome', 'geneset', 'sequence', 'crossgenome']
    ls_zoonomia_mode = ['geneset', 'crossgenome']

    # Check GenArk
    if mode in ls_genark_mode:
        genark_meta_url = GENARK_URL + 'UCSC_GI.assemblyHubList.txt'
        genark_meta = requests.get(genark_meta_url, verify=False).text.split('\n')

        dic_genark_meta = {}
        dic_genark_meta_species = {}
        for line in genark_meta:
            tmp = line.split()
            if len(tmp):
                if tmp[0] != '#':
                    # accession : "assembly"\t"scientific name"
                    dic_genark_meta[line.split('\t')[0]] = line.split('\t')[1]
                    dic_genark_meta_species[line.split('\t')[0]] = line.split('\t')[2]

        df[['GenArk']] = ''
        print(f'  UCSC GenArk  : {len(dic_genark_meta)} genomes across {len(list(set(dic_genark_meta_species.values())))} species')

    """
    # Check Ensembl Rapid Release
    if mode in ls_ensembl_mode:
        ensembl_meta_url = ENSEMBL_RAPID_FTP_URL + 'species_metadata.json'
        ensembl_meta = requests.get(ensembl_meta_url).json()

        dic_ensembl_meta = {}
        for acc in ensembl_meta:
            dic_ensembl_meta[acc['assembly_accession']] = acc['species']

        df[['Ensembl']] = ''
        print(f'  Ensembl Rapid: {len(dic_ensembl_meta)} genomes across {len(list(set(dic_ensembl_meta.values())))} species')
    """
    # Check Ensembl Beta
    if mode in ls_ensembl_mode:
        global ensembl_meta
        ensembl_meta = requests.get(ENSEMBL_TREE_JSON).json()

        dic_ensembl_meta = {}
        for species in ensembl_meta:
            for acc in ensembl_meta[species].keys():
                dic_ensembl_meta[acc] = species

        df[['Ensembl']] = ''
        print(f'  Ensembl Beta : {len(dic_ensembl_meta)} genomes across {len(list(set(dic_ensembl_meta.values())))} species')

    # Check Zoonomia
    if mode in ls_zoonomia_mode:
        # Download the csv bytes
        resp = requests.get(ZOONOMIA_META, verify=False)
        resp.raise_for_status()

        # Read into pandas from a BytesIO buffer
        csv_text = resp.content.decode("utf-8")
        df_zoonomia = pd.read_csv(StringIO(csv_text))

        # add an empty 'Zoonomia' column
        df[['Zoonomia']] = ''
        # now print the TOGA summary
        print(f'  Zoonomia TOGA: {df_zoonomia["Assembly name"].nunique()} genomes across {df_zoonomia["Species"].nunique()} species')

    print('')

    # Check availability of the searched genomes in GenArk and Ensembl
    for idx in df.index:
        refseq_id = df.loc[idx]['RefSeq']
        genbank_id = df.loc[idx]['GenBank']

        if mode in ls_genark_mode and (refseq_id in dic_genark_meta or genbank_id in dic_genark_meta):
            df.loc[idx, 'GenArk'] = 'v'
        if mode in ls_ensembl_mode and (genbank_id in dic_ensembl_meta or refseq_id in dic_ensembl_meta):
            df.loc[idx, 'Ensembl'] = 'v'
        if mode in ls_zoonomia_mode and (genbank_id in df_zoonomia["NCBI accession / source"].tolist() or refseq_id in df_zoonomia["NCBI accession / source"].tolist()):
            df.loc[idx, 'Zoonomia'] = 'v'

    # Remove space in assembly name
    df['Assembly name'] = df['Assembly name'].str.replace(' ', '_', regex=False)

    if mode == 'genome':
        return df, dic_genark_meta, dic_ensembl_meta
    elif mode == 'geneset':
        return df, dic_genark_meta, dic_ensembl_meta, df_zoonomia
    elif mode == 'sequence':
        return df, dic_ensembl_meta
    elif mode == 'annotation':
        return df, dic_genark_meta
    elif mode == 'crossgenome':
        return df, dic_ensembl_meta, df_zoonomia


## ---------------------------------------------
## gencube genome
## ---------------------------------------------
# Download genome data and assembly report
def download_genome (df, types, dic_genark_meta, dic_ensembl_meta, recursive):

    ls_type_raw = types.split(',')

    ls_types = []
    # Change the order priority (refseq -> genbank -> others)
    for type_tmp in ['refseq', 'genbank', 'genark', 'ensembl']:
        if type_tmp in ls_type_raw:
            ls_types.append(type_tmp)

    dic_download = {}
    for idx in df.index:
        assembly_id = df.loc[idx]['Assembly name']
        genbank_id = df.loc[idx]['GenBank']
        refseq_id = df.loc[idx]['RefSeq']
        check_ensembl = df.loc[idx]['Ensembl']
        check_genark = df.loc[idx]['GenArk']
        # Renaming ex. "Homo sapiens (human)" -> "Homo_sapiens"
        organism = re.sub(r'\s*\([^)]*\)', '', df.loc[idx]['Organism']).replace(' ', '_')

        ls_download = []
        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')

        report = False
        for type in ls_types:
            print(f'- {type}')
            # RefSeq
            if refseq_id:
                ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}")
                refseq_assembly_id = ls_source[0]

                refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_assembly_id}'
                refseq_fa = f'{refseq_dir}/{refseq_assembly_id}_genomic.fna.gz'
                refseq_md5sum = f'{refseq_dir}/md5checksums.txt'
                refseq_rp = f'{refseq_dir}/{refseq_assembly_id}_assembly_report.txt'
                out_fa_name = f'{organism}-{assembly_id}-refseq.sm.fa.gz'
                out_rp_name = f'{organism}-{assembly_id}_assembly_report.txt'

                # Check and download assembly report.
                if check_url(refseq_md5sum, show_output=False):

                    if not report:

                        if check_url(refseq_rp):
                            download_genome_url(refseq_rp, out_rp_name, refseq_md5sum, recursive=recursive)
                            report = True
                    # Check and download RefSeq genome.
                    if type == 'refseq':

                        if check_url(refseq_fa):
                            download_genome_url(refseq_fa, out_fa_name, refseq_md5sum, recursive=recursive)
                            ls_download.append('refseq.sm')
                            continue
                        else:
                            print('  !! RefSeq genome is not available. Try to download in GenBank database')
                else:
                    print('  !! RefSeq data is not available. Try to download in GenBank database')

            elif type == 'refseq':
                print('  !! There is no RefSeq accession. Try to download using GenBank accession')

            # GenBank
            if not report or type == 'genbank':
                ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{genbank_id[0:3]}/{genbank_id[4:7]}/{genbank_id[7:10]}/{genbank_id[10:13]}")
                genbank_assembly_id = ls_source[0]

                genbank_dir = f'{NCBI_FTP_URL}/{genbank_id[0:3]}/{genbank_id[4:7]}/{genbank_id[7:10]}/{genbank_id[10:13]}/{genbank_assembly_id}'
                genbank_fa = f'{genbank_dir}/{genbank_assembly_id}_genomic.fna.gz'
                genbank_md5sum = f'{genbank_dir}/md5checksums.txt'
                genbank_rp = f'{genbank_dir}/{genbank_assembly_id}_assembly_report.txt'
                out_fa_name = f'{organism}-{assembly_id}-genbank.sm.fa.gz'
                out_rp_name = f'{organism}-{assembly_id}_assembly_report.txt'

                if check_url(genbank_md5sum, show_output=False):
                    if check_url(genbank_rp):
                        # Check and download assembly report.
                        if not report:
                            download_genome_url(genbank_rp, out_rp_name, genbank_md5sum, recursive=recursive)
                            report = True
                        # Check and download GenBank genome.
                        if type in ['genbank', 'refseq']:
                            if check_url(genbank_fa):
                                download_genome_url(genbank_fa, out_fa_name, genbank_md5sum, recursive=recursive)
                                ls_download.append('genbank.sm')
                                continue
                    else:
                        if not report:
                            print('  !! Assembly_report is not found in RefSeq and GenBank')

                        if type == 'refseq':
                            print('  !! GenBank genome is also not available')
                        elif type == 'genbank':
                                print('  !! GenBank genome is not available')
                else:
                    if type == 'refseq':
                        print('  !! GenBank genome is also not available')
                    elif type == 'genbank':
                        print('  !! GenBank genome is not available')

            # GenArk
            if type == 'genark':
                if check_genark:
                    if refseq_id in dic_genark_meta:
                        genark_id = refseq_id
                    else:
                        genark_id = genbank_id

                    genark_dif = f'{GENARK_URL}/{genark_id[0:3]}/{genark_id[4:7]}/{genark_id[7:10]}/{genark_id[10:13]}/{genark_id}'
                    genark_fa = f'{genark_dif}/{genark_id}.fa.gz'
                    out_fa_name = f'{organism}-{assembly_id}-genark.sm.fa.gz'

                    if check_url(genark_fa):
                        download_genome_url(genark_fa, out_fa_name, recursive=recursive)
                        ls_download.append('genark.sm')
                        continue
                    else:
                        continue

                else:
                    print('  GenArk genome is not available')
                    continue

            # Ensembl Beta
            if type == 'ensembl':
                if check_ensembl:
                    if genbank_id in dic_ensembl_meta:
                        ensembl_acc = genbank_id
                        organism_ens = dic_ensembl_meta[genbank_id].replace(' ', '_')
                    elif refseq_id in dic_ensembl_meta:
                        ensembl_acc = refseq_id
                        organism_ens = dic_ensembl_meta[refseq_id].replace(' ', '_')

                    """ Rapid Release
                    ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}'
                    ls_source = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)
                    source = ''
                    ls_excluded = ['rnaseq', 'brake', 'statistics']
                    if ls_source:
                        for tmp in ls_excluded:
                            if tmp in ls_source:
                                ls_source.remove(tmp)

                    for source in ls_source:
                        ensembl_dir = f'{ENSEMBL_RAPID_FTP_URL}/species/{organism_ens}/{ensembl_acc}/{source}/genome'

                        ensembl_fa = f'{ensembl_dir}/{organism}-{ensembl_acc}-softmasked.fa.gz'
                        ensembl_md5sum = f'{ensembl_dir}/md5sum.txt'
                        out_fa_name = f'{organism}-{assembly_id}-ensembl_{source}.sm.fa.gz'
                    """
                    ensembl_dir = f'{ENSEMBL_BETA_FTP_URL}/{organism_ens}/{ensembl_acc}/genome'

                    ensembl_fa = f'{ensembl_dir}/softmasked.fa.gz'
                    ensembl_md5sum = f'{ensembl_dir}/md5sum.txt'
                    out_fa_name = f'{organism}-{assembly_id}-ensembl.sm.fa.gz'

                    if check_url(ensembl_fa):
                        download_genome_url(ensembl_fa, out_fa_name, ensembl_md5sum, recursive=recursive)
                        ls_download.append('ensembl.sm')
                    else:
                        continue


                else:
                    print('  Ensembl genome is not available')
                    continue
        print('')
        dic_download[genbank_id] = ls_download
    return dic_download

# Change chromosome names of genome file
def convert_chr_label_genome (df, dic_download, style, masking, compresslevel, recursive):

    dic_out_mask = {'soft' : 'soft-masked', 'hard' : 'hard-masked', 'none' : 'unmasked'}
    print(f'# Change chromosome names and masking method: {style}-style & {dic_out_mask[masking]}')

    for idx in df.index:
        assembly_id = df.loc[idx]['Assembly name']
        organism = re.sub(r'\s*\([^)]*\)', '', df.loc[idx]['Organism']).replace(' ', '_')
        genbank_id = df.loc[idx]['GenBank']
        refseq_id = df.loc[idx]['RefSeq']

        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')

        if not os.path.exists(f'{DOWNLOAD_FOLDER_NAME}/{organism}-{assembly_id}_assembly_report.txt'):
            print('  !! Assembly report is not found in RefSeq and GenBank')
            print("     Chromosome names can't be changed\n")
            break

        # Make integrative dataframe of chromosome label of databases
        df_report_edit = make_chr_dataframe (organism, assembly_id)

        # Check the number of UCSC 'na'
        ucsc_na_count = df_report_edit['UCSC'].str.contains('na', na=False).sum()
        # UCSC check
        if style == 'ucsc' and ucsc_na_count > 1:
            print(f'  {ucsc_na_count} chromosome(s) have not ucsc name')
            break

        # Convert the DataFrame to dictionaries for faster lookup
        genbank_to_ensembl = df_report_edit.set_index('GenBank')['Ensembl'].to_dict()
        refseq_to_ensembl = df_report_edit.set_index('RefSeq')['Ensembl'].to_dict()
        genbank_to_gencode = df_report_edit.set_index('GenBank')['Gencode'].to_dict()
        refseq_to_gencode = df_report_edit.set_index('RefSeq')['Gencode'].to_dict()
        genbank_to_ucsc = df_report_edit.set_index('GenBank')['UCSC'].to_dict()
        refseq_to_ucsc = df_report_edit.set_index('RefSeq')['UCSC'].to_dict()
        ensembl_to_gencode = df_report_edit.set_index('Ensembl')['Gencode'].to_dict()
        ensembl_to_ucsc = df_report_edit.set_index('Ensembl')['UCSC'].to_dict()

        # Check downloaded files
        ls_download = dic_download[genbank_id]
        print(f'  Downloaded file(s): {ls_download}')
        for db in ls_download:
            print(f'  - {db}')

            dic_out_mask = {'soft' : 'sm', 'hard' : 'hm', 'none' : 'um'}
            dic_out_suffix = {'ensembl' : 'ens-id', 'gencode' : 'gc-id', 'ucsc' : 'ucsc-id'}

            # Input and output name
            db_name = db.split('.')[0]

            in_file = f'{organism}-{assembly_id}-{db}.fa.gz'
            out_file = f'{organism}-{assembly_id}-{db_name}.{dic_out_mask[masking]}.{dic_out_suffix[style]}.fa.gz'

            # Remove output file (--reculsive)
            if recursive:
                delete_file(f'{out_subfolder}/{out_file}')

            if 'ensembl' in db_name and style == 'ensembl' and masking == 'soft':
                if not os.path.exists(f'{out_subfolder}/{out_file}'):
                    print('    The downloaded file already contains Ensembl-style names and is soft-masked')
                    print('    Just copy to output folder')
                    subprocess.run(['cp', f'{DOWNLOAD_FOLDER_NAME}/{in_file}', f'{out_subfolder}/{out_file}'], check=True)
                    continue
                else:
                    print('    The converted file already exists')
                    continue

            # File check in working directory
            ls_download_files = os.listdir(DOWNLOAD_FOLDER_NAME)
            ls_output_folder_files = os.listdir(out_subfolder)
            ls_read = []
            ls_write = []
            if in_file in ls_download_files and out_file not in ls_output_folder_files:

                start_time = time.time() # record start times

                print('    Modify chromosome names')
                try:
                    with gzip.open(f'{DOWNLOAD_FOLDER_NAME}/{in_file}', 'rt') as f_in:
                        for line in f_in:
                            ls_read.append(line)
                except gzip.BadGzipFile:
                    print('    !! The file is not a valid gzip file\n')

                for line in ls_read:
                    if line.startswith('>'):
                        # Chromosome name and extra
                        chr_name = line.strip().split()[0][1:]
                        extra = ' '.join(line.strip().split()[1:])
                        # Edit masking information in extra part (only Ensembl)
                        if 'ensembl' in db_name:
                            if masking == 'hard':
                                extra = extra.replace('softmasked', 'hardmasked')
                            elif masking == 'none':
                                extra = extra.replace('softmasked', 'unmasked')

                        # Perform the lookup based on the db and style
                        if db_name in ['genbank', 'refseq', 'genark']:
                            if style == 'ensembl':
                                chage_chr_name = genbank_to_ensembl.get(chr_name, refseq_to_ensembl.get(chr_name, chr_name))
                            elif style == 'gencode':
                                chage_chr_name = genbank_to_gencode.get(chr_name, refseq_to_gencode.get(chr_name, chr_name))
                            elif style == 'ucsc':
                                chage_chr_name = genbank_to_ucsc.get(chr_name, refseq_to_ucsc.get(chr_name, chr_name))
                        elif 'ensembl' in db_name:
                            if style == 'gencode':
                                chage_chr_name = ensembl_to_gencode.get(chr_name, chr_name)
                            elif style == 'ucsc':
                                chage_chr_name = ensembl_to_ucsc.get(chr_name, chr_name)

                        if 'ensembl' in db_name and style == 'ensembl':
                            ls_write.append(f'>{chr_name} {extra}\n')
                        else:
                            ls_write.append(f'>{chage_chr_name} {extra}\n')
                    else:
                        # Convert masking method in sequence information
                        # - Hard-masking
                        if masking == 'hard':
                            line_edit = ''
                            for char in line:
                                if char.islower():
                                    line_edit += 'N'
                                else:
                                    line_edit += char

                            ls_write.append(line_edit)
                        # - Unmasking
                        elif masking == 'none':
                            ls_write.append(line.upper())

                        # - Soft-masking
                        elif masking == 'soft':
                            ls_write.append(line)


                # Write and compressed fasta file
                print(f'    Write compressed fasta file (compresslevel: {compresslevel})')
                with gzip.open(f'{out_subfolder}/{out_file}', 'wt', compresslevel=int(compresslevel)) as f_out:
                    for line in ls_write:
                        f_out.write(line)

                end_time = time.time() # record end time
                elapsed_time = end_time - start_time
                print(f'    Processing time: {int(elapsed_time)} seconds')

            elif out_file in ls_output_folder_files:
                print('    The converted file already exists')
        print('\n  !! If the file appears to have any problems, please delete it and retry the process\n')


## ---------------------------------------------
## gencube geneset
## ---------------------------------------------
# Check full accessibility
def process_row_geneset(idx, row, dic_genark_meta, dic_ensembl_meta, df_zoonomia):
    result = {
        'index': idx,
        'RefSeq': '',
        'GenArk': '',
        'Ensembl': '',
        'Zoonomia': ''
    }
    genbank_id = row['GenBank']
    refseq_id = row['RefSeq']
    check_ensembl = row['Ensembl']
    check_genark = row['GenArk']
    check_zoonomia = row['Zoonomia']

    # RefSeq
    if refseq_id:
        ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}")
        refseq_assembly_id = ls_source[0]

        refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_assembly_id}'
        url_md5sum = f'{refseq_dir}/md5checksums.txt'
        refseq_gtf = f'{refseq_assembly_id}_genomic.gtf.gz'
        refseq_gff = f'{refseq_assembly_id}_genomic.gff.gz'
        refseq_genomon = f'{refseq_assembly_id}_gnomon_model.gff.gz'
        refseq_cross = f'{refseq_assembly_id}_cross_species_tx_alns.gff.gz'
        refseq_same = f'{refseq_assembly_id}_same_species_tx_alns.gff.gz'

        if check_url(url_md5sum, show_output=False):
            df_md5 = get_md5(url_md5sum)  # Read md5sum information
            ls_filename = df_md5['Filename'].str.split('/').str[-1].tolist()

            if refseq_gtf in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'gtf')
            if refseq_gff in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'gff')
            if refseq_genomon in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'gnomon')
            if refseq_cross in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'cross')
            if refseq_same in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'same')

    # Genark
    if check_genark:
        genark_id = refseq_id if refseq_id in dic_genark_meta else genbank_id
        genark_dir = f'{GENARK_URL}/{genark_id[0:3]}/{genark_id[4:7]}/{genark_id[7:10]}/{genark_id[10:13]}/{genark_id}'

        ls_gene_files = list_http_files(genark_dir + '/genes')
        for file in ls_gene_files:
            if genark_id in file:
                type = file.split('.')[-3]
                if type == 'augustus':
                    result['GenArk'] = add_string(result['GenArk'], 'augustus')
                if type == 'xenoRefGene':
                    result['GenArk'] = add_string(result['GenArk'], 'xenoref')
                if type == 'ncbiRefSeq':
                    result['GenArk'] = add_string(result['GenArk'], 'ref')

    # Ensembl Beta
    if check_ensembl:
        if genbank_id in dic_ensembl_meta:
            ensembl_acc = genbank_id
            organism_ens = dic_ensembl_meta[genbank_id].replace(' ', '_')
        elif refseq_id in dic_ensembl_meta:
            ensembl_acc = refseq_id
            organism_ens = dic_ensembl_meta[refseq_id].replace(' ', '_')

        """ Rapid Release
        ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}'
        ls_source = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)

        source = ''
        ls_excluded = ['rnaseq', 'brake', 'statistics']
        if ls_source:
            for tmp in ls_source:
                if tmp not in ls_excluded:
                    if not source:
                        source = tmp
                    else:
                        source += f', {tmp}'
        """
        # Check source directories
        ls_source = list(ensembl_meta[organism_ens][ensembl_acc].keys())

        source = ''
        ls_excluded = ['genome']
        if ls_source:
            for tmp in ls_excluded:
                if tmp in ls_source:
                    ls_source.remove(tmp)

        ls_input = []
        for source in ls_source:
            ls_tmp = [source]
            # Check the geneset folder name
            genomebuild = list(ensembl_meta[organism_ens][ensembl_acc][source]['geneset'].keys())[0]

            ls_filename = list(ensembl_meta[organism_ens][ensembl_acc][source]['geneset'][genomebuild].keys())
            ensembl_gtf = 'genes.gtf.gz'
            ensembl_gff = 'genes.gff3'

            if ensembl_gtf in ls_filename:
                ls_tmp.append('gtf')
            if ensembl_gff in ls_filename:
                ls_tmp.append('gff')

            ls_input.append(ls_tmp)

        if ls_input:
            str_input = ''
            for i in range(len(ls_input)):
                if i != 0:
                    str_input += ' / '
                for j in range(len(ls_input[i])):
                    if j == 0:
                        str_input += f'{ls_input[i][j]}:'
                    else:
                        str_input += f' {ls_input[i][j]}'
                    if j != 0 and j != len(ls_input[i]) - 1:
                        str_input += ','

            result['Ensembl'] = add_string(result['Ensembl'], str_input)

        """
        # Check source directories
        ls_source = list(ensembl_meta[organism_ens][ensembl_acc].keys())

        source = ''
        if ls_source:
            for tmp in ls_source:
                if tmp != 'genome':
                    if not source:
                        source = tmp
                    else:
                        source += f', {tmp}'

        result['Ensembl'] = add_string(result['Ensembl'], source)
        """

    # Zoonomia
    if check_zoonomia:
        ls_reference = df_zoonomia[
            (df_zoonomia["NCBI accession / source"] == genbank_id) |
            (df_zoonomia["NCBI accession / source"] == refseq_id)
        ]['reference'].tolist()
        reference = ', '.join(list(set(ls_reference)))

        result['Zoonomia'] = add_string(result['Zoonomia'], reference)

    return result

def check_access_full_geneset(df, dic_genark_meta, dic_ensembl_meta, df_zoonomia):
    df_full_annotation = df[['Assembly name', 'UCSC']].copy()
    ls_geneset = ['RefSeq', 'GenArk', 'Ensembl', 'Zoonomia']

    for label in ls_geneset:
        df_full_annotation[label] = ''

    print('# Check accessible data in databases')

    results = []
    with ThreadPoolExecutor() as executor:
        futures = []
        for idx in df.index:
            row = df.loc[idx]
            futures.append(executor.submit(process_row_geneset, idx, row, dic_genark_meta, dic_ensembl_meta, df_zoonomia))

        for future in futures:
            results.append(future.result())

    for result in results:
        idx = result.pop('index')
        for label in ls_geneset:
            df_full_annotation.loc[idx, label] = result[label]

    # Searched result
    print(tabulate(df_full_annotation, headers='keys', tablefmt='grid'))
    print('')

    return df_full_annotation

# Download geneset data
def download_geneset(df, df_genome, dic_ensembl_meta, dic_genark_meta, df_zoonomia, types, recursive):

    ls_types = types.split(',')

    print('# Download geneset data')
    dic_download = {}
    for idx in df.index:
        # Accession or name
        assembly_id = df_genome.loc[idx]['Assembly name']
        genbank_id = df_genome.loc[idx]['GenBank']
        refseq_id = df_genome.loc[idx]['RefSeq']
        organism = re.sub(r'\s*\([^)]*\)', '', df_genome.loc[idx]['Organism']).replace(' ', '_')

        check_refseq = df.loc[idx]['RefSeq']
        check_genark = df.loc[idx]['GenArk']
        check_ensembl = df.loc[idx]['Ensembl']
        check_zoonomia = df.loc[idx]['Zoonomia']

        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')

        ls_download = []
        # RefSeq
        if len(list(set(['refseq_gtf', 'refseq_gff', 'gnomon', 'cross_species', 'same_species']) & set(ls_types))) > 0:

            ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}")
            refseq_assembly_id = ls_source[0]
            refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_assembly_id}'

            ls_search = check_refseq.replace(' ', '').split(',')
            url_md5sum = f'{refseq_dir}/md5checksums.txt'
            refseq_gtf = f'{refseq_dir}/{refseq_assembly_id}_genomic.gtf.gz'
            refseq_gff = f'{refseq_dir}/{refseq_assembly_id}_genomic.gff.gz'
            refseq_genomon = f'{refseq_dir}/Gnomon_models/{refseq_assembly_id}_gnomon_model.gff.gz'
            refseq_cross = f'{refseq_dir}/Evidence_alignments/{refseq_assembly_id}_cross_species_tx_alns.gff.gz'
            refseq_same = f'{refseq_dir}/Evidence_alignments/{refseq_assembly_id}_same_species_tx_alns.gff.gz'

            if check_refseq:
                if 'refseq_gtf' in ls_types:
                    if 'gtf' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq.gtf.gz'
                        download_url(refseq_gtf, out_name, url_md5sum=url_md5sum, recursive=recursive)
                        ls_download.append('refseq_gtf')

                if 'refseq_gff' in ls_types:
                    if 'gff' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq.gff.gz'
                        download_url(refseq_gff, out_name, url_md5sum=url_md5sum, recursive=recursive)
                        ls_download.append('refseq_gff')

                if 'gnomon' in ls_types:
                    if 'gnomon' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq_genomon.gff.gz'
                        download_url(refseq_genomon, out_name, url_md5sum=url_md5sum, recursive=recursive)
                        ls_download.append('gnomon')

                if 'cross' in ls_types:
                    if 'cross' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq_cross.gff.gz'
                        download_url(refseq_cross, out_name, url_md5sum=url_md5sum, recursive=recursive)
                        ls_download.append('cross_species')

                if 'same' in ls_types:
                    if 'same' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq_same.gff.gz'
                        download_url(refseq_same, out_name, url_md5sum=url_md5sum, recursive=recursive)
                        ls_download.append('same_species')

        # Genark
        if len(list(set(['augustus', 'xenoref', 'genark_ref']) & set(ls_types))) > 0:

            if refseq_id in dic_genark_meta:
                genark_id = refseq_id
            else:
                genark_id = genbank_id
            genark_dif = f'{GENARK_URL}/{genark_id[0:3]}/{genark_id[4:7]}/{genark_id[7:10]}/{genark_id[10:13]}/{genark_id}'

            ls_search = check_genark.replace(' ', '').split(',')
            genark_augustus = f'{genark_dif}/genes/{genark_id}_{dic_genark_meta[genark_id]}.augustus.gtf.gz'
            genark_xeno = f'{genark_dif}/genes/{genark_id}_{dic_genark_meta[genark_id]}.xenoRefGene.gtf.gz'
            genark_ref = f'{genark_dif}/genes/{genark_id}_{dic_genark_meta[genark_id]}.ncbiRefSeq.gtf.gz'

            if check_genark:

                if 'augustus' in ls_types:
                    if 'augustus' in ls_search:
                        out_name = f'{organism}-{assembly_id}-genark_augustus.gtf.gz'
                        download_url(genark_augustus, out_name, verify=False, recursive=recursive)
                        ls_download.append('augustus')

                if 'xenoref' in ls_types:
                    if 'xenoref' in ls_search:
                        out_name = f'{organism}-{assembly_id}-genark_xenoref.gtf.gz'
                        download_url(genark_xeno, out_name, verify=False, recursive=recursive)
                        ls_download.append('xenoref')

                if 'genark_ref' in ls_types:
                    if 'ref' in ls_search:
                        out_name = f'{organism}-{assembly_id}-genark_refseq.gtf.gz'
                        download_url(genark_ref, out_name, verify=False, recursive=recursive)
                        ls_download.append('genark_ref')

        # Ensembl Beta
        if len(list(set(['ensembl_gtf', 'ensembl_gff']) & set(ls_types))) > 0:

            if genbank_id in dic_ensembl_meta:
                ensembl_acc = genbank_id
                organism_ens = dic_ensembl_meta[genbank_id].replace(' ', '_')
            elif refseq_id in dic_ensembl_meta:
                ensembl_acc = refseq_id
                organism_ens = dic_ensembl_meta[refseq_id].replace(' ', '_')

            #if check_ensembl:
            #    ls_source = check_ensembl.replace(' ', '').split(',')
            # 지울 예정

            if check_ensembl:
                ls_source = []
                for search_in_source in check_ensembl.replace(' ', '').split('/'):
                    source = search_in_source.strip().split(':')[0]
                    ls_source.append(source)
                    ls_search = search_in_source.split(':')[1].split(',')

                """ Rapid Release
                for source in ls_source:
                    ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}/{source}/geneset'
                    # Check the geneset folder name
                    geneset = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)[0]

                    ensembl_file_dir = f'{ENSEMBL_RAPID_FTP_URL}/species/{organism_ens}/{ensembl_acc}/{source}/geneset/{geneset}'

                    url_md5sum = f'{ensembl_file_dir}/md5sum.txt'
                    ensembl_gtf = f'{ensembl_file_dir}/{organism_ens}-{ensembl_acc}-{geneset}-genes.gtf.gz'
                    ensembl_gff = f'{ensembl_file_dir}/{organism_ens}-{ensembl_acc}-{geneset}-genes.gff3.gz'

                    if 'ensembl_gtf' in ls_types:
                        out_name = f'{organism}-{assembly_id}-ensembl_{source}.gtf.gz'
                        if check_url(ensembl_gtf, file_name=out_name):
                            download_url(ensembl_gtf, out_name, url_md5sum=url_md5sum)
                            ls_download.append(f'ensembl-{source}-gtf')

                    if 'ensembl_gff' in ls_types:
                        out_name = f'{organism}-{assembly_id}-ensembl_{source}.gff.gz'
                        if check_url(ensembl_gff, file_name=out_name):
                            download_url(ensembl_gff, out_name, url_md5sum=url_md5sum)
                            ls_download.append(f'ensembl-{source}-gff')
                """
                for source in ls_source:
                    # Check the geneset folder name
                    genomebuild = list(ensembl_meta[organism_ens][ensembl_acc][source]['geneset'].keys())[0]
                    ensembl_file_dir = f'{ENSEMBL_BETA_FTP_URL}/{organism_ens}/{ensembl_acc}/{source}/geneset/{genomebuild}'
                    # url_md5sum = f'{ensembl_file_dir}/md5sum.txt'
                    ensembl_gtf = f'{ensembl_file_dir}/genes.gtf.gz'
                    ensembl_gff = f'{ensembl_file_dir}/genes.gff3'

                    if 'ensembl_gtf' in ls_types:
                        out_name = f'{organism}-{assembly_id}-ensembl_{source}.gtf.gz'
                        if check_url(ensembl_gtf, file_name=out_name):
                            download_url(ensembl_gtf, out_name, recursive=recursive)
                            ls_download.append(f'ensembl-{source}-gtf')

                    if 'ensembl_gff' in ls_types:
                        out_name = f'{organism}-{assembly_id}-ensembl_{source}.gff'
                        if check_url(ensembl_gff, file_name=out_name):
                            download_url(ensembl_gff, out_name, recursive=recursive)
                            ls_download.append(f'ensembl-{source}-gff')

        # Zoonomia
        if len(list(set(['toga_gtf', 'toga_bed', 'toga_pseudo']) & set(ls_types))) > 0:

            if check_zoonomia:
                ls_reference = check_zoonomia.replace(' ', '').split(',')

                for reference in ls_reference:

                    zoonomia_dir = f'{ZOONOMIA_URL}/{DIC_ZOONOMIA[reference]}'

                    df_tmp = df_zoonomia[
                        ((df_zoonomia['NCBI accession / source'] == genbank_id) |
                        (df_zoonomia['NCBI accession / source'] == refseq_id)) &
                        (df_zoonomia['reference'] == reference)
                    ]

                    for i in range(len(df_tmp.index)):
                        taxo = df_tmp['Taxonomic Lineage'].values[i]
                        species = df_tmp['Species'].values[i].replace(' ', '_')
                        name = df_tmp['Common name'].values[i].replace(' ', '_')
                        assembly = df_tmp['Assembly name'].values[i]

                        if reference in ['human', 'mouse', 'chicken']:
                            ls_folders = list_http_folders(zoonomia_dir)

                            for folder in ls_folders:
                                if folder in taxo:
                                    category = folder
                                    break

                            zoonomia_file_dir = f'{zoonomia_dir}/{category}/{species}__{name}__{assembly}'
                        else:
                            zoonomia_file_dir = f'{zoonomia_dir}/{species}__{name}__{assembly}'

                        zoonomia_gtf = f'{zoonomia_file_dir}/geneAnnotation.gtf.gz'
                        zoonomia_bed = f'{zoonomia_file_dir}/geneAnnotation.bed.gz'
                        zoonomia_pseudo = f'{zoonomia_file_dir}/processedPseudogeneAnnotation.bed.gz'

                        if 'toga_gtf' in ls_types:
                            out_name = f'{organism}-{assembly_id}-toga_{reference}.gtf.gz'
                            if check_url(zoonomia_gtf, verify=False, file_name=out_name):
                                download_url(zoonomia_gtf, out_name, verify=False, recursive=recursive)
                                ls_download.append(f'zoonomia-{reference}-gtf')

                        if 'toga_bed' in ls_types:
                            out_name = f'{organism}-{assembly_id}-toga_{reference}.bed.gz'
                            if check_url(zoonomia_bed, verify=False, file_name=out_name):
                                download_url(zoonomia_bed, out_name, verify=False, recursive=recursive)
                                ls_download.append(f'zoonomia-{reference}-bed')

                        if 'toga_pseudo' in ls_types:
                            out_name = f'{organism}-{assembly_id}-toga_{reference}_pseudo.bed.gz'
                            if check_url(zoonomia_pseudo, verify=False, file_name=out_name):
                                download_url(zoonomia_pseudo, out_name, verify=False, recursive=recursive)
                                ls_download.append(f'zoonomia-{reference}-pseudobed')

        dic_download[genbank_id] = ls_download
        print('')
    return dic_download

# Change chromosome names of geneset file
def convert_chr_label_geneset (df, dic_download, style, recursive):

    print(f'# Change chromosome names: {style}-style')

    for idx in df.index:
        assembly_id = df.loc[idx]['Assembly name']
        organism = re.sub(r'\s*\([^)]*\)', '', df.loc[idx]['Organism']).replace(' ', '_')
        genbank_id = df.loc[idx]['GenBank']
        refseq_id = df.loc[idx]['RefSeq']

        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')

        # Download assembly report
        print('  Download assembly report')
        report = False
        # RefSeq
        if refseq_id:
            ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}")
            refseq_assembly_id = ls_source[0]
            refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_assembly_id}'
            refseq_md5sum = f'{refseq_dir}/md5checksums.txt'
            refseq_rp = f'{refseq_dir}/{refseq_assembly_id}_assembly_report.txt'
            out_rp_name = f'{organism}-{assembly_id}_assembly_report.txt'
            # Check and download assembly report.

            if check_url(refseq_rp):
                download_genome_url(refseq_rp, out_rp_name, refseq_md5sum)
                report = True

        # GenBank
        if not report:
            ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{genbank_id[0:3]}/{genbank_id[4:7]}/{genbank_id[7:10]}/{genbank_id[10:13]}")
            genbank_assembly_id = ls_source[0]
            genbank_dir = f'{NCBI_FTP_URL}/{genbank_id[0:3]}/{genbank_id[4:7]}/{genbank_id[7:10]}/{genbank_id[10:13]}/{genbank_assembly_id}'
            genbank_md5sum = f'{genbank_dir}/md5checksums.txt'
            genbank_rp = f'{genbank_dir}/{genbank_assembly_id}_assembly_report.txt'
            out_rp_name = f'{organism}-{assembly_id}_assembly_report.txt'

            if check_url(genbank_rp):
                # Check and download assembly report.
                download_genome_url(genbank_rp, out_rp_name, genbank_md5sum)

            else:
                print('  !! Assembly report is not found in RefSeq and GenBank')
                print("     Chromosome names can't be changed\n")
                break

        # Make integrative dataframe of chromosome label of databases
        df_report_edit = make_chr_dataframe (organism, assembly_id)

        # Check the number of UCSC 'na'
        ucsc_na_count = df_report_edit['UCSC'].str.contains('na', na=False).sum()
        # UCSC check
        if style == 'ucsc' and ucsc_na_count > 1:
            print(f'  {ucsc_na_count} chromosome(s) have not ucsc name')
            break

        df_report_edit['genbank_noversion'] = df_report_edit['GenBank'].apply(lambda x: x.split('.')[0])
        df_report_edit['refseq_noversion'] = df_report_edit['RefSeq'].apply(lambda x: x.split('.')[0])

        # Convert the DataFrame to dictionaries for faster lookup
        genbank_to_ensembl = df_report_edit.set_index('GenBank')['Ensembl'].to_dict()
        refseq_to_ensembl = df_report_edit.set_index('RefSeq')['Ensembl'].to_dict()
        genbank_to_gencode = df_report_edit.set_index('GenBank')['Gencode'].to_dict()
        refseq_to_gencode = df_report_edit.set_index('RefSeq')['Gencode'].to_dict()
        genbank_to_ucsc = df_report_edit.set_index('GenBank')['UCSC'].to_dict()
        refseq_to_ucsc = df_report_edit.set_index('RefSeq')['UCSC'].to_dict()
        ensembl_to_gencode = df_report_edit.set_index('Ensembl')['Gencode'].to_dict()
        ensembl_to_ucsc = df_report_edit.set_index('Ensembl')['UCSC'].to_dict()
        ucsc_to_ensembl = df_report_edit.set_index('UCSC')['Ensembl'].to_dict()
        ucsc_to_gencode = df_report_edit.set_index('UCSC')['Gencode'].to_dict()
        # For Zoonomia data
        genbank_noversion_to_ensembl = df_report_edit.set_index('genbank_noversion')['Ensembl'].to_dict()
        refseq_noversion_to_ensembl = df_report_edit.set_index('refseq_noversion')['Ensembl'].to_dict()
        genbank_noversion_to_gencode = df_report_edit.set_index('genbank_noversion')['Gencode'].to_dict()
        refseq_noversion_to_gencode = df_report_edit.set_index('refseq_noversion')['Gencode'].to_dict()
        genbank_noversion_to_ucsc = df_report_edit.set_index('genbank_noversion')['UCSC'].to_dict()
        refseq_noversion_to_ucsc = df_report_edit.set_index('refseq_noversion')['UCSC'].to_dict()

        # Check the number of UCSC 'na'
        ucsc_na_count = df_report_edit['UCSC'].str.contains('na', na=False).sum()
        # UCSC check
        if style == 'ucsc' and ucsc_na_count > 1:
            print(f'  {ucsc_na_count} chromosome(s) have not ucsc name')
            break

        # Check downloaded files
        ls_download = dic_download[genbank_id]
        print(f'  Downloaded file(s): {ls_download}')


        dic_db_suffix = {
            'refseq_gtf' : 'refseq.gtf.gz',
            'refseq_gff' : 'refseq.gff.gz',
            'gnomon' : 'refseq_genomon.gff.gz',
            'cross_species' : 'refseq_cross.gff.gz',
            'same_species' : 'refseq_same.gff.gz',
            'augustus' : 'genark_augustus.gtf.gz',
            'xenoref' : 'genark_xenoref.gtf.gz',
            'genark_ref' : 'genark_refseq.gtf.gz',
        }
        dic_out_suffix = {'ensembl' : 'ens-id', 'gencode' : 'gc-id', 'ucsc' : 'ucsc-id'}

        for db in ls_download:
            print(f'  - {db}')

            # Input and output name
            if 'ensembl' in db:
                ls_tmp = db.split('-')
                if ls_tmp[2] == 'gtf':
                    in_file = f'{organism}-{assembly_id}-ensembl_{ls_tmp[1]}.gtf.gz'
                    out_file = f'{organism}-{assembly_id}-ensembl_{ls_tmp[1]}.{dic_out_suffix[style]}.gtf'
                elif ls_tmp[2] == 'gff':
                    in_file = f'{organism}-{assembly_id}-ensembl_{ls_tmp[1]}.gff'
                    out_file = f'{organism}-{assembly_id}-ensembl_{ls_tmp[1]}.{dic_out_suffix[style]}.gff'
            elif 'zoonomia' in db:
                ls_tmp = db.split('-')
                if ls_tmp[2] == 'gtf':
                    in_file = f'{organism}-{assembly_id}-toga_{ls_tmp[1]}.gtf.gz'
                    out_file = f'{organism}-{assembly_id}-toga_{ls_tmp[1]}.{dic_out_suffix[style]}.gtf'
                elif ls_tmp[2] == 'bed':
                    in_file = f'{organism}-{assembly_id}-toga_{ls_tmp[1]}.bed.gz'
                    out_file = f'{organism}-{assembly_id}-toga_{ls_tmp[1]}.{dic_out_suffix[style]}.bed'
                elif ls_tmp[2] == 'pseudobed':
                    in_file = f'{organism}-{assembly_id}-toga_{ls_tmp[1]}_pseudo.bed.gz'
                    out_file = f'{organism}-{assembly_id}-toga_{ls_tmp[1]}_pseudo.{dic_out_suffix[style]}.bed'
            else:
                db_suffix = dic_db_suffix[db].split('.')

                in_file = f'{organism}-{assembly_id}-{dic_db_suffix[db]}'
                out_file = f'{organism}-{assembly_id}-{db_suffix[0]}.{dic_out_suffix[style]}.{db_suffix[1]}'

            # Remove output file (--reculsive)
            if recursive:
                delete_file(f'{out_subfolder}/{out_file}')

            if db in ['ensembl_gtf', 'ensembl_gff'] and style == 'ensembl':
                # print('  !! The file downloaded from the Ensembl database already has ensembl-style chromosome names')
                continue

            # File check in working directory
            ls_download_files = os.listdir(DOWNLOAD_FOLDER_NAME)
            ls_output_folder_files = os.listdir(out_subfolder)

            print('    Modify chromosome names')
            if in_file in ls_download_files and out_file not in ls_output_folder_files:

                ls_read = []
                ls_write = []
                try:
                    if in_file.split('.')[-1] == 'gz':
                        with gzip.open(f'{DOWNLOAD_FOLDER_NAME}/{in_file}', 'rt') as f_in:
                            for line in f_in:
                                ls_read.append(line)
                    else:
                        with open(f'{DOWNLOAD_FOLDER_NAME}/{in_file}', 'r') as f_in:
                            for line in f_in:
                                ls_read.append(line)

                except gzip.BadGzipFile:
                    print('  !! The file is not a valid gzip file\n')

                count = 0
                for line in ls_read:
                    if not line.startswith('#'):
                        ls_tmp = line.strip().split('\t')
                        chr_name = ls_tmp[0]

                        # Remove 'gene' block in refeq geneset data
                        category = ls_tmp[2]
                        if category == 'gene' or category == 'region':
                                continue

                        # Change the chromosome name
                        if db in ['refseq_gtf', 'refseq_gff', 'gnomon', 'cross_species', 'same_species', 'augustus', 'xenoref', 'genark_ref']:
                            if style == 'ensembl':
                                ls_tmp[0] = genbank_to_ensembl.get(chr_name, refseq_to_ensembl.get(chr_name, chr_name))
                            elif style == 'gencode':
                                ls_tmp[0] = genbank_to_gencode.get(chr_name, refseq_to_gencode.get(chr_name, chr_name))
                            elif style == 'ucsc':
                                ls_tmp[0] = genbank_to_ucsc.get(chr_name, refseq_to_ucsc.get(chr_name, chr_name))
                        elif 'ensembl' in db:
                            if style == 'gencode':
                                ls_tmp[0] = ensembl_to_gencode.get(chr_name, chr_name)
                            elif style == 'ucsc':
                                ls_tmp[0] = ensembl_to_ucsc.get(chr_name, chr_name)
                        elif 'zoonomia' in db:
                            if style == 'ensembl':
                                ls_tmp[0] = ucsc_to_ensembl.get(chr_name, genbank_noversion_to_ensembl.get(chr_name, refseq_noversion_to_ensembl.get(chr_name, chr_name)))
                            elif style == 'gencode':
                                ls_tmp[0] = ucsc_to_gencode.get(chr_name, genbank_noversion_to_gencode.get(chr_name, refseq_noversion_to_gencode.get(chr_name, chr_name)))
                            elif style == 'ucsc':
                                ls_tmp[0] = genbank_noversion_to_ucsc.get(chr_name, refseq_noversion_to_ucsc.get(chr_name, chr_name))

                        line_edit = '\t'.join(ls_tmp) + '\n'
                        ls_write.append(line_edit)

                        count += 1
                    else:
                        ls_write.append(line)

                format = out_file.split('.')[-1]
                # bed format
                if format == 'bed':
                    # Write gtf, gff, and bed file
                    with open(f'{out_subfolder}/{out_file}_tmp', 'w') as f_out:
                        for line in ls_write:
                            f_out.write(line)

                    # Sort file
                    print('    Sort file')
                    with open(f'{out_subfolder}/{out_file}', 'w') as f_out:
                        subprocess.run(['sort', '-k1,1', '-k2,2n', f'{out_subfolder}/{out_file}_tmp'], stdout=f_out, check=True)

                    # Remove temporary file
                    delete_file(f'{out_subfolder}/{out_file}_tmp')

                    # Compress file
                    subprocess.run(['gzip', '--best', '--force', f'{out_subfolder}/{out_file}'], check=True)

                # gtf or gff format
                else:
                    # Write and compressed gtf, gff, and bed file
                    with gzip.open(f'{out_subfolder}/{out_file}.gz', 'wt', compresslevel=9) as f_out:
                        for line in ls_write:
                            f_out.write(line)

                """
                # gtf or gff format
                else:
                    with open(f'{OUT_FOLDER_NAME}/{out_file}', 'w') as f_out:
                        subprocess.run(['sort', '-k1,1', '-k4,4n', f'{OUT_FOLDER_NAME}/{out_file}_tmp'], stdout=f_out, check=True)

                # Make index file
                print(f'  Make index file: {OUT_FOLDER_NAME}/{out_file}.gz')
                if format == 'bed': # bed format
                    #bed = pybedtools.BedTool(f'{OUT_FOLDER_NAME}/{out_file}.gz')
                    pysam.tabix_index(f'{OUT_FOLDER_NAME}/{out_file}.gz', force=True, preset='bed')
                elif format =='gtf':
                    #gtf = pybedtools.BedTool(f'{OUT_FOLDER_NAME}/{out_file}.gz')
                    pysam.tabix_index(f'{OUT_FOLDER_NAME}/{out_file}.gz', force=True, preset='gff')
                elif format == 'gff':
                    #gff = pybedtools.BedTool(f'{OUT_FOLDER_NAME}/{out_file}.gz')
                    pysam.tabix_index(f'{OUT_FOLDER_NAME}/{out_file}.gz', force=True, preset='gff')
                else:
                    print(f'debug point: format == {format}')
                """

                print('')

            else:
                print('  The output file already exists')
                print('  !! If the file appears to have any problems, please delete it and retry the process\n')


## ---------------------------------------------
## gencube annotation
## ---------------------------------------------
# Check full accessibility
def process_row_annotation(idx, row, dic_genark_meta):
    result = {
        'index': idx,
        'GenArk': ''
    }

    genbank_id = row['GenBank']
    refseq_id = row['RefSeq']
    check_genark = row['GenArk']

    # GenArk
    if check_genark:
        genark_id = refseq_id if refseq_id in dic_genark_meta else genbank_id
        genark_dir = f'{GENARK_URL}/{genark_id[0:3]}/{genark_id[4:7]}/{genark_id[7:10]}/{genark_id[10:13]}/{genark_id}'

        ls_bb_files = list_http_files(genark_dir + '/bbi')
        gap = False
        rmsk = False
        cpg_island = False
        for file in ls_bb_files:
            if genark_id in file:
                type = file.split('.')[-2]
                if type == 'gc5Base':
                    result['GenArk'] = add_string(result['GenArk'], 'gc')
                if type == 'simpleRepeat':
                    result['GenArk'] = add_string(result['GenArk'], 'sr')
                if type == 'tandemDups':
                    result['GenArk'] = add_string(result['GenArk'], 'td')
                if type == 'windowMasker':
                    result['GenArk'] = add_string(result['GenArk'], 'wm')
                if ('.allGaps.bb' in file or '.gap.bb' in file) and not gap:
                    result['GenArk'] = add_string(result['GenArk'], 'gap')
                    gap = True
                if '.rmsk.' in file and not rmsk:
                    result['GenArk'] = add_string(result['GenArk'], 'rm')
                    rmsk = True
                if '.cpgIslandExt' in file and not cpg_island:
                    result['GenArk'] = add_string(result['GenArk'], 'cpgisland')
                    cpg_island = True

    return result

def check_access_full_annotation(df, dic_genark_meta):
    df_full_annotation = df[['Assembly name', 'UCSC']].copy()
    ls_annotation = ['GenArk']

    for label in ls_annotation:
        df_full_annotation[label] = ''

    print('# Check accessible data in databases')

    results = []
    with ThreadPoolExecutor() as executor:
        futures = []
        for idx in df.index:
            row = df.loc[idx]
            futures.append(executor.submit(process_row_annotation, idx, row, dic_genark_meta))

        for future in futures:
            results.append(future.result())

    for result in results:
        idx = result.pop('index')
        for label in ls_annotation:
            df_full_annotation.loc[idx, label] = result[label]

    # Searched result
    print(tabulate(df_full_annotation, headers='keys', tablefmt='grid'))
    print('')

    return df_full_annotation

# Download annotation data
def download_annotation(df, df_genome, dic_genark_meta, types, recursive):

    ls_types = types.split(',')

    print('# Download annotation data')
    dic_download = {}
    for idx in df.index:
        # Accession or name
        assembly_id = df_genome.loc[idx]['Assembly name']
        genbank_id = df_genome.loc[idx]['GenBank']
        refseq_id = df_genome.loc[idx]['RefSeq']

        #check_genbank = df.loc[idx]['genbank']
        check_genark = df.loc[idx]['GenArk']
        organism = re.sub(r'\s*\([^)]*\)', '', df_genome.loc[idx]['Organism']).replace(' ', '_')

        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')


        ls_download = []
        """
        # RefSeq
        if len(list(set(['ontology', 'repeatmasker']) & set(ls_types))) > 0:

            ls_search = check_refseq.replace(' ', '').split(',')
            refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_id}_{assembly_id}'
            url_md5sum = f'{refseq_dir}/md5checksums.txt'
            refseq_rm = f'{refseq_dir}/{refseq_id}_{assembly_id}_rm.out.gz'

            print(refseq_dir)

            if check_refseq:
                if 'repeatmasker' in ls_types:
                    if 'repeatmasker' in ls_search:
                        if check_url(refseq_rm):
                            out_name = f'{organism}-{assembly_id}-refseq.repeatMasker.out.gz'
                            download_url(refseq_rm, out_name, url_md5sum=url_md5sum)
                            ls_download.append('refseq_rm')

                if 'ontology' in ls_types:
                    if 'ontology' in ls_search:
                        if check_url(url_md5sum):
                            df_md5 = get_md5 (url_md5sum) # Read md5sum information
                            ls_filename = df_md5['Filename'].str.split('/').str[-1].tolist()

                            for file in ls_filename:
                                if '_gene_ontology.gaf.gz' in file:
                                    out_name = f'{organism}-{assembly_id}-refseq.gene_ontology.gaf.gz'
                                    download_url(f'{refseq_dir}/{file}' , out_name, url_md5sum=url_md5sum)
                                    ls_download.append('refseq_ontolgy')
                                    break
        """

        # GenArk
        if len(list(set(['gc', 'sr', 'td', 'wm', 'cb', 'gap', 'rm', 'cpgisland']) & set(ls_types))) > 0:
            if refseq_id in dic_genark_meta:
                genark_id = refseq_id
            else:
                genark_id = genbank_id

            ls_search = check_genark.replace(' ', '').split(',')

            genark_dir = f'{GENARK_URL}/{genark_id[0:3]}/{genark_id[4:7]}/{genark_id[7:10]}/{genark_id[10:13]}/{genark_id}'

            #genark_rmask = f'{genark_dir}/{genark_id}.repeatMasker.out.gz'
            #genark_rmodel = f'{genark_dir}/{genark_id}.repeatModeler.out.gz'
            genark_gc = f'{genark_dir}/bbi/{genark_id}_{assembly_id}.gc5Base.bw'
            genark_sr = f'{genark_dir}/bbi/{genark_id}_{assembly_id}.simpleRepeat.bb'
            genark_td = f'{genark_dir}/bbi/{genark_id}_{assembly_id}.tandemDups.bb'
            genark_wm = f'{genark_dir}/bbi/{genark_id}_{assembly_id}.windowMasker.bb'
            genark_chrsize = f'{genark_dir}/{genark_id}.chrom.sizes.txt'

            if 'gc' in ls_types:
                if 'gc' in ls_search:
                    out_name = f'{organism}-{assembly_id}-genark.gc5Base.bw'
                    download_url(genark_gc, out_name, recursive=recursive)
                    ls_download.append('genark.gc5Base')

            if 'sr' in ls_types:
                if 'sr' in ls_search:
                    out_name = f'{organism}-{assembly_id}-genark.simpleRepeat.bb'
                    download_url(genark_sr, out_name, recursive=recursive)
                    ls_download.append('genark.simpleRepeat')
            if 'td' in ls_types:
                if 'td' in ls_search:
                    out_name = f'{organism}-{assembly_id}-genark.tandemDups.bb'
                    download_url(genark_td, out_name, recursive=recursive)
                    ls_download.append('genark.tandemDups')
            if 'wm' in ls_types:
                if 'wm' in ls_search:
                    out_name = f'{organism}-{assembly_id}-genark.windowMasker.bb'
                    download_url(genark_wm, out_name, recursive=recursive)
                    ls_download.append('genark.windowMasker')

            ls_bb_files = list_http_files(genark_dir + '/bbi')
            ls_gaps = []
            ls_rmsk = []
            ls_cpg_island = []
            for file in ls_bb_files:
                if '.gap.bb' in file or '.allGaps.bb' in file:
                    ls_gaps.append(file)
                elif '.rmsk.' in file:
                    ls_rmsk.append(file)
                elif '.cpgIslandExt' in file:
                    ls_cpg_island.append(file)

            # gap
            if 'gap' in ls_types:
                if 'gap' in ls_search:
                    for in_name in ls_gaps:
                        if '.gap.bb' in in_name:
                            out_name = f'{organism}-{assembly_id}-genark.gap.bb'
                            download_url(f'{genark_dir}/bbi/{in_name}', out_name, recursive=recursive)
                            ls_download.append('genark.gap')
                        elif '.allGaps.bb' in in_name:
                            out_name = f'{organism}-{assembly_id}-genark.allGaps.bb'
                            download_url(f'{genark_dir}/bbi/{in_name}', out_name, recursive=recursive)
                            ls_download.append('genark.allGaps')

            # repeatMasker
            if 'rmsk' in ls_types:
                if 'rmsk' in ls_search:
                    for in_name in ls_rmsk:
                        if '.rmsk.bb' in in_name:
                            out_name = f'{organism}-{assembly_id}-genark.rmsk.bb'
                            download_url(f'{genark_dir}/bbi/{in_name}', out_name, recursive=recursive)
                            ls_download.append('genark.rmsk')
                        else:
                            suffix = in_name.split('.')[-2]
                            out_name = f'{organism}-{assembly_id}-genark.rmsk.{suffix}.bb'
                            download_url(f'{genark_dir}/bbi/{in_name}', out_name, recursive=recursive)
                            ls_download.append(f'genark.rmsk.{suffix}')
            # CpG island
            if 'cpgisland' in ls_types:
                if 'cpgisland' in ls_search:
                    for in_name in ls_cpg_island:
                        suffix = in_name.split('.')[-2]
                        out_name = f'{organism}-{assembly_id}-genark.{suffix}.bb'
                        download_url(f'{genark_dir}/bbi/{in_name}', out_name, recursive=recursive)
                        ls_download.append(f'genark.{suffix}')

            # Chrom size
            out_chrsize_name = f'{organism}-{assembly_id}-genark.chrom.sizes.txt'
            download_url(genark_chrsize, out_chrsize_name, recursive=recursive)

        dic_download[genbank_id] = ls_download
        print('')

    return dic_download

# Change chromosome names of annotation file
def convert_chr_label_annotation (df, dic_download, style, recursive):

    print(f'# Change chromosome names: {style}-style')

    # Download UCSC genome browser applications
    # & Change the file system permissions of files
    path_bw2bdg = f'{DOWNLOAD_FOLDER_NAME}/bigWigToBedGraph'
    path_bdg2bw = f'{DOWNLOAD_FOLDER_NAME}/bedGraphToBigWig'
    path_bb2bed = f'{DOWNLOAD_FOLDER_NAME}/bigBedToBed'
    path_bed2bb = f'{DOWNLOAD_FOLDER_NAME}/bedToBigBed'

    if not os.path.exists(path_bw2bdg) or not os.path.exists(path_bdg2bw) or \
        not os.path.exists(path_bb2bed) or not os.path.exists(path_bed2bb) or \
        recursive:

        print('  Download UCSC genome browser applications')

        # Determine the system's OS and architecture
        system = platform.system()
        arch = platform.machine()
        if system == 'Linux' and arch == 'x86_64':
            base_url = UCSC_KENT_URL + '/linux.x86_64'
        elif system == 'Darwin' and arch == 'x86_64':
            base_url = UCSC_KENT_URL + '/macOSX.x86_64'
        elif system == 'Darwin' and arch == 'arm64':
            base_url = UCSC_KENT_URL + '/macOSX.arm64'
        else:
            raise ValueError(f'Unsupported system: {system} {arch}')

        print(f'  Confirmed system: {system} {arch}')

        # UCSC-Kent applications
        bw2bdg_url = f'{base_url}/bigWigToBedGraph'
        bdg2bw_url = f'{base_url}/bedGraphToBigWig'
        bb2bed_url = f'{base_url}/bigBedToBed'
        bed2bb_url = f'{base_url}/bedToBigBed'

        if not os.path.exists(path_bw2bdg) or recursive:
            download_url(bw2bdg_url, recursive=recursive, warning=False)
            subprocess.run(['chmod', '+x', path_bw2bdg], check=True)
        if not os.path.exists(path_bdg2bw) or recursive:
            download_url(bdg2bw_url, recursive=recursive, warning=False)
            subprocess.run(['chmod', '+x', path_bdg2bw], check=True)
        if not os.path.exists(path_bb2bed) or recursive:
            download_url(bb2bed_url, recursive=recursive, warning=False)
            subprocess.run(['chmod', '+x', path_bb2bed], check=True)
        if not os.path.exists(path_bed2bb) or recursive:
            download_url(bed2bb_url, recursive=recursive, warning=False)
            subprocess.run(['chmod', '+x', path_bed2bb], check=True)
        print('')

    for idx in df.index:
        assembly_id = df.loc[idx]['Assembly name']
        organism = re.sub(r'\s*\([^)]*\)', '', df.loc[idx]['Organism']).replace(' ', '_')
        genbank_id = df.loc[idx]['GenBank']
        refseq_id = df.loc[idx]['RefSeq']

        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')

        # Download assembly report
        print('  Download assembly report')
        report = False
        # RefSeq
        if refseq_id:
            ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}")
            refseq_assembly_id = ls_source[0]
            refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_assembly_id}'
            refseq_md5sum = f'{refseq_dir}/md5checksums.txt'
            refseq_rp = f'{refseq_dir}/{refseq_assembly_id}_assembly_report.txt'
            out_rp_name = f'{organism}-{assembly_id}_assembly_report.txt'
            # Check and download assembly report.

            if check_url(refseq_rp):
                download_genome_url(refseq_rp, out_rp_name, refseq_md5sum)
                report = True

        # GenBank
        if not report:
            ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{genbank_id[0:3]}/{genbank_id[4:7]}/{genbank_id[7:10]}/{genbank_id[10:13]}")
            genbank_assembly_id = ls_source[0]
            genbank_dir = f'{NCBI_FTP_URL}/{genbank_id[0:3]}/{genbank_id[4:7]}/{genbank_id[7:10]}/{genbank_id[10:13]}/{genbank_assembly_id}'
            genbank_md5sum = f'{genbank_dir}/md5checksums.txt'
            genbank_rp = f'{genbank_dir}/{genbank_assembly_id}_assembly_report.txt'
            out_rp_name = f'{organism}-{assembly_id}_assembly_report.txt'

            if check_url(genbank_rp):
                # Check and download assembly report.
                download_genome_url(genbank_rp, out_rp_name, genbank_md5sum)

            else:
                print('  !! Assembly report is not found in RefSeq and GenBank')
                print("     Chromosome names can't be changed\n")
                break

        # Make integrative dataframe of chromosome label of databases
        df_report_edit = make_chr_dataframe (organism, assembly_id)

        # Check the number of UCSC 'na'
        ucsc_na_count = df_report_edit['UCSC'].str.contains('na', na=False).sum()
        # UCSC check
        if style == 'ucsc' and ucsc_na_count > 1:
            print(f'  {ucsc_na_count} chromosome(s) have not ucsc name')
            break

        # Convert the DataFrame to dictionaries for faster lookup
        genbank_to_ensembl = df_report_edit.set_index('GenBank')['Ensembl'].to_dict()
        refseq_to_ensembl = df_report_edit.set_index('RefSeq')['Ensembl'].to_dict()
        genbank_to_gencode = df_report_edit.set_index('GenBank')['Gencode'].to_dict()
        refseq_to_gencode = df_report_edit.set_index('RefSeq')['Gencode'].to_dict()
        genbank_to_ucsc = df_report_edit.set_index('GenBank')['UCSC'].to_dict()
        refseq_to_ucsc = df_report_edit.set_index('RefSeq')['UCSC'].to_dict()

        # Check the number of UCSC 'na'
        ucsc_na_count = df_report_edit['UCSC'].str.contains('na', na=False).sum()
        # UCSC check
        if style == 'ucsc' and ucsc_na_count > 1:
            print(f'  {ucsc_na_count} chromosome(s) have not ucsc name')
            break

        # Check downloaded files
        ls_download = dic_download[genbank_id]
        print(f'  Downloaded file(s): {ls_download}')

        dic_out_suffix = {'ensembl' : 'ens-id', 'gencode' : 'gc-id', 'ucsc' : 'ucsc-id'}

        for download in ls_download:
            print(f'  - {download}')

            # Input and output name
            if download == 'genark.gc5Base':
                in_raw_file = f'{organism}-{assembly_id}-genark.gc5Base.bw'
                in_file = f'{organism}-{assembly_id}-genark.gc5Base.bedgraph'
                out_file = f'{organism}-{assembly_id}-genark.gc5Base.{dic_out_suffix[style]}.bedgraph'
                out_file_final = f'{organism}-{assembly_id}-genark.gc5Base.{dic_out_suffix[style]}.bw'
            else:
                in_raw_file = f'{organism}-{assembly_id}-{download}.bb'
                in_file = f'{organism}-{assembly_id}-{download}.bed'
                out_file = f'{organism}-{assembly_id}-{download}.{dic_out_suffix[style]}.bed'
                out_file_final = f'{organism}-{assembly_id}-{download}.{dic_out_suffix[style]}.bb'

            # Remove output file (--reculsive)
            if recursive:
                delete_file(f'{out_subfolder}/{out_file_final}')
                delete_file(f'{out_subfolder}/{out_file}')

            # File check in working directory
            ls_download_files = os.listdir(DOWNLOAD_FOLDER_NAME)
            ls_output_folder_files = os.listdir(out_subfolder)

            if in_raw_file in ls_download_files and out_file_final not in ls_output_folder_files:

                # Change file format (binary to readable format)
                if download == 'genark.gc5Base':
                    # bigWig to bedGraph
                    print('    Convert bigwig to bedgraph')
                    subprocess.run([path_bw2bdg, f'{DOWNLOAD_FOLDER_NAME}/{in_raw_file}', f'{DOWNLOAD_FOLDER_NAME}/{in_file}'], check=True)
                else:
                    # bigBed to bed
                    print('    Convert bigBed to bed')
                    subprocess.run([f'{path_bb2bed}', f'{DOWNLOAD_FOLDER_NAME}/{in_raw_file}', f'{DOWNLOAD_FOLDER_NAME}/{in_file}'], check=True)

                print('    Modify chromosome names')
                with open(f'{DOWNLOAD_FOLDER_NAME}/{in_file}', 'r') as f_in, open(f'{out_subfolder}/{out_file}', 'w') as f_out:
                    for line in f_in:

                        if not line.startswith('#'):
                            ls_tmp = line.strip().split('\t')
                            chr_name = ls_tmp[0]

                            # Change the chromosome name
                            if style == 'ensembl':
                                ls_tmp[0] = genbank_to_ensembl.get(chr_name, refseq_to_ensembl.get(chr_name, chr_name))
                            elif style == 'gencode':
                                ls_tmp[0] = genbank_to_gencode.get(chr_name, refseq_to_gencode.get(chr_name, chr_name))
                            elif style == 'ucsc':
                                ls_tmp[0] = genbank_to_ucsc.get(chr_name, refseq_to_ucsc.get(chr_name, chr_name))

                            line_edit = '\t'.join(ls_tmp) + '\n'
                            f_out.write(line_edit)

                        else:
                            f_out.write(line)

            else:
                print('    The output file already exists')
                print('    !! If the file appears to have any problems, please delete it and retry the process\n')
                continue

            # Chrom size file
            in_genark_chrsizee = f'{organism}-{assembly_id}-genark.chrom.sizes.txt'
            out_genark_chrsize = f'{organism}-{assembly_id}-genark.chrom.sizes.{dic_out_suffix[style]}.tmp.txt'
            out_genark_chrsize_sorted = f'{organism}-{assembly_id}-genark.chrom.sizes.{dic_out_suffix[style]}.txt'

            if in_genark_chrsizee in ls_download_files and out_genark_chrsize not in ls_output_folder_files:

                with open(f'{DOWNLOAD_FOLDER_NAME}/{in_genark_chrsizee}', 'r') as f_in, open(f'{out_subfolder}/{out_genark_chrsize}', 'w') as f_out:
                    for line in f_in:

                        if not line.startswith('#'):
                            ls_tmp = line.strip().split('\t')
                            chr_name = ls_tmp[0]

                            # Change the chromosome name
                            if style == 'ensembl':
                                ls_tmp[0] = genbank_to_ensembl.get(chr_name, refseq_to_ensembl.get(chr_name, chr_name))
                            elif style == 'gencode':
                                ls_tmp[0] = genbank_to_gencode.get(chr_name, refseq_to_gencode.get(chr_name, chr_name))
                            elif style == 'ucsc':
                                ls_tmp[0] = genbank_to_ucsc.get(chr_name, refseq_to_ucsc.get(chr_name, chr_name))

                            line_edit = '\t'.join(ls_tmp) + '\n'
                            f_out.write(line_edit)

                        else:
                            f_out.write(line)

                # Sort chrom size file
                with open(f'{out_subfolder}/{out_genark_chrsize_sorted}', 'w') as f_out:
                    subprocess.run(['sort', '-k1,1', '-k2,2n', f'{out_subfolder}/{out_genark_chrsize}'], stdout=f_out, check=True)
                # Remove temporary chrom size file
                delete_file(f'{out_subfolder}/{out_genark_chrsize}')

            # Change file format (binary to readable format)
            if download == 'genark.gc5Base':
                print(f'    Convert bedgraph to bigwig: {out_file_final}')

                # bigWig to bedGraph
                subprocess.run([f'{path_bdg2bw}', f'{out_subfolder}/{out_file}', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)

                # Remove temporary files
                delete_file(f'{DOWNLOAD_FOLDER_NAME}/{in_file}')
                delete_file(f'{out_subfolder}/{out_file}')

            else:
                print(f'    Convert bed to bigbed: {out_file_final}')

                # Sort bed file
                with open(f'{out_subfolder}/{out_file}_sorted', 'w') as f_out:
                    subprocess.run(['sort', '-k1,1', '-k2,2n', f'{out_subfolder}/{out_file}'], stdout=f_out, check=True)

                # bigBed to bed
                ls_download_factor = download.split('.')
                autodql_dir = script_dir + '/autosql'
                if download == 'genark.simpleRepeat':
                    subprocess.run([path_bed2bb, '-tab', '-type=bed4+12', f'-as={autodql_dir}/simpleRepeat.as', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)
                elif download == 'genark.tandemDups':
                    subprocess.run([path_bed2bb, '-type=bed12+1', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)
                elif download == 'genark.windowMasker':
                    subprocess.run([path_bed2bb, '-type=bed3', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)
                elif download == 'genark.allGaps':
                    subprocess.run([path_bed2bb, '-type=bed3', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)
                elif download == 'genark.gap':
                    subprocess.run([path_bed2bb, '-extraIndex=name', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)

                elif download == 'genark.rmsk':
                    subprocess.run([path_bed2bb, '-tab', '-type=bed9+5', f'-as={autodql_dir}/bigRmskBed.as', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)
                elif len(ls_download_factor) == 3 and ls_download_factor[1] == 'rmsk':
                    subprocess.run([path_bed2bb, '-tab', '-type=bed6+10', f'-as={autodql_dir}/rmskBed6+10.as', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)
                elif 'cpgIslandExt' in download:
                    subprocess.run([path_bed2bb, '-tab', '-type=bed4+6', f'-as={autodql_dir}/cpgIslandExt.as', '-verbose=0', f'{out_subfolder}/{out_file}_sorted', f'{out_subfolder}/{out_genark_chrsize_sorted}', f'{out_subfolder}/{out_file_final}'], check=True)

                # Remove temporary files
                delete_file(f'{DOWNLOAD_FOLDER_NAME}/{in_file}')
                delete_file(f'{out_subfolder}/{out_file}')
                delete_file(f'{out_subfolder}/{out_file}_sorted')

            print('')


## ---------------------------------------------
## gencube sequence
## ---------------------------------------------
# Check full accessibility
def process_row_sequence(idx, row, dic_ensembl_meta):
    result = {
        'index': idx,
        'RefSeq': '',
        'Ensembl': ''
    }

    assembly_id = row['Assembly name']
    genbank_id = row['GenBank']
    refseq_id = row['RefSeq']
    check_ensembl = row['Ensembl']

    # RefSeq
    if refseq_id:
        ls_source = list_ftp_directory(NCBI_FTP_HOST, f"genomes/all/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}")
        refseq_assembly_id = ls_source[0]
        refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_assembly_id}'
        url_md5sum = f'{refseq_dir}/md5checksums.txt'
        refseq_rna = f'{refseq_assembly_id}_rna.fna.gz'
        refseq_rna_geno = f'{refseq_assembly_id}_rna_from_genomic.fna.gz'
        refseq_cds_geno = f'{refseq_assembly_id}_cds_from_genomic.fna.gz'
        refseq_pseudo = f'{refseq_assembly_id}_pseudo_without_product.fna.gz'
        refseq_protein = f'{refseq_assembly_id}_protein.faa.gz'
        refseq_trans_cds = f'{refseq_assembly_id}_translated_cds.faa.gz'

        if check_url(url_md5sum, show_output=False):
            df_md5 = get_md5(url_md5sum)  # Read md5sum information
            ls_filename = df_md5['Filename'].str.split('/').str[-1].tolist()

            if refseq_rna in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'rna')
            if refseq_rna_geno in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'rna_genomic')
            if refseq_cds_geno in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'cds_genomic')
            if refseq_pseudo in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'pseudo')
            if refseq_protein in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'pep')
            if refseq_trans_cds in ls_filename:
                result['RefSeq'] = add_string(result['RefSeq'], 'pep_cds')

    # Ensembl Beta
    if check_ensembl:
        if genbank_id in dic_ensembl_meta:
            ensembl_acc = genbank_id
            organism_ens = dic_ensembl_meta[genbank_id].replace(' ', '_')
        elif refseq_id in dic_ensembl_meta:
            ensembl_acc = refseq_id
            organism_ens = dic_ensembl_meta[refseq_id].replace(' ', '_')

        """ Rapid Release
        ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}'
        ls_source = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)

        source = ''
        ls_excluded = ['rnaseq', 'brake', 'statistics']
        if ls_source:
            for tmp in ls_excluded:
                if tmp in ls_source:
                    ls_source.remove(tmp)

        ls_input = []
        for source in ls_source:
            ls_tmp = [source]
            ensembl_geneset_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}/{source}/geneset'
            geneset = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_geneset_dir)[0]
            ensembl_file_dir = f'{ENSEMBL_RAPID_FTP_URL}/species/{organism_ens}/{ensembl_acc}/{source}/geneset/{geneset}'

            url_md5sum = f'{ensembl_file_dir}/md5sum.txt'
            ensembl_cdna = f'{organism}-{ensembl_acc}-{geneset}-cdna.fa.gz'
            ensembl_cds = f'{organism}-{ensembl_acc}-{geneset}-cds.fa.gz'
            ensembl_pep = f'{organism}-{ensembl_acc}-{geneset}-pep.fa.gz'

            if check_url(url_md5sum):
                df_md5 = get_md5(url_md5sum)  # Read md5sum information
                ls_filename = df_md5['Filename'].str.split('/').str[-1].tolist()

                if ensembl_cdna in ls_filename:
                    ls_tmp.append('cdna')
                if ensembl_cds in ls_filename:
                    ls_tmp.append('cds')
                if ensembl_pep in ls_filename:
                    ls_tmp.append('pep')

                ls_input.append(ls_tmp)
            else:
                print('There is not md5sum file')
        """
        # Check source directories
        ls_source = list(ensembl_meta[organism_ens][ensembl_acc].keys())

        source = ''
        ls_excluded = ['genome']
        if ls_source:
            for tmp in ls_excluded:
                if tmp in ls_source:
                    ls_source.remove(tmp)

        ls_input = []
        for source in ls_source:
            ls_tmp = [source]
            # Check the geneset folder name
            genomebuild = list(ensembl_meta[organism_ens][ensembl_acc][source]['geneset'].keys())[0]

            ls_filename = list(ensembl_meta[organism_ens][ensembl_acc][source]['geneset'][genomebuild].keys())
            ensembl_cdna = 'cdna.fa.gz'
            ensembl_pep = 'pep.fa.gz'

            if ensembl_cdna in ls_filename:
                ls_tmp.append('cdna')
            if ensembl_pep in ls_filename:
                ls_tmp.append('pep')

            ls_input.append(ls_tmp)

        if ls_input:
            str_input = ''
            for i in range(len(ls_input)):
                if i != 0:
                    str_input += ' / '
                for j in range(len(ls_input[i])):
                    if j == 0:
                        str_input += f'{ls_input[i][j]}:'
                    else:
                        str_input += f' {ls_input[i][j]}'
                    if j != 0 and j != len(ls_input[i]) - 1:
                        str_input += ','

            result['Ensembl'] = add_string(result['Ensembl'], str_input)

    return result

def check_access_full_sequence(df, dic_ensembl_meta):
    df_full_annotation = df[['Assembly name', 'UCSC']].copy()
    ls_sequence = ['RefSeq', 'Ensembl']

    for label in ls_sequence:
        df_full_annotation[label] = ''

    print('# Check accessible data in databases')

    results = []
    with ThreadPoolExecutor() as executor:
        futures = []
        for idx in df.index:
            row = df.loc[idx]
            futures.append(executor.submit(process_row_sequence, idx, row, dic_ensembl_meta))

        for future in futures:
            results.append(future.result())

    for result in results:
        idx = result.pop('index')
        for label in ls_sequence:
            df_full_annotation.loc[idx, label] = result[label]

    # Searched result
    print(tabulate(df_full_annotation, headers='keys', tablefmt='grid'))
    print('')

    return df_full_annotation

# Download sequence data
def download_sequence(df, df_genome, dic_ensembl_meta, types, recursive):

    ls_types = types.split(',')

    print('# Download sequence data')
    for idx in df.index:
        # Accession or name
        assembly_id = df_genome.loc[idx]['Assembly name']
        genbank_id = df_genome.loc[idx]['GenBank']
        refseq_id = df_genome.loc[idx]['RefSeq']
        organism = re.sub(r'\s*\([^)]*\)', '', df_genome.loc[idx]['Organism']).replace(' ', '_')

        check_refseq = df.loc[idx]['RefSeq']
        check_ensembl = df.loc[idx]['Ensembl']
        #check_ensembl_repeat = df.loc[idx]['ensembl_repeat']

        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')

        # RefSeq
        if len(list(set(['refseq_rna', 'refseq_rna_genomic', 'refseq_cds_genomic', 'refseq_pseudo', 'refseq_pep', 'refseq_pep_cds']) & set(ls_types))) > 0:

            ls_search = check_refseq.replace(' ', '').split(',')
            refseq_dir = f'{NCBI_FTP_URL}/{refseq_id[0:3]}/{refseq_id[4:7]}/{refseq_id[7:10]}/{refseq_id[10:13]}/{refseq_id}_{assembly_id}'
            url_md5sum = f'{refseq_dir}/md5checksums.txt'
            refseq_rna = f'{refseq_dir}/{refseq_id}_{assembly_id}_rna.fna.gz'
            refseq_rna_gen = f'{refseq_dir}/{refseq_id}_{assembly_id}_rna_from_genomic.fna.gz'
            refseq_cds_gen = f'{refseq_dir}/{refseq_id}_{assembly_id}_cds_from_genomic.fna.gz'
            refseq_pseudo = f'{refseq_dir}/{refseq_id}_{assembly_id}_pseudo_without_product.fna.gz'
            refseq_protein = f'{refseq_dir}/{refseq_id}_{assembly_id}_protein.faa.gz'
            refseq_trans_cds = f'{refseq_dir}/{refseq_id}_{assembly_id}_translated_cds.faa.gz'

            if check_refseq:
                if 'refseq_rna' in ls_types:
                    if 'rna' in ls_search:
                        if check_url(refseq_rna):
                            out_name = f'{organism}-{assembly_id}-refseq.rna.fna.gz'
                            download_url(refseq_rna, out_name, url_md5sum=url_md5sum, recursive=recursive, out_path=out_subfolder)

                if 'refseq_rna_genomic' in ls_types:
                    if 'rna_genomic' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq.rna_from_genomic.fna.gz'
                        download_url(refseq_rna_gen, out_name, url_md5sum=url_md5sum, recursive=recursive, out_path=out_subfolder)

                if 'refseq_cds_genomic' in ls_types:
                    if 'cds_genomic' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq.cds_from_genomic.fna.gz'
                        download_url(refseq_cds_gen, out_name, url_md5sum=url_md5sum, recursive=recursive, out_path=out_subfolder)

                if 'refseq_pseudo' in ls_types:
                    if 'pseudo' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq.pseudo_without_product.fna.gz'
                        download_url(refseq_pseudo, out_name, url_md5sum=url_md5sum, recursive=recursive, out_path=out_subfolder)

                if 'refseq_pep' in ls_types:
                    if 'pep' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq.protein.faa.gz'
                        download_url(refseq_protein, out_name, url_md5sum=url_md5sum, recursive=recursive, out_path=out_subfolder)

                if 'refseq_pep_cds' in ls_types:
                    if 'pep_cds' in ls_search:
                        out_name = f'{organism}-{assembly_id}-refseq.translated_cds.faa.gz'
                        download_url(refseq_trans_cds, out_name, url_md5sum=url_md5sum, recursive=recursive, out_path=out_subfolder)

        # Ensembl Beta
        #if len(list(set(['ensembl_cdna', 'ensembl_cds', 'ensembl_pep']) & set(ls_types))) > 0:
        if len(list(set(['ensembl_cdna', 'ensembl_pep']) & set(ls_types))) > 0:

            if genbank_id in dic_ensembl_meta:
                ensembl_acc = genbank_id
                organism_ens = dic_ensembl_meta[genbank_id].replace(' ', '_')
            elif refseq_id in dic_ensembl_meta:
                ensembl_acc = refseq_id
                organism_ens = dic_ensembl_meta[refseq_id].replace(' ', '_')

            if check_ensembl:
                ls_source = []
                for search_in_source in check_ensembl.replace(' ', '').split('/'):
                    source = search_in_source.strip().split(':')[0]
                    ls_source.append(source)
                    ls_search = search_in_source.split(':')[1].split(',')

                """ Rapid Release
                for source in ls_source:
                    ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}/{source}/geneset'
                    # Check the geneset folder nameW
                    geneset = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)[0]

                    ensembl_file_dir = f'{ENSEMBL_RAPID_FTP_URL}/species/{organism_ens}/{ensembl_acc}/{source}/geneset/{geneset}'

                    url_md5sum = f'{ensembl_file_dir}/md5sum.txt'
                    ensembl_cdna = f'{ensembl_file_dir}/{organism}-{ensembl_acc}-{geneset}-cdna.fa.gz'
                    ensembl_cds = f'{ensembl_file_dir}/{organism}-{ensembl_acc}-{geneset}-cds.fa.gz'
                    ensembl_pep = f'{ensembl_file_dir}/{organism}-{ensembl_acc}-{geneset}-pep.fa.gz'

                    if 'ensembl_cdna' in ls_types:
                        if 'cdna' in ls_search:
                            out_name = f'{organism}-{assembly_id}-ensembl_{source}.cdna.fna.gz'
                            download_url(ensembl_cdna, out_name, url_md5sum=url_md5sum, recursive=recursive)

                    if 'ensembl_cds' in ls_types:
                        if 'cds' in ls_search:
                            out_name = f'{organism}-{assembly_id}-ensembl_{source}.cds.fna.gz'
                            download_url(ensembl_cds, out_name, url_md5sum=url_md5sum, recursive=recursive)


                    if 'ensembl_pep' in ls_types:
                        if 'pep' in ls_search:
                            out_name = f'{organism}-{assembly_id}-ensembl_{source}.pep.faa.gz'
                            download_url(ensembl_pep, out_name, url_md5sum=url_md5sum, recursive=recursive)
                """
                for source in ls_source:
                    # Check the geneset folder name
                    genomebuild = list(ensembl_meta[organism_ens][ensembl_acc][source]['geneset'].keys())[0]
                    ensembl_file_dir = f'{ENSEMBL_BETA_FTP_URL}/{organism_ens}/{ensembl_acc}/{source}/geneset/{genomebuild}'
                    # url_md5sum = f'{ensembl_file_dir}/md5sum.txt'
                    ensembl_cdna = f'{ensembl_file_dir}/cdna.fa.gz'
                    ensembl_pep = f'{ensembl_file_dir}/pep.fa.gz'

                    if 'ensembl_cdna' in ls_types:
                        if 'cdna' in ls_search:
                            out_name = f'{organism}-{assembly_id}-ensembl_{source}.cdna.fna.gz'
                            download_url(ensembl_cdna, out_name, recursive=recursive, out_path=out_subfolder)

                    if 'ensembl_pep' in ls_types:
                        if 'pep' in ls_search:
                            out_name = f'{organism}-{assembly_id}-ensembl_{source}.pep.faa.gz'
                            download_url(ensembl_pep, out_name, recursive=recursive, out_path=out_subfolder)

        """
        # Ensembl Repeatmodeler
        if check_ensembl_repeat:
            ensembl_repeat = f'{ENSEMBL_RM_FTP_URL}/{organism.lower()}/{genbank_id}.repeatmodeler.fa'

            if 'ensembl_repeat' in ls_types:
                out_name = f'{organism}-{assembly_id}-ensembl.repeatmodeler.fa'
                download_url(ensembl_repeat, out_name, recursive=recursive)
        """

        print('')


## ---------------------------------------------
## gencube crossgenome
## ---------------------------------------------
# Check full accessibility
def process_row_crossgenome(idx, row, dic_ensembl_meta, df_zoonomia):
    result = {
        'index': idx,
        'Ensembl': '',
        'Zoonomia': ''
    }

    genbank_id = row['GenBank']
    refseq_id = row['RefSeq']
    check_ensembl = row['Ensembl']
    check_zoonomia = row['Zoonomia']

    # Ensembl Beta
    if check_ensembl:
        if genbank_id in dic_ensembl_meta:
            ensembl_acc = genbank_id
            organism_ens = dic_ensembl_meta[genbank_id].replace(' ', '_')
        elif refseq_id in dic_ensembl_meta:
            ensembl_acc = refseq_id
            organism_ens = dic_ensembl_meta[refseq_id].replace(' ', '_')

        """ Rapid Release
        ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}'
        ls_source = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)

        source = ''
        ls_excluded = ['rnaseq', 'brake', 'statistics']
        if ls_source:
            for tmp in ls_source:
                if tmp not in ls_excluded:
                    ensembl_homology_dir = f'{ENSEMBL_RAPID_FTP_URL}/species/{organism_ens}/{ensembl_acc}/{tmp}/homology'
                    if check_url(ensembl_homology_dir, show_output=False):
                        ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}/{tmp}/homology'
                        geneset = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)[0]

                        homology_file = f'{ensembl_homology_dir}/{geneset}/{organism_ens}-{ensembl_acc}-{geneset}-homology.tsv.gz'
                        if check_url(homology_file, show_output=False):
                            if not source:
                                source = tmp
                            else:
                                source += f', {tmp}'
                        else:
                            continue
                    else:
                        continue
        """
        # Check source directories
        ls_source = list(ensembl_meta[organism_ens][ensembl_acc].keys())

        ls_excluded = ['genome']
        if ls_source:
            for tmp in ls_excluded:
                if tmp in ls_source:
                    ls_source.remove(tmp)

        for source in ls_source:
            # Check the geneset folder name
            if 'homology' in list(ensembl_meta[organism_ens][ensembl_acc][source]):
                genomebuild = list(ensembl_meta[organism_ens][ensembl_acc][source]['homology'].keys())[0]

                ls_filename = list(ensembl_meta[organism_ens][ensembl_acc][source]['homology'][genomebuild].keys())
                ensembl_homology = f'{organism_ens}-{ensembl_acc}-{genomebuild}-homology.tsv.gz'

                if ensembl_homology in ls_filename:
                    result['Ensembl'] = add_string(result['Ensembl'], source)

    # Zoonomia
    if check_zoonomia:
        ls_reference = df_zoonomia[
            (df_zoonomia["NCBI accession / source"] == genbank_id) |
            (df_zoonomia["NCBI accession / source"] == refseq_id)
        ]['reference'].tolist()
        reference = ', '.join(list(set(ls_reference)))

        result['Zoonomia'] = add_string(result['Zoonomia'], reference)

    return result

def check_access_full_crossgenome(df, dic_ensembl_meta, df_zoonomia):
    df_full_annotation = df[['Assembly name', 'UCSC']].copy()
    ls_geneset = ['Ensembl', 'Zoonomia']

    for label in ls_geneset:
        df_full_annotation[label] = ''

    print('# Check accessible data in databases')

    results = []
    with ThreadPoolExecutor() as executor:
        futures = []
        for idx in df.index:
            row = df.loc[idx]
            futures.append(executor.submit(process_row_crossgenome, idx, row, dic_ensembl_meta, df_zoonomia))

        for future in futures:
            results.append(future.result())

    for result in results:
        idx = result.pop('index')
        for label in ls_geneset:
            df_full_annotation.loc[idx, label] = result[label]

    # Searched result
    print(tabulate(df_full_annotation, headers='keys', tablefmt='grid'))
    print('')

    return df_full_annotation

# Download crossgenome data
def download_crossgenome (df, df_genome, dic_ensembl_meta, df_zoonomia, types, recursive):
    ls_types = types.split(',')

    print('# Download geneset data')
    for idx in df.index:
        # Accession or name
        assembly_id = df_genome.loc[idx]['Assembly name']
        genbank_id = df_genome.loc[idx]['GenBank']
        refseq_id = df_genome.loc[idx]['RefSeq']
        organism = re.sub(r'\s*\([^)]*\)', '', df_genome.loc[idx]['Organism']).replace(' ', '_')

        check_ensembl = df.loc[idx]['Ensembl']
        check_zoonomia = df.loc[idx]['Zoonomia']

        # Print assembly
        if refseq_id:
            print(f'[{genbank_id} / {refseq_id} / {assembly_id}]')
        else:
            print(f'[{genbank_id} / {assembly_id}]')

        # Ensembl Beta
        if 'ensembl_homology' in ls_types:

            if genbank_id in dic_ensembl_meta:
                ensembl_acc = genbank_id
                organism_ens = dic_ensembl_meta[genbank_id].replace(' ', '_')
            elif refseq_id in dic_ensembl_meta:
                ensembl_acc = refseq_id
                organism_ens = dic_ensembl_meta[refseq_id].replace(' ', '_')


            if check_ensembl:
                ls_source = check_ensembl.replace(' ', '').split(',')

                """ Rapid Release
                for source in ls_source:
                    ensembl_dir = f'/pub/rapid-release/species/{organism_ens}/{ensembl_acc}/{source}/homology'
                    if ensembl_dir:
                        # Check the geneset folder name
                        geneset = list_ftp_directory(ENSEMBL_FTP_HOST, ensembl_dir)[0]

                        ensembl_file_dir = f'{ENSEMBL_RAPID_FTP_URL}/species/{organism_ens}/{ensembl_acc}/{source}/homology/{geneset}'
                        ensembl_homology = f'{ensembl_file_dir}/{organism_ens}-{ensembl_acc}-{geneset}-homology.tsv.gz'

                        if 'ensembl_homology' in ls_types:
                            if check_url(ensembl_homology):
                                out_name = f'{organism}-{assembly_id}-ensembl_{source}_homology.tsv.gz'
                                download_url(ensembl_homology, out_name, recursive=recursive)
                """
                for source in ls_source:
                    # Check the geneset folder name
                    genomebuild = list(ensembl_meta[organism_ens][ensembl_acc][source]['homology'].keys())[0]
                    ensembl_file_dir = f'{ENSEMBL_BETA_FTP_URL}/{organism_ens}/{ensembl_acc}/{source}/homology/{genomebuild}'
                    url_md5sum = f'{ensembl_file_dir}/md5sum.txt'
                    ensembl_homology = f'{ensembl_file_dir}/{organism_ens}-{ensembl_acc}-{genomebuild}-homology.tsv.gz'

                    out_name = f'{organism}-{assembly_id}-ensembl_{source}_homology.tsv.gz'
                    download_url(ensembl_homology, out_name, url_md5sum=url_md5sum, recursive=recursive, out_path=out_subfolder)

        # Zoonomia
        if len(list(set(['toga_homology', 'toga_align_codon', 'toga_align_protein', 'toga_inact_mut']) & set(ls_types))) > 0:

            if check_zoonomia:
                ls_reference = check_zoonomia.replace(' ', '').split(',')

                for reference in ls_reference:

                    zoonomia_dir = f'{ZOONOMIA_URL}/{DIC_ZOONOMIA[reference]}'

                    df_tmp = df_zoonomia[
                        ((df_zoonomia['NCBI accession / source'] == genbank_id) |
                        (df_zoonomia['NCBI accession / source'] == refseq_id)) &
                        (df_zoonomia['reference'] == reference)
                    ]

                    for i in range(len(df_tmp.index)):
                        taxo = df_tmp['Taxonomic Lineage'].values[i]
                        species = df_tmp['Species'].values[i].replace(' ', '_')
                        name = df_tmp['Common name'].values[i].replace(' ', '_')
                        assembly = df_tmp['Assembly name'].values[i]

                        if reference in ['human', 'mouse', 'chicken']:
                            ls_folders = list_http_folders(zoonomia_dir)

                            for folder in ls_folders:
                                if folder in taxo:
                                    category = folder
                                    break

                            zoonomia_file_dir = f'{zoonomia_dir}/{category}/{species}__{name}__{assembly}'
                        else:
                            zoonomia_file_dir = f'{zoonomia_dir}/{species}__{name}__{assembly}'

                        zoonomia_orth = f'{zoonomia_file_dir}/orthologsClassification.tsv.gz'
                        zoonomia_align_codon = f'{zoonomia_file_dir}/codonAlignments.fa.gz'
                        zoonomia_align_codoncesar = f'{zoonomia_file_dir}/codonAlignments.allCESARexons.fa.gz'
                        zoonomia_align_protein = f'{zoonomia_file_dir}/proteinAlignments.fa.gz'
                        zoonomia_align_proteincesar = f'{zoonomia_file_dir}/proteinAlignments.allCESARexons.fa.gz'
                        zoonomia_inact_mut = f'{zoonomia_file_dir}/loss_summ_data.tsv.gz'

                        if 'toga_homology' in ls_types:
                            if check_url(zoonomia_orth, verify=False):
                                out_gtf_name = f'{organism}-{assembly_id}-toga_{reference}_homology.tsv.gz'
                                download_url(zoonomia_orth, out_gtf_name, verify=False, recursive=recursive, out_path=out_subfolder)
                        if 'toga_align_codon' in ls_types:
                            if check_url(zoonomia_align_codon, verify=False):
                                out_gtf_name = f'{organism}-{assembly_id}-toga_{reference}_align_codon.fa.gz'
                                download_url(zoonomia_align_codon, out_gtf_name, verify=False, recursive=recursive, out_path=out_subfolder)
                            if check_url(zoonomia_align_codoncesar, verify=False):
                                out_gtf_name = f'{organism}-{assembly_id}-toga_{reference}_align_codon_cesar.fa.gz'
                                download_url(zoonomia_align_codoncesar, out_gtf_name, verify=False, recursive=recursive, out_path=out_subfolder)
                        if 'toga_align_protein' in ls_types:
                            if check_url(zoonomia_align_protein, verify=False):
                                out_gtf_name = f'{organism}-{assembly_id}-toga_{reference}_align_protein.fa.gz'
                                download_url(zoonomia_align_protein, out_gtf_name, verify=False, recursive=recursive, out_path=out_subfolder)
                            if check_url(zoonomia_align_proteincesar, verify=False):
                                out_gtf_name = f'{organism}-{assembly_id}-toga_{reference}_align_protein_cesar.fa.gz'
                                download_url(zoonomia_align_proteincesar, out_gtf_name, verify=False, recursive=recursive, out_path=out_subfolder)
                        if 'toga_inact_mut' in ls_types:
                            if check_url(zoonomia_inact_mut, verify=False):
                                out_gtf_name = f'{organism}-{assembly_id}-toga_{reference}_loss_summ_data.tsv.gz'
                                download_url(zoonomia_inact_mut, out_gtf_name, verify=False, recursive=recursive, out_path=out_subfolder)

        print('')


## ---------------------------------------------
## gencube seqmeta
## ---------------------------------------------
# Make query for searching
def make_query(
    organism, strategy, source, platform, selection,
    filter, layout, access, bioproject, biosample, accession,
    title, author, publication, modification,
    properties, readlength, mbases, textword, keywords, exclude):

    print('# Make query for searching')
    query = ''
    dic_query_pars = {}
    dic_query_kwds = {}
    dic_query_excs = {}
    # Organism
    if organism:
        query = add_query (organism, 'Organism', query, dic_query_pars)
    # Strategy
    if strategy:
        query = add_query (strategy, 'Strategy', query, dic_query_pars, phrase=True)
    # Source
    if source:
        query = add_query (source, 'Source', query, dic_query_pars, phrase=True)
    # Platform
    if platform:
        query = add_query (platform, 'Platform', query, dic_query_pars, phrase=True)
    # Selection
    if selection:
        query = add_query (selection, 'Selection', query, dic_query_pars, phrase=True)
    # Filter
    if filter:
        query = add_query (filter, 'Filter', query, dic_query_pars, phrase=True)
    # Properties
    if properties:
        query = add_query (properties, 'Properties', query, dic_query_pars, phrase=True)
    # Layout
    if layout:
        query = add_query (layout, 'Layout', query, dic_query_pars, phrase=True)
    # Access
    if access:
        query = add_query (access, 'Access', query, dic_query_pars, phrase=True)
    # BioProject
    if bioproject:
        query = add_query (bioproject, 'BioProject', query, dic_query_pars)
    # BioSample
    if biosample:
        query = add_query (biosample, 'BioSample', query, dic_query_pars)
    # Accession
    if accession:
        query = add_query (accession, 'Accession', query, dic_query_pars)
    # Title
    if title:
        query = add_query (title, 'Title', query, dic_query_pars)
    # Aligned
    #if aligned:
    #    query = add_query (aligned, 'Aligned', query, dic_query_pars)
    # Author
    if author:
        query = add_query (author, 'Author', query, dic_query_pars)

    # Publication Date
    if publication:
        query = add_query_with_date (publication, 'Publication Date', query, dic_query_pars)
    # Modification Date
    if modification:
        query = add_query_with_date (modification, 'Modification Date', query, dic_query_pars)

    # ReadLength
    if readlength:
        query = add_query (readlength, 'ReadLength', query, dic_query_pars)
    # Mbases
    if mbases:
        query = add_query (mbases, 'Mbases', query, dic_query_pars)
    # Text Word
    if textword:
        query = add_query (textword, 'Text Word', query, dic_query_pars)

    # For search in steps
    if len(dic_query_pars.keys()) > 1:
        dic_query_pars['Intersection'] = query

    # Keywords
    tmp_show = ''
    if keywords:
        query_tmp_merge = ''
        count = 0
        print(f'  Keywords: {keywords}')
        for i in range(len(keywords)):

            # OR operator
            ls_keywords = keywords[i].split(',')

            count += len(ls_keywords)
            for j in range(len(ls_keywords)):

                kwd_replace = ls_keywords[j].replace('_', ' ')
                # For phrase search
                if kwd_replace[-1] == '^':
                    kwd_replace = f'"{kwd_replace[:-1]}"'

                if j == 0:
                    query_tmp = f'{kwd_replace}'
                else:
                    query_tmp += f' OR {kwd_replace}'

                # For search in steps
                if dic_query_pars:
                    query_step_tmp = f'{dic_query_pars["Intersection"]} AND {kwd_replace}'
                else:
                    query_step_tmp = kwd_replace
                dic_query_kwds[ls_keywords[j]] = query_step_tmp

            # For search in steps
            if ls_keywords:
                if dic_query_pars:
                    query_step_tmp = f'{dic_query_pars["Intersection"]} AND ({query_tmp})'
                else:
                    query_step_tmp = f'({query_tmp})'
                tmp_key = keywords[i].replace('_', ' ').replace(',', '|')
                dic_query_kwds[tmp_key] = query_step_tmp
                if not tmp_show:
                    tmp_show = tmp_key
                else:
                    tmp_show += f' & {tmp_key}'
            if len(keywords) > 1:
                dic_query_kwds[f'space{i}'] = ' '

            # AND operator
            if not query_tmp_merge:
                query_tmp_merge = f'({query_tmp})'
            else:
                query_tmp_merge += f' AND ({query_tmp})'

        if query:
            query = f'({query} AND ({query_tmp_merge}))'
        else:
            query = query_tmp_merge

        # For search in steps
        if len(keywords) > 1:
            dic_query_kwds[f'Intersection ({tmp_show})'] = query

    # Excluded keywords (NOT operator)
    if exclude:
        ls_exclude = exclude.split(',')
        print(f'  Excluded: {ls_exclude}')

        for i in range(len(ls_exclude)):
            exclude_replace = ls_exclude[i].replace('_', ' ')
            # For phrase search
            if exclude_replace[-1] == '^':
                exclude_replace = f'"{exclude_replace[:-1]}"'

            query_tmp = f'{exclude_replace}'
            query = f'{query} NOT {query_tmp}'

            if i == 0:
                query_tmp_or = f'{exclude_replace}'
            else:
                query_tmp_or += f' OR {exclude_replace}'

            # For search in steps
            if dic_query_pars:
                query_step_tmp = f'{dic_query_pars["Intersection"]} AND {query_tmp}'
            else:
                query_step_tmp = {query_tmp}
            dic_query_excs[exclude_replace] = query_step_tmp
        # For search in steps
        if dic_query_pars:
            query_step_tmp = f'{dic_query_pars["Intersection"]} AND ({query_tmp_or})'
        else:
            query_step_tmp = f'{query_tmp_or}'
        dic_query_excs[exclude.replace('_', ' ').replace(',', '|')] = query_step_tmp

    print('')

    return query, dic_query_pars, dic_query_kwds, dic_query_excs

# Make query for searching - Add query
def add_query (input, option, query, dic, phrase=False):
    dic[option] = ' '
    list = input.split(',')
    print(f'  {option}: {list}')

    for i in range(len(list)):
        tmp = list[i].replace('_', ' ')
        # For phrase search
        if phrase:
            tmp = f'"{tmp}"'
        elif tmp[-1] == '^':
            tmp = f'"{tmp[:-1]}"'

        if i == 0:
            query_tmp = f'{tmp}[{option}]'
        else:
            query_tmp += f' OR {tmp}[{option}]'

        dic[tmp] = f'{tmp}[{option}]'
    if not query:
        return f'({query_tmp})'
    else:
        return f'({query} AND ({query_tmp}))'

# Make query for searching - Add query with date
def add_query_with_date (input, option, query, dic):
    dic[option] = ' '
    list = input.split(':')
    print(f'  {option}: {list}')

    for i in range(len(list)):
        tmp = list[i].replace('.', '/')

        if i == 0:
            if len(list) == 2:
                query_tmp = f'"{tmp}"[{option}]'
            elif len(list) == 1:
                query_tmp = f'"{tmp}"[{option}] : "3000"[{option}]'
        else:
            query_tmp += f' : "{tmp}"[{option}]'

    dic[option] = query_tmp
    return f'({query} AND ({query_tmp}))'

# ESearch
def search_sra(query):
    handle = Entrez.esearch(db="sra",  # Database to search
                            term=query,  # Search term
                            retmax=10000  # Number of results to return
                            )
    record = Entrez.read(handle)

    return record

def fetch_single_meta(id, api_key=None, rate_limit=1):
    time.sleep(rate_limit)  # Rate limiting
    handle = Entrez.efetch(db="sra",
                           id=id, rettype="gb",
                           retmode="text",
                           api_key=api_key)
    xml_data = handle.read()
    dict_data = xmltodict.parse(xml_data)
    json_data = json.dumps(dict_data)
    df_tmp = pd.json_normalize(json.loads(json_data))
    df_tmp.columns = df_tmp.columns.str.replace('EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.', '', regex=True)
    df_tmp.index = [id]
    return df_tmp

def fetch_meta(ls_id):
    if len(ls_id) == 0:
        print("No searched id")
        return pd.DataFrame()

    print('# Fetch metadata')
    start_time = time.time()  # record start time

    df = pd.DataFrame()

    # Determine thread number and rate limit based on the presence of an API key
    num_cores  = os.cpu_count()
    if api_key:
        if num_cores >= 10:
            thread_num = 10
        else:
            thread_num = num_cores
        print(f'  Threads: {thread_num} (NCBI API key applied - 10 requests/sec)')
    else:
        if num_cores >= 3:
            thread_num = 3
        else:
            thread_num = num_cores
        print(f'  Threads: {thread_num} (NCBI API key not applied - 3 requests/sec)')
    #rate_limit = 1/(thread_num)
    rate_limit = 0.3

    # Wrapper function to pass rate_limit to fetch_single_meta
    def fetch_single_meta_with_limit(id):
        if api_key:
            return fetch_single_meta(id, api_key=api_key, rate_limit=rate_limit)
        else:
            return fetch_single_meta(id, rate_limit=rate_limit)

    with ThreadPoolExecutor(max_workers=thread_num) as executor:
        results = list(tqdm(executor.map(fetch_single_meta_with_limit, ls_id),
                            desc="  Fetching metadata",
                            unit="id",
                            total=len(ls_id)
                            )
                       )

    df = pd.concat(results, axis=0)

    end_time = time.time()  # record end time
    elapsed_time = end_time - start_time
    print(f'  Total fetching time: {int(elapsed_time)} seconds\n')

    return df

# Make output dataframe format
def convert_format(df, query):
    if df.shape[0] != 0:

        # Columns of output dataframe
        ls_out_study_label = LS_SRA_META_STUDY_LABEL
        ls_out_experiment_label = LS_SRA_META_EXPERIMENT_LABEL

        # Check intersected values
        ls_study_ovlp = list(set(df.columns.tolist()) & set(LS_SRA_META_STUDY_KEY))
        ls_experiment_ovlp = list(set(df.columns.tolist()) & set(LS_SRA_META_EXPERIMENT_KEY))

        ls_label_study = []
        ls_label_experiment = []
        for id in ls_study_ovlp:
            label = LS_SRA_META_STUDY_LABEL[LS_SRA_META_STUDY_KEY.index(id)]
            ls_label_study.append(label)
        for id in ls_experiment_ovlp:
            label = LS_SRA_META_EXPERIMENT_LABEL[LS_SRA_META_EXPERIMENT_KEY.index(id)]
            ls_label_experiment.append(label)

        # Extract experiment info for re-formatted experiment table
        df_experiment = df[ls_study_ovlp + ls_experiment_ovlp]
        df_experiment.columns = ls_label_study + ls_label_experiment

        df_experiment.reset_index(drop=True, inplace=True)
        df_experiment_edit = df_experiment.copy() # copy

        # Remove unnecessary info
        if 'Submission' in df_experiment.columns:
            df_experiment_edit['Submission'] = df_experiment_edit['Submission'].replace('GEO', '')
        if 'GSE' in df_experiment.columns:
            df_experiment_edit['GSE'] = df_experiment['GSE'].apply(lambda x: x if isinstance(x, str) and x.startswith('GSE') else '') # remove except of GSE ids
        if 'GSM' in df_experiment.columns:
            df_experiment_edit['GSM'] = df_experiment['GSM'].apply(lambda x: x if isinstance(x, str) and x.startswith('GSM') else '') # remove except of GSM ids
        if 'Published' in df_experiment.columns:
            df_experiment_edit['Published'] = df_experiment['Published'].apply(lambda x: x if pd.isna(x) else x.split(' ')[0]) # Leave only YYYY-MM-DD & remove time

        # Merge columns
        # 1. BioProject
        if 'BioProject_alt1' not in df_experiment_edit.columns:
            df_experiment_edit['BioProject_alt1'] = None
        else:
            ls_label_study.remove('BioProject_alt1')
            ls_out_study_label.remove('BioProject_alt1')
        if 'BioProject_alt2' not in df_experiment_edit.columns:
            df_experiment_edit['BioProject_alt2'] = None
        else:
            ls_label_study.remove('BioProject_alt2')
            ls_out_study_label.remove('BioProject_alt2')
        if 'BioProject_alt3' not in df_experiment_edit.columns:
            df_experiment_edit['BioProject_alt3'] = None
        else:
            ls_label_study.remove('BioProject_alt3')
            ls_out_study_label.remove('BioProject_alt3')
        ls_label_study.append('BioProject')

        ls_out_study_label.insert(2, 'BioProject')
        ls_out_experiment_label.insert(2, 'Study')
        # 2. Country
        if 'Country_alt1' not in df_experiment_edit.columns:
            df_experiment_edit['Country_alt1'] = None
        else:
            ls_label_study.remove('Country_alt1')
            ls_out_study_label.remove('Country_alt1')
        if 'Country_alt2' not in df_experiment_edit.columns:
            df_experiment_edit['Country_alt2'] = None
        else:
            ls_label_study.remove('Country_alt2')
            ls_out_study_label.remove('Country_alt2')
        ls_out_study_label.insert(-2, 'Country')

        # Combine the three columns into a single list, dropping NaNs and ensuring uniqueness
        # Drop the original columns if needed
        # 1. BioProject
        df_experiment_edit['BioProject'] = df_experiment_edit[['BioProject_alt1', 'BioProject_alt2', 'BioProject_alt3']].apply(lambda row: ','.join(map(str, pd.unique(row.dropna()))), axis=1)
        df_experiment_edit.drop(columns=['BioProject_alt1', 'BioProject_alt2', 'BioProject_alt3'], axis=1, inplace=True)
        # 2. Country
        df_experiment_edit['Country'] = df_experiment_edit[['Country_alt1', 'Country_alt2']].apply(lambda row: ','.join(map(str, pd.unique(row.dropna()))), axis=1)
        df_experiment_edit.drop(columns=['Country_alt1', 'Country_alt2'], axis=1, inplace=True)

        ## Table for study ------------------------
        df_study = df_experiment_edit[ls_label_study]
        df_study_dropdup = df_study.drop_duplicates(subset='Study', keep='first') # remove duplicates
        df_study_dropdup.reset_index(drop=True, inplace=True) # remove index ids
        # Count the number of experiments
        df_study_edit = df_study_dropdup.copy()
        df_study_edit['# Experiment'] = [df_experiment_edit['Study'].tolist().count(val) for val in df_study_dropdup['Study']]

        # Make final output dataframe
        df_out_study = pd.DataFrame(columns=ls_out_study_label) # Make empty dataframe for final output
        df_out_experiment = pd.DataFrame(columns=ls_out_experiment_label) # Make empty dataframe for final output
        df_out_study = pd.concat([df_out_study, df_study_edit])

        # Convert the boolean Series to a list of column labels that are in df_experiment_edit
        valid_columns = list(set(df_experiment_edit.columns.tolist()) & set(ls_out_experiment_label))
        df_out_experiment = pd.concat([df_out_experiment, df_experiment_edit[valid_columns]])

        # Expand tree structure information to table format
        # Process 'Sample attribute' column
        df_out_experiment = expand_attributes(df_out_experiment, 'Sample attribute')
        # Process 'Experiment attribute' column
        df_out_experiment = expand_attributes(df_out_experiment, 'Experiment attribute')
        # Process 'File information' column
        df_out_experiment = expand_attributes(df_out_experiment, 'File information', tag_key=None, value_key=None)
        # Remove columns
        df_out_experiment = df_out_experiment.drop(columns=['Sample attribute', 'Experiment attribute', 'File information'])

        # Add query info.
        df_out_study['Query'] = ''
        df_out_experiment['Query'] = ''
        df_out_study.at[0, 'Query'] = query
        df_out_experiment.at[0, 'Query'] = query

        ## Print the number of study and experiment
        print('# Confirmed total numbers')
        print(f'  Study     : {len(df_out_study.index)}')
        print(f'  Experiment: {len(df_out_experiment.index)}\n')

        return df_out_study, df_out_experiment

def expand_attributes(df, column_name, tag_key='TAG', value_key='VALUE'):
    """
    Expands a column with complex data structures (lists/dicts) into multiple columns.
    """
    # Check if the column exists and has any non-null values
    if column_name not in df.columns or df[column_name].isnull().all():
        print(f"Column '{column_name}' does not exist or is entirely empty. Skipping expansion.")
        return df

    # Create a set to hold all unique keys
    all_tags = set()

    # Function to parse each cell
    def parse_cell(cell):
        # Check if cell is scalar
        if pd.api.types.is_scalar(cell):
            if pd.isnull(cell) or (isinstance(cell, str) and cell == ''):
                return []
            if isinstance(cell, str):
                try:
                    # Safely evaluate the string to a Python object
                    return ast.literal_eval(cell)
                except Exception as e:
                    print(f"Error parsing cell:\n{cell}\n{e}")
                    return []
            else:
                return []
        # If cell is array-like (list, tuple, np.ndarray, pd.Series)
        elif isinstance(cell, (list, tuple, np.ndarray, pd.Series)):
            return cell
        else:
            return []

    # Updated flatten_dict function
    def flatten_dict(d, parent_key='', sep='_'):
        items = []
        if isinstance(d, dict):
            for k, v in d.items():
                new_key = f"{parent_key}{sep}{k}" if parent_key else k
                if isinstance(v, dict):
                    items.extend(flatten_dict(v, new_key, sep=sep).items())
                elif isinstance(v, list):
                    for i, item in enumerate(v):
                        if isinstance(item, dict):
                            items.extend(flatten_dict(item, f"{new_key}{sep}{i}", sep=sep).items())
                        else:
                            items.append((f"{new_key}{sep}{i}", item))
                else:
                    items.append((new_key, v))
        elif isinstance(d, list):
            for i, item in enumerate(d):
                new_key = f"{parent_key}{sep}{i}" if parent_key else str(i)
                if isinstance(item, dict):
                    items.extend(flatten_dict(item, new_key, sep=sep).items())
                else:
                    items.append((new_key, item))
        else:
            items.append((parent_key, d))
        return dict(items)

    # Parse the column and collect all unique tags
    parsed_column = df[column_name].apply(parse_cell)
    for items in parsed_column:
        for item in items:
            if isinstance(item, dict):
                if tag_key and value_key and tag_key in item and value_key in item:
                    tag = item[tag_key]
                    all_tags.add(tag)
                else:
                    # Flatten the dictionary to get all keys
                    flat_item = flatten_dict(item)
                    all_tags.update(flat_item.keys())

    # Create new columns for each unique tag
    for tag in all_tags:
        new_col_name = f"{column_name}_{tag}"
        if new_col_name not in df.columns:
            df[new_col_name] = pd.Series(dtype='object')  # Initialize column with object dtype

    # Populate the new columns with the corresponding values
    for idx, items in parsed_column.items():
        if not items:
            continue
        flat_items = {}
        for item in items:
            if isinstance(item, dict):
                if tag_key and value_key and tag_key in item and value_key in item:
                    tag = item[tag_key]
                    value = item[value_key]
                    flat_items[tag] = value
                else:
                    # Flatten nested dictionaries
                    flat_item = flatten_dict(item)
                    flat_items.update(flat_item)
        for key, value in flat_items.items():
            new_col_name = f"{column_name}_{key}"
            if new_col_name in df.columns:
                df.at[idx, new_col_name] = value
            else:
                # In case a new key appears that wasn't in all_tags
                df[new_col_name] = pd.Series(dtype='object')
                df.at[idx, new_col_name] = value

    return df

# Save metadata
def save_seq_metadata (df_study, df_experiment, organism, now):
    # Make output folder
    if not os.path.exists(out_subfolder):
        os.mkdir(out_subfolder)

    # Add scientific name in output name if len(organism.split(',')) == 1
    if len(organism.split(',')) == 1:
        sci_name = organism.replace(' ', '_').lower()
        out_name_study = f'Meta_seq_{sci_name}_{now}_study_n{df_study.shape[0]}.txt'
        out_name_experiment = f'Meta_seq_{sci_name}_{now}_experiment_n{df_experiment.shape[0]}.txt'
        out_name_url = f'Meta_seq_{sci_name}_{now}_experiment_n{df_experiment.shape[0]}_urls.txt'
    else:
        out_name_study = f'Meta_seq_{now}_study_n{df_study.shape[0]}.txt'
        out_name_experiment = f'Meta_seq_{now}_experiment_n{df_experiment.shape[0]}.txt'
        out_name_url = f'Meta_seq_{now}_experiment_n{df_experiment.shape[0]}_urls.txt'

    # Export the outputs
    df_study.to_csv(f'{out_subfolder}/{out_name_study}', sep='\t', header=True, index=False)
    df_experiment.to_csv(f'{out_subfolder}/{out_name_experiment}'
                     , sep='\t', header=True, index=False)

    print('# Metadata are saved')
    print(f'  {out_name_study}')
    print(f'  {out_name_experiment}')
    print('')

    return out_name_url

# Get fastq links from EBI ENA
def get_fastq_dataframe(df_experiment, out_name_url):

    print('# Retrieve the URL address of the raw data')
    """
    Takes a list of experiment accession IDs (SRX, ERX, DRX) and returns a pandas DataFrame
    with columns: 'Accession', 'Type' ('paired' or 'single'), 'URL'.
    Duplicate URLs are removed. Progress is displayed during execution.

    Parameters:
        exp_ids (list): A list of experiment accession IDs (strings).

    Returns:
        pd.DataFrame: A DataFrame with the requested columns, without duplicate URLs.
    """
    exp_ids = df_experiment['Experiment'].tolist()

    results = []
    total_tasks = len(exp_ids)
    completed_tasks = 0

    def fetch_fastq(exp_id):
        ena_url = (
            f'https://www.ebi.ac.uk/ena/portal/api/filereport'
            f'?accession={exp_id}&result=read_run&fields=fastq_ftp&format=json'
        )

        try:
            response = requests.get(ena_url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                for entry in data:
                    if 'fastq_ftp' in entry and entry['fastq_ftp']:
                        ftp_links = entry['fastq_ftp'].split(';')
                        ftp_links = ['ftp://' + url for url in ftp_links]
                        # Determine if paired or single based on number of files
                        file_type = 'paired' if len(ftp_links) == 2 else 'single'
                        for url in ftp_links:
                            results.append({
                                'Experiment': exp_id,
                                'Type': file_type,
                                'URL': url
                            })
                    else:
                        print(f"No FastQ URLs found for {exp_id}")
            else:
                print(f"Failed to fetch data for {exp_id}: HTTP {response.status_code}")
        except Exception as e:
            print(f"Error fetching {exp_id}: {e}")

    # Use ThreadPoolExecutor to fetch data in parallel
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_exp = {executor.submit(fetch_fastq, exp_id): exp_id for exp_id in exp_ids}
        for future in as_completed(future_to_exp):
            exp_id = future_to_exp[future]
            try:
                future.result()
            except Exception as e:
                print(f"Error processing {exp_id}: {e}")
            # Update progress
            completed_tasks += 1
            progress = (completed_tasks / total_tasks) * 100
            print(f"  Progress: {completed_tasks}/{total_tasks} ({progress:.2f}%) completed", end='\r')

    print()  # Move to the next line after progress is complete

    # Convert results to DataFrame
    df = pd.DataFrame(results, columns=['Experiment', 'Type', 'URL'])

    # Remove duplicate URLs
    df = df.drop_duplicates(subset='URL', keep='first').reset_index(drop=True)

    df.to_csv(f'{out_subfolder}/{out_name_url}'
              , sep='\t', header=True, index=False)
    print(f'  {out_name_url}')
    print('')

# Join variables with newlines
def join_variables_with_newlines(items, max_line_length=100):
    lines = []
    current_line = ""

    for item in items:
        if current_line:  # If current_line is not empty
            if len(current_line) + len(item) + 2 <= max_line_length:  # 2 for ', '
                current_line += ", " + item
            else:
                lines.append(current_line)
                current_line = item
        else:
            current_line = item

    if current_line:  # Append any remaining items in current_line to lines
        lines.append(current_line)

    return '\n'.join(lines)
