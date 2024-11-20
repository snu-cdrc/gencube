import pytest
import pandas as pd
from gencube.gencube_genome import genome
from unittest.mock import patch
from tabulate import tabulate

@pytest.fixture
def mock_utils(mocker):
    # Mock the utility functions used in genome
    mocker.patch('gencube.utils.check_now', return_value='2024-05-27')
    mocker.patch('gencube.utils.search_assembly', return_value=[{
        'assembly_accession': 'GCF_000001405.39',
        'refseq_category': 'representative genome',
        'taxid': 9606,
        'species_taxid': 9606,
        'organism_name': 'Homo sapiens',
        'infraspecific_name': 'GRCh38.p13',
        'isolate': '',
        'version_status': 'latest',
        'assembly_level': 'Chromosome',
        'release_type': 'major',
        'genome_representation': 'full',
        'seq_rel_date': '2024-05-01',
        'asm_name': 'GRCh38.p13',
        'submitter': 'Genome Reference Consortium',
        'gbrs_paired_asm': 'GCA_000001405.28',
        'paired_asm_comp': 'identical'
    }])
    mocker.patch('gencube.utils.json_to_dataframe', return_value=pd.DataFrame({
        'RefSeq': ['GCF_000001405.39'],
        'Genbank': ['GCA_000001405.28'],
        'Assembly name': ['GRCh38.p13'],
        'UCSC': ['hg38'],
        'Release': ['2024-05-01'],
        'Taxid': [9606],
        'Organism': ['Homo sapiens'],
        'Assembly type': ['Primary Assembly'],
        'Level': ['Chromosome'],
        'ContigN50': [50000000],
        'ScaffoldN50': [50000000],
        'Coverage': ['50x'],
        'ReleaseLevel': ['Full'],
        'ReleaseType': ['Major'],
        'AsmReleaseDate_GenBank': ['2024-05-01'],
        'AsmUpdateDate': ['2024-05-15'],
        'SubmitterOrganization': ['Genome Reference Consortium'],
        'Biosource': ['blood'],
        'PropertyList': ['Human Genome']
    }))
    mocker.patch('gencube.utils.check_access_database', return_value=(pd.DataFrame({
        'RefSeq': ['GCF_000001405.39'],
        'Genbank': ['GCA_000001405.28'],
        'Assembly name': ['GRCh38.p13'],
        'UCSC': ['hg38'],
        'Release': ['2024-05-01'],
        'Taxid': [9606],
        'Organism': ['Homo sapiens'],
        'Assembly type': ['Primary Assembly'],
        'Level': ['Chromosome'],
        'ContigN50': [50000000],
        'ScaffoldN50': [50000000],
        'Coverage': ['50x'],
        'ReleaseLevel': ['Full'],
        'ReleaseType': ['Major'],
        'AsmReleaseDate_GenBank': ['2024-05-01'],
        'AsmUpdateDate': ['2024-05-15'],
        'SubmitterOrganization': ['Genome Reference Consortium'],
        'Biosource': ['blood'],
        'PropertyList': ['Human Genome']
    }), {}, {}))
    mocker.patch('gencube.utils.save_genome_metadata')
    mocker.patch('gencube.utils.download_genome', return_value={'GCA_000001405.28': ['refseq']})
    mocker.patch('gencube.utils.convert_chr_label_genome')

def test_genome(mock_utils, capsys):
    keywords = ['Homo_sapiens']
    level = 'Chromosome'
    refseq = True
    ucsc = True
    latest = True
    metadata = True
    download = True
    fasta = 'refseq'
    chr_style = 'ensembl'
    compresslevel = 9

    genome(keywords, level, refseq, ucsc, latest, metadata, download, fasta, chr_style, compresslevel)
    
    # Capture the printed output
    captured = capsys.readouterr()
    expected_output = tabulate(pd.DataFrame({
        'RefSeq': ['GCF_000001405.39'],
        'Genbank': ['GCA_000001405.28'],
        'Assembly name': ['GRCh38.p13'],
        'UCSC': ['hg38'],
        'Release': ['2024-05-01'],
        'Taxid': [9606],
        'Organism': ['Homo sapiens'],
        'Assembly type': ['Primary Assembly'],
        'Level': ['Chromosome'],
        'ContigN50': [50000000],
        'ScaffoldN50': [50000000],
        'Coverage': ['50x'],
        'ReleaseLevel': ['Full'],
        'ReleaseType': ['Major'],
        'AsmReleaseDate_GenBank': ['2024-05-01'],
        'AsmUpdateDate': ['2024-05-15'],
        'SubmitterOrganization': ['Genome Reference Consortium'],
        'Biosource': ['blood'],
        'PropertyList': ['Human Genome']
    })[['Assembly name', 'Taxid', 'Release', 'RefSeq', 'Genbank', 'UCSC', 'Coverage', 'SubmitterOrganization', 'Biosource', 'PropertyList']], headers='keys', tablefmt='grid') + '\n\n'
    assert captured.out == expected_output

if __name__ == "__main__":
    pytest.main()
