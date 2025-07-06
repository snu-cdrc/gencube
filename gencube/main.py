import argparse
import sys
import requests

# Sub-modules
from .__init__ import __version__
from .gencube_genome import genome
from .gencube_geneset import geneset
from .gencube_annotation import annotation
from .gencube_sequence import sequence
from .gencube_crossgenome import crossgenome
from .gencube_seqmeta import seqmeta
from .gencube_info import info

# Custom functions
from .utils import (
    get_entrez_info,
    join_variables_with_newlines,
    )
# Constant variables
from .constants import (
    SRA_SEARCH_INDEX_JSON,
    )

def main():
    # Get email and api key information for E-utility
    info_save = get_entrez_info()

    # Define parent parser
    parser = argparse.ArgumentParser(
        prog='gencube', description=f"gencube v{__version__}"
        )
    # Initiate subparsers
    subparsers = parser.add_subparsers(
        dest='command'
        )

    # If the user inputs the "seqmeta" subcommand, fetch the search index information
    LS_STRATEGY = ''
    LS_SOURCE = ''
    LS_PLATFORM = ''
    LS_SELECTION = ''
    LS_FILTER = ''
    LS_PROPERTIES = ''
    if len(sys.argv) > 1 and sys.argv[1] == "seqmeta":
        response = requests.get(SRA_SEARCH_INDEX_JSON)
        if response.status_code == 200:
            data = response.json()  # Convert to JSON format
            LS_STRATEGY = data.get("LS_STRATEGY", [])
            LS_SOURCE = data.get("LS_SOURCE", [])
            LS_PLATFORM = data.get("LS_PLATFORM", [])
            LS_SELECTION = data.get("LS_SELECTION", [])
            LS_FILTER = data.get("LS_FILTER", [])
            LS_PROPERTIES = data.get("LS_PROPERTIES", [])

    ## ---------------------------------------------
    ## gencube genome
    ## ---------------------------------------------
    # gencube genome subparser
    genome_desc = 'Search, download, and modify chromosome labels for genome assemblies'
    parser_genome = subparsers.add_parser(
        'genome', description=genome_desc, help=genome_desc, add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
        )
    # gencube genome arguments
    parser_genome.add_argument(
        'keywords',
        nargs='*',
        help=(
            'Taxonomic names to search for genomes\n'
            'You can provide various forms such as species names or accession numbers\n'
            "Examples: homo_sapiens, human, GCF_000001405.40, GCA_000001405.29, GRCh38, hg38\n"
            '\n'
            'Multiple names can be combined and will be merged in the search results\n'
            'To specify multiple names, separate them with spaces'
            )
        )
    parser_genome.add_argument(
        '-v', '--level',
        metavar='level',
        default='complete,chromosome',
        required=False,
        help=(
            'Specify the genome assembly level (default: complete,chromosome)\n'
            'complete   : Fully assembled genomes\n'
            'chromosome : Assembled at the chromosome level\n'
            'scaffold   : Assembled into scaffolds, but not to the chromosome level\n'
            'contig     : Contiguous sequences without gaps\n'
            '\n'
            )
        )
    parser_genome.add_argument(
        '-r', '--refseq',
        action='store_true',
        required=False,
        help='Show genomes that have RefSeq accession (GCF_* format)'
        )
    parser_genome.add_argument(
        '-u', '--ucsc',
        action='store_true',
        required=False,
        help='Show genomes that have UCSC name'
        )
    parser_genome.add_argument(
        '-l', '--latest',
        action='store_true',
        required=False,
        help='Show genomes corresponding to the latest version'
        )
    parser_genome.add_argument(
        '-m', '--metadata',
        action='store_true',
        required=False,
        help='Save metadata for the searched genomes'
        )
    parser_genome.add_argument(
        '-d', '--download',
        action='store_true',
        required=False,
        help='Download "fasta" formatted genome file'
        )
    parser_genome.add_argument(
        '-db', '--database',
        metavar='types',
        default='refseq',
        required=False,
        help=(
            'Database where genome file is downloaded (default: refseq)\n'
            'Default is from the RefSeq database\n'
            'If not available, download from the GenBank database\n'
            'genbank : by NCBI GenBank\n'
            'refseq  : by NCBI RefSeq\n'
            'genark  : by UCSC GenArk\n'
            'ensembl : by Ensembl Beta\n'
            )
        )
    parser_genome.add_argument(
        '-c','--chr_style',
        metavar='type',
        default='ensembl',
        choices=['ensembl', 'gencode', 'ucsc', 'raw'],
        required=False,
        help=(
            'Chromosome label style used in the download file (default: ensembl)\n'
            'ensembl : 1, 2, X, MT & unknowns (GenBank IDs)\n'
            'gencode : chr1, chr2, chrX, chrM & unknowns (GenBank IDs)\n'
            'ucsc    : chr1, chr2, chrX, chrM & unknowns (UCSC-specific IDs)\n'
            '          !! Limited use if UCSC IDs are not issued\n'
            'raw     : Uses raw file labels without modification\n'
            '         - NCBI GenBank: CM_* or other-form IDs\n'
            '         - NCBI RefSeq : NC_*, NW_* or other-form IDs\n'
            '         - GenArk      : GenBank or RefSeq IDs\n'
            '         - Ensembl     : Ensembl IDs\n'
            )
        )
    parser_genome.add_argument(
        '-mk','--masking',
        metavar='type',
        default='soft',
        choices=['soft', 'hard', 'none'],
        required=False,
        help=(
            'Masking type for output data (default: soft)\n'
            'soft : soft-masked\n'
            'hard : hard-masked\n'
            'none : unmasked\n'
            )
        )
    parser_genome.add_argument(
        '-cl','--compresslevel',
        metavar='1-9',
        default='6',
        choices=['1', '2', '3', '4', '5', '6', '7', '8', '9'],
        required=False,
        help=(
            'Compression level for output data (default: 6)\n'
            'Lower numbers are faster but have lower compression\n '
            )
        )
    parser_genome.add_argument(
        '--recursive',
        action='store_true',
        required=False,
        help='Download files regardless of their presence only if integrity check is not possible\n '
        )

    ## ---------------------------------------------
    ## gencube geneset
    ## ---------------------------------------------
    # gencube geneset subparser
    geneset_desc = 'Search, download, and modify chromosome labels for genesets (gene annotations)'
    parser_geneset = subparsers.add_parser(
        'geneset', description=geneset_desc, help=geneset_desc, add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
        )
    # gencube geneset arguments
    parser_geneset.add_argument(
        'keywords',
        nargs='*',
        help=(
            'Taxonomic names to search for genomes\n'
            'You can provide various forms such as species names or accession numbers\n'
            "Examples: homo_sapiens, human, GCF_000001405.40, GCA_000001405.29, GRCh38, hg38\n"
            '\n'
            'Multiple names can be combined and will be merged in the search results\n'
            'To specify multiple names, separate them with spaces'
            )
        )
    parser_geneset.add_argument(
        '-v', '--level',
        metavar='level',
        default='complete,chromosome',
        required=False,
        help=(
            'Specify the genome assembly level (default: complete,chromosome)\n'
            'complete   : Fully assembled genomes\n'
            'chromosome : Assembled at the chromosome level\n'
            'scaffold   : Assembled into scaffolds, but not to the chromosome level\n'
            'contig     : Contiguous sequences without gaps\n'
            '\n'
            )
        )
    parser_geneset.add_argument(
        '-r', '--refseq',
        action='store_true',
        required=False,
        help='Show genomes that have RefSeq accession (GCF_* format)'
        )
    parser_geneset.add_argument(
        '-u', '--ucsc',
        action='store_true',
        required=False,
        help='Show genomes that have UCSC name'
        )
    parser_geneset.add_argument(
        '-l', '--latest',
        action='store_true',
        required=False,
        help='Show genomes corresponding to the latest version'
        )
    parser_geneset.add_argument(
        '-m', '--metadata',
        action='store_true',
        required=False,
        help='Save metadata for the searched genesets'
        )
    parser_geneset.add_argument(
        '-d', '--download',
        metavar='types',
        required=False,
        help=(
            'Type of gene set\n'
            'refseq_gtf    : RefSeq gene set (GTF format)\n'
            'refseq_gff    : RefSeq gene set (GFF)\n'
            'gnomon        : RefSeq Gnomon gene prediction (GFF)\n'
            'cross         : RefSeq Cross-species alignments (GFF)\n'
            'same          : RefSeq Same-species alignments (GFF)\n'
            'augustus      : GenArk Augustus gene prediction (GFF)\n'
            'xenoref       : GenArk XenoRefGene (GFF)\n'
            'genark_ref    : GenArk RefSeq gene models (GFF)\n'
            'ensembl_gtf   : Ensembl Beta gene set (GTF)\n'
            'ensembl_gff   : Ensembl Beta gene set (GFF)\n'
            'toga_gtf      : Zoonomia TOGA gene set (GTF)\n'
            'toga_bed      : Zoonomia TOGA gene set (BED)\n'
            'toga_pseudo   : Zoonomia TOGA processed pseudogenes (BED)\n'
            )
        )
    parser_geneset.add_argument(
        '-c','--chr_style',
        metavar='type',
        default='ensembl',
        choices=['ensembl', 'gencode', 'ucsc', 'raw'],
        required=False,
        help=(
            'Chromosome label style used in the download file (default: ensembl)\n'
            'ensembl : 1, 2, X, MT & unknowns (GenBank IDs)\n'
            'gencode : chr1, chr2, chrX, chrM & unknowns (GenBank IDs)\n'
            'ucsc    : chr1, chr2, chrX, chrM & unknowns (UCSC-specific IDs)\n'
            '          !! Limited use if UCSC IDs are not issued\n'
            'raw     : Uses raw file labels without modification\n'
            '         - NCBI GenBank: CM_* or other-form IDs\n'
            '         - NCBI RefSeq : NC_*, NW_* or other-form IDs\n'
            '         - GenArk      : GenBank or RefSeq IDs\n'
            '         - Ensembl     : Ensembl IDs\n '
            )
        )
    parser_geneset.add_argument(
        '--recursive',
        action='store_true',
        required=False,
        help='Download files regardless of their presence only if integrity check is not possible\n '
        )

    ## ---------------------------------------------
    ## gencube annotation
    ## ---------------------------------------------
    # gencube annotation subparser
    annotation_desc = 'Search, download, and modify chromosome labels for various genome annotations, such as gaps and repeats'
    parser_annotation = subparsers.add_parser(
        'annotation', description=annotation_desc, help=annotation_desc, add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
        )
    # gencube annotation arguments
    parser_annotation.add_argument(
        'keywords',
        nargs='*',
        help=(
            'Taxonomic names to search for genomes\n'
            'You can provide various forms such as species names or accession numbers\n'
            "Examples: homo_sapiens, human, GCF_000001405.40, GCA_000001405.29, GRCh38, hg38\n"
            '\n'
            'Multiple names can be combined and will be merged in the search results\n'
            'To specify multiple names, separate them with spaces'
            )
        )
    parser_annotation.add_argument(
        '-v', '--level',
        metavar='level',
        default='complete,chromosome',
        required=False,
        help=(
            'Specify the genome assembly level (default: complete,chromosome)\n'
            'complete   : Fully assembled genomes\n'
            'chromosome : Assembled at the chromosome level\n'
            'scaffold   : Assembled into scaffolds, but not to the chromosome level\n'
            'contig     : Contiguous sequences without gaps\n'
            '\n'
            )
        )
    parser_annotation.add_argument(
        '-r', '--refseq',
        action='store_true',
        required=False,
        help='Show genomes that have RefSeq accession (GCF_* format)'
        )
    parser_annotation.add_argument(
        '-u', '--ucsc',
        action='store_true',
        required=False,
        help='Show genomes that have UCSC name'
        )
    parser_annotation.add_argument(
        '-l', '--latest',
        action='store_true',
        required=False,
        help='Show genomes corresponding to the latest version'
        )
    parser_annotation.add_argument(
        '-m', '--metadata',
        action='store_true',
        required=False,
        help='Save metadata for the searched annotations'
        )
    parser_annotation.add_argument(
        '-d', '--download',
        metavar='types',
        required=False,
        help=(
            'Download annotation file.\n'
            'gap : Genomic gaps - AGP defined (bigBed format)\n'
            'sr   : Simple tandem repeats by TRF (bigBed)\n'
            'td   : Tandem duplications (bigBed)\n'
            'wm   : Genomic intervals masked by WindowMasker + SDust (bigBed)\n'
            'rmsk : Repeated elements annotated by RepeatMasker (bigBed)\n'
            'cpg  : CpG Islands - Islands < 300 bases are light green (bigBed)\n'
            'gc   : GC percent in 5-Base window (bigWig)\n'
            )
        )
    parser_annotation.add_argument(
        '-c','--chr_style',
        metavar='type',
        default='ensembl',
        choices=['ensembl', 'gencode', 'ucsc', 'raw'],
        required=False,
        help=(
            'Chromosome label style used in the download file (default: ensembl)\n'
            'ensembl : 1, 2, X, MT & unknowns (GenBank IDs)\n'
            'gencode : chr1, chr2, chrX, chrM & unknowns (GenBank IDs)\n'
            'ucsc    : chr1, chr2, chrX, chrM & unknowns (UCSC-specific IDs)\n'
            '          !! Limited use if UCSC IDs are not issued\n'
            'raw     : Uses raw file labels without modification\n'
            '         - NCBI GenBank: CM_* or other-form IDs\n'
            '         - NCBI RefSeq : NC_*, NW_* or other-form IDs\n'
            '         - GenArk      : GenBank or RefSeq IDs\n'
            '         - Ensembl     : Ensembl IDs\n '
            )
        )
    parser_annotation.add_argument(
        '--recursive',
        action='store_true',
        required=False,
        help='Download files regardless of their presence only if integrity check is not possible\n '
        )

    ## ---------------------------------------------
    ## gencube sequence
    ## ---------------------------------------------
    # gencube sequence subparser
    sequence_desc = 'Search and download sequence data of genesets'
    parser_sequence = subparsers.add_parser(
        'sequence', description=sequence_desc, help=sequence_desc, add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
        )
    # gencube sequence arguments
    parser_sequence.add_argument(
        'keywords',
        nargs='*',
        help=(
            'Taxonomic names to search for genomes\n'
            'You can provide various forms such as species names or accession numbers\n'
            "Examples: homo_sapiens, human, GCF_000001405.40, GCA_000001405.29, GRCh38, hg38\n"
            '\n'
            'Multiple names can be combined and will be merged in the search results\n'
            'To specify multiple names, separate them with spaces'
            )
        )
    parser_sequence.add_argument(
        '-v', '--level',
        metavar='level',
        default='complete,chromosome',
        required=False,
        help=(
            'Specify the genome assembly level (default: complete,chromosome)\n'
            'complete   : Fully assembled genomes\n'
            'chromosome : Assembled at the chromosome level\n'
            'scaffold   : Assembled into scaffolds, but not to the chromosome level\n'
            'contig     : Contiguous sequences without gaps\n'
            '\n'
            )
        )
    parser_sequence.add_argument(
        '-r', '--refseq',
        action='store_true',
        required=False,
        help='Show genomes that have RefSeq accession (GCF_* format)'
        )
    parser_sequence.add_argument(
        '-u', '--ucsc',
        action='store_true',
        required=False,
        help='Show genomes that have UCSC name'
        )
    parser_sequence.add_argument(
        '-l', '--latest',
        action='store_true',
        required=False,
        help='Show genomes corresponding to the latest version'
        )
    parser_sequence.add_argument(
        '-m', '--metadata',
        action='store_true',
        required=False,
        help='Save metadata for the searched sequence data'
        )
    parser_sequence.add_argument(
        '-d', '--download',
        metavar='types',
        required=False,
        help=(
            'Download "fasta" formatted sequence file\n'
            '1. Nucleotide sequences:\n'
            '   refseq_rna         : Accessioned RNA sequences annotated on the genome assembly\n'
            '   refseq_rna_genomic : RNA features based on the genome sequence\n'
            '   refseq_cds_genomic : CDS features based on the genome sequence\n'
            '   refseq_pseudo      : Pseudogene and other gene regions without transcribed RNA or translated protein products\n'
            '   ensembl_cdna       : Ensembl Beta cDNA sequences of transcripts\n'
            # '   ensembl_cds        : Ensembl Beta coding sequences (CDS)\n'
            # '   ensembl_repeat     : Ensembl repeat modeler sequences\n'
            '2. Protein sequences:\n'
            '   refseq_pep         : Accessioned protein sequences annotated on the genome assembly\n'
            '   refseq_pep_cds     : CDS features translated into protein sequences\n'
            '   ensembl_pep        : Ensembl Beta protein sequences\n '
            )
        )
    parser_sequence.add_argument(
        '--recursive',
        action='store_true',
        required=False,
        help='Download files regardless of their presence only if integrity check is not possible\n '
        )

    ## ---------------------------------------------
    ## gencube crossgenome
    ## ---------------------------------------------
    # gencube crossgenome subparser
    crossgenome_desc = (
        'Search and download comparative genomics data, such as homology, and codon or protein alignments'

    )
    parser_crossgenome = subparsers.add_parser(
        'crossgenome', description=crossgenome_desc, help=crossgenome_desc, add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
        )
    # gencube geneset arguments
    parser_crossgenome.add_argument(
        'keywords',
        nargs='*',
        help=(
            'Taxonomic names to search for genomes\n'
            'You can provide various forms such as species names or accession numbers\n'
            "Examples: homo_sapiens, human, GCF_000001405.40, GCA_000001405.29, GRCh38, hg38\n"
            '\n'
            'Multiple names can be combined and will be merged in the search results\n'
            'To specify multiple names, separate them with spaces'
            )
        )
    parser_crossgenome.add_argument(
        '-v', '--level',
        metavar='level',
        default='complete,chromosome',
        required=False,
        help=(
            'Specify the genome assembly level (default: complete,chromosome)\n'
            'complete   : Fully assembled genomes\n'
            'chromosome : Assembled at the chromosome level\n'
            'scaffold   : Assembled into scaffolds, but not to the chromosome level\n'
            'contig     : Contiguous sequences without gaps\n'
            '\n'
            )
        )
    parser_crossgenome.add_argument(
        '-r', '--refseq',
        action='store_true',
        required=False,
        help='Show genomes that have RefSeq accession (GCF_* format)'
        )
    parser_crossgenome.add_argument(
        '-u', '--ucsc',
        action='store_true',
        required=False,
        help='Show genomes that have UCSC name'
        )
    parser_crossgenome.add_argument(
        '-l', '--latest',
        action='store_true',
        required=False,
        help='Show genomes corresponding to the latest version'
        )
    parser_crossgenome.add_argument(
        '-m', '--metadata',
        action='store_true',
        required=False,
        help='Save metadata for the searched comparative genomics data'
        )
    parser_crossgenome.add_argument(
        '-d', '--download',
        metavar='types',
        required=False,
        help=(
            'ensembl_homology   : Homology data from Ensembl Beta,\n'
            '                     detailing gene orthology relationships across species\n'
            'toga_homology      : Homology data from TOGA, providing predictions of\n'
            '                     orthologous genes based on genome alignments\n'
            'toga_align_codon   : Codon alignment data from TOGA, showing aligned codon\n'
            '                     sequences between reference and query species\n'
            'toga_align_protein : Protein alignment data from TOGA, detailing aligned\n'
            '                     protein sequences between reference and query species\n'
            'toga_inact_mut     : List of inactivating mutations from TOGA, identifying\n'
            '                     mutations that disrupt gene function\n '
            )
        )
    parser_crossgenome.add_argument(
        '--recursive',
        action='store_true',
        required=False,
        help='Download files regardless of their presence only if integrity check is not possible\n '
        )

    ## ---------------------------------------------
    ## gencube seqmeta
    ## ---------------------------------------------
    # gencube seqmeta subparser
    seq_desc = (
        'Search, retrive, and integrate metadata of experimental sequencing data'
        )
    parser_seqmeta = subparsers.add_parser(
        'seqmeta', description=seq_desc, help=seq_desc, add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser_seqmeta.add_argument(
        'keywords',
        nargs='*',
        help=(
            'Keywords to search for sequencing-based experimental data. You can provide various forms\n'
            'Examples: liver, k562, cancer, breast_cancer, etc\n'
            '\n'
            'Multiple keywords can be combined\n'
            'Keywords separated by commas will combine their results\n'
            'Keywords separated by spaces will intersect their results\n'
            'Example: liver,lung cancer,tumor'
            )
        )
    # gencube seqmeta arguments
    parser_seqmeta.add_argument(
        '-o', '--organism',
        metavar='string',
        default='',
        help=(
            'Scientific name or common name (as found in the NCBI Taxonomy Browser)\n' +
            'Example: homo_sapiens or human'
            )
        )
    parser_seqmeta.add_argument(
        '-st', '--strategy',
        metavar='string',
        default='',
        help=(
            'Sequencing strategy:\n' +
            join_variables_with_newlines(LS_STRATEGY)
            )
        )
    parser_seqmeta.add_argument(
        '-sr', '--source',
        metavar='string',
        default='',
        help=(
            'Source of the biological data:\n' +
            join_variables_with_newlines(LS_SOURCE)
            )
        )
    parser_seqmeta.add_argument(
        '-pl', '--platform',
        metavar='string',
        default='',
        help=(
            'Name of the sequencing platform:\n' +
            join_variables_with_newlines(LS_PLATFORM)
            )
        )
    parser_seqmeta.add_argument(
        '-sl', '--selection',
        metavar='string',
        default='',
        help=(
            'Library selection methodology:\n' +
            join_variables_with_newlines(LS_SELECTION)
            )
        )
    parser_seqmeta.add_argument(
        '-fi', '--filter',
        metavar='string',
        default='',
        help=(
            'Option to find SRA records that are cross-referenced with other NCBI databases\n' +
            '(PubMed, PubMed Central (PMC), Nucleotide, Assembly, and others):\n' +
            join_variables_with_newlines(LS_FILTER)
            )
        )
    parser_seqmeta.add_argument(
        '-pr', '--properties',
        metavar='string',
        default='',
        help=(
            "Option to narrow search results by controlled-vocabulary library's annotations:\n" +
            join_variables_with_newlines(LS_PROPERTIES)
            )
        )
    parser_seqmeta.add_argument(
        '-ly', '--layout',
        metavar='string',
        default='',
        help=(
            'Library layout of the sequencing data:\n' +
            'paired, single'
            )
        )
    parser_seqmeta.add_argument(
        '-ac', '--access',
        metavar='string',
        default='',
        help=(
            'Data accessibility:\n' +
            'public, controlled'
            )
        )
    parser_seqmeta.add_argument(
        '-bp', '--bioproject',
        metavar='string',
        default='',
        help=(
            'BioProject accession in the form of PRJNA#, PRJEB#, or PRJDB#'
            )
        )
    parser_seqmeta.add_argument(
        '-bs', '--biosample',
        metavar='string',
        default='',
        help=(
            'BioSample accession in the form of SAMN#, SAMEA#, or SAMD#'
            )
        )
    parser_seqmeta.add_argument(
        '-as', '--accession',
        metavar='string',
        default='',
        help=(
            'SRA/ENA/DDBJ accession\n'
            'Study with accessions in the form of SRP#, ERP#, or DRP#\n'
            'Sample with accessions in the form of SRS#, ERS#, or DRS#\n'
            'Experiment with accessions in the form of SRX#, ERX#, or DRX#\n'
            'Run with accessions in the form of SRR#, ERR#, or DRR#'
            )
        )
    parser_seqmeta.add_argument(
        '-ti', '--title',
        metavar='string',
        default='',
        help=(
            'Descriptive name of the dataset'
            )
        )
    parser_seqmeta.add_argument(
        '-at', '--author',
        metavar='string',
        default='',
        help=(
            'Researcher or group that submitted the data\n '
            'Example: SON_KH'
            )
        )
    parser_seqmeta.add_argument(
        '-pd', '--publication',
        metavar='range',
        default='',
        help=(
            'Publication Date\n'
            'YYYY.MM.DD : YYYY.MM.DD format\n'
            'Example: 2016, 2016.07, 2016.07.01, 2016.07:2023.02'
            )
        )
    parser_seqmeta.add_argument(
        '-md', '--modification',
        metavar='range',
        default='',
        help=(
            'Modification Date\n'
            'YYYY.MM.DD : YYYY.MM.DD format\n'
            'Example: 2016, 2016.07, 2016.07.01, 2016.07:2023.02'
            )
        )
    parser_seqmeta.add_argument(
        '-rl', '--readlength',
        metavar='range',
        default='',
        help=(
            'Length of the sequencing reads'
            'Example: 100 or 100:500'
            )
        )
    parser_seqmeta.add_argument(
        '-mb', '--mbases',
        metavar='string',
        default='',
        help=(
            'Number of mega bases in the SRA Runs'
            )
        )
    parser_seqmeta.add_argument(
        '-tw', '--textword',
        metavar='string',
        default='',
        help=(
            'General search term for finding datasets by specific words in metadata'
            )
        )
    parser_seqmeta.add_argument(
        '-ex', '--exclude',
        metavar='keywords',
        default='',
        help=(
            'Exclude the results for the keywords used in this option\n'
            'Example: cell_line,normal,crispr\n '
            )
        )
    parser_seqmeta.add_argument(
        '-d',
        '--detail',
        action='store_true',
        help='Show the number of searched results for each option and keyword\n '
        )
    parser_seqmeta.add_argument(
        '-m',
        '--metadata',
        action='store_true',
        help='Save integrated metadata\n '
        )
    parser_seqmeta.add_argument(
        '-u',
        '--url',
        action='store_true',
        help='Save the file, including the URL address of the raw data (.fastq).\n '
        )

    ## ---------------------------------------------
    ## gencube info
    ## ---------------------------------------------
    # gencube info subparser
    info_desc = (
        "Resubmit email and NCBI API key for use with NCBI's Entrez Utilities (E-Utilities)"
    )
    parser_info = subparsers.add_parser(
        'info', description=info_desc, help=info_desc, add_help=True,
        formatter_class=argparse.RawTextHelpFormatter
        )

    ## Define return values
    args = parser.parse_args()

    # Gencube genome
    if args.command == 'genome':
        if not args.keywords:
            parser_genome.print_help()
        else:
            genome(
                keywords=args.keywords,
                level=args.level,
                refseq=args.refseq,
                ucsc=args.ucsc,
                latest=args.latest,
                metadata=args.metadata,
                download=args.download,
                database=args.database,
                masking=args.masking,
                chr_style=args.chr_style,
                compresslevel=args.compresslevel,
                recursive=args.recursive,
                )
    # Gencube geneset
    elif args.command == 'geneset':
        if not args.keywords:
            parser_geneset.print_help()
        else:
            geneset(
                keywords=args.keywords,
                level=args.level,
                refseq=args.refseq,
                ucsc=args.ucsc,
                latest=args.latest,
                metadata=args.metadata,
                download=args.download,
                chr_style=args.chr_style,
                recursive=args.recursive,
                )
    # Gencube annotation
    elif args.command == 'annotation':
        if not args.keywords:
            parser_annotation.print_help()
        else:
            annotation(
                keywords=args.keywords,
                level=args.level,
                refseq=args.refseq,
                ucsc=args.ucsc,
                latest=args.latest,
                metadata=args.metadata,
                download=args.download,
                chr_style=args.chr_style,
                recursive=args.recursive,
                )
    # Gencube sequence
    elif args.command == 'sequence':
        if not args.keywords:
            parser_sequence.print_help()
        else:
            sequence(
                keywords=args.keywords,
                level=args.level,
                refseq=args.refseq,
                ucsc=args.ucsc,
                latest=args.latest,
                metadata=args.metadata,
                download=args.download,
                recursive=args.recursive,
                )
    # Gencube crossgenome
    elif args.command == 'crossgenome':
        if not args.keywords:
            parser_crossgenome.print_help()
        else:
            crossgenome(
                keywords=args.keywords,
                level=args.level,
                refseq=args.refseq,
                ucsc=args.ucsc,
                latest=args.latest,
                metadata=args.metadata,
                download=args.download,
                recursive=args.recursive,
                )
    # Gencube seqmeta
    elif args.command == 'seqmeta':
        if not args.keywords and not args.organism and not args.strategy and not args.source and not args.platform and \
        not args.selection and not args.filter and not args.properties and not args.layout and not args.access and \
        not args.bioproject and not args.biosample and not args.accession and not args.title and not args.author and \
        not args.publication and not args.modification and not args.readlength and not args.mbases and \
        not args.textword and not args.exclude and not args.detail and not args.metadata and not args.url:
            parser_seqmeta.print_help()
        else:
            seqmeta(
                keywords=args.keywords,
                organism=args.organism,
                strategy=args.strategy,
                source=args.source,
                platform=args.platform,
                selection=args.selection,
                filter=args.filter,
                properties=args.properties,
                layout=args.layout,
                access=args.access,
                bioproject=args.bioproject,
                biosample=args.biosample,
                accession=args.accession,
                title=args.title,
                author=args.author,
                publication=args.publication,
                modification=args.modification,
                readlength=args.readlength,
                mbases=args.mbases,
                textword=args.textword,
                exclude=args.exclude,
                detail=args.detail,
                metadata=args.metadata,
                url=args.url,
                LS_STRATEGY=LS_STRATEGY,
                LS_SOURCE=LS_SOURCE,
                LS_PLATFORM=LS_PLATFORM,
                LS_SELECTION=LS_SELECTION,
                LS_FILTER=LS_FILTER,
                LS_PROPERTIES=LS_PROPERTIES,
        )
    # Gencube info
    elif args.command == 'info':
        info(info_save)

    else:
        parser.print_help()
