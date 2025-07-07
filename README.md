# <img src="https://github.com/snu-cdrc/gencube/blob/main/figures/logo.png?raw=true" alt="gencube logo" width="200">


<!-- 1. GitHub Version Badge: -->
<!-- 2. PyPI Version Badge: -->
<!-- 3. Supported Python Versions Badge: -->
<!-- 6. PyPI Downloads Badge: -->
<!-- 8. License Badge: -->
[![pypi version](https://img.shields.io/pypi/v/gencube)](https://pypi.org/project/gencube/)
![python versions](https://img.shields.io/pypi/pyversions/gencube)
[![license](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
<br>
![linux](https://img.shields.io/badge/platform-linux--x86__64-orange)
![macos x86_64](https://img.shields.io/badge/platform-macOSX--x86__64-orange)
![macos arm64](https://img.shields.io/badge/platform-macOSX--arm64-orange)

<!-- 4. GitHub Actions CI Status Badge:
![status](https://github.com/keun-hong/gencube/workflows/CI/badge.svg) -->
<!-- 5. Codecov Badge:
[![codecov](https://codecov.io/gh/keun-hong/gencube/branch/master/graph/badge.svg)](https://codecov.io/gh/keun-hong/gencube) -->
<!-- 7. Documentation Badge:
[![docs](https://readthedocs.org/projects/gencube/badge/?version=latest)](https://gencube.readthedocs.io/en/latest/?badge=latest) -->

### Centralized retrieval and integration of multi-omics resources from leading databases

[**Keun Hong Son**](https://keun-hong.github.io/)<sup>1,2,3</sup> and [**Je-Yoel Cho**](https://vetbio.snu.ac.kr/)<sup>1,2,3</sup>

<sup>1</sup> Department of Biochemistry, College of Veterinary Medicine, Seoul National University, Seoul, Korea<br>
<sup>2</sup> Comparative Medicine and Disease Research Center (CDRC), Science Research Center (SRC), Seoul National University, Seoul, Korea<be>
### Publication and Resources
📑 **Paper**:&nbsp;&nbsp;*[**Bioinformatics, 2025**](https://doi.org/10.1093/bioinformatics/btaf128)*<br><br>
🎬 **Video**:&nbsp;&nbsp;[**Save Hundreds of Hours in Your Multi-omics Research!**](https://youtu.be/ftQP0_Dv9yI?si=6_j0FCQxTgPXatqY)<br>
📌 **Blog**:&nbsp;&nbsp;&nbsp;&nbsp;[**Step-by-step guide to developing a Python-based tool**](https://keun-hong.github.io/computerscience/python-tool-dev/)<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[**Step-by-Step guide to registering a package on Bioconda**](https://keun-hong.github.io/computerscience/bioconda-upload/)<br><br>
🛠️**Update**: [**Version History**](https://github.com/snu-cdrc/gencube/blob/main/HISTORY.md)


---

**Your interest and contributions help Gencube evolve to meet your needs!**<br>
If you have any questions, ideas, or suggestions, please share them on our [**Issues page**](https://github.com/snu-cdrc/gencube/issues)🔥.<br>
We'd love to hear from you!

---
`gencube` enables researchers to search for, download, and unify genome assemblies and diverse types of annotations, and retrieve integrated metadata for high-throughput sequencing-based multi-omics resources suitable for specific requirements.

![gencube_overview](https://github.com/snu-cdrc/gencube/blob/main/figures/gencube_overview.jpg?raw=true)

### Databases accessed from gencube
- [GenBank](https://www.ncbi.nlm.nih.gov/genbank/): NCBI GenBank Nucleotide Sequence Database
- [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/): NCBI Reference Sequence Database
- [GenArk](https://hgdownload.soe.ucsc.edu/hubs/): UCSC Genome Archive
- [Ensembl Beta (formerly Rapid Release)](https://beta.ensembl.org/): Ensembl genome browser that provides frequent updates for newly sequenced species
- [Zoonomia TOGA](https://zoonomiaproject.org/the-data/): Tool to infer Orthologs from Genome Alignments
- [INSDC](https://www.insdc.org/): International Nucleotide Sequence Database Collaboration
- [SRA](https://www.ncbi.nlm.nih.gov/sra): NCBI Sequence Read Archive
- [ENA](https://www.ebi.ac.uk/ena/browser/home): EMBL-EBI European Nucleotide Archive
- [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html): DNA Data Bank of Japan

### Detailed information of each database
- [GenBank & RefSeq README.txt](https://ftp.ncbi.nlm.nih.gov/genomes/all/README.txt) - `genome`, `geneset`, `sequence`
- [UCSC GenArk paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03057-x) - `genome`, `geneset`, `annotation`
- [Ensembl Beta Help](https://beta.ensembl.org/help) & [Ensembl 2023 paper](https://academic.oup.com/nar/article/51/D1/D933/6786199?login=false) - `genome`, `geneset`, `sequence`, `crossgenome`
- [Zoonomia TOGA README.txt](https://genome.senckenberg.de/download/TOGA/README.txt) & [Paper](https://www.science.org/doi/10.1126/science.abn3107) - `geneset`, `crossgenome`
- [Search in SRA Entrez](https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/), [Entrez Help](https://www.ncbi.nlm.nih.gov/books/NBK3837/) & [SRA Advanced Search Builder](https://www.ncbi.nlm.nih.gov/sra/advanced) - `seqmeta`
<br>

## Installation
The latest release can be installed with
```bash
$ pip install gencube
```
Alternative
```bash
$ conda install bioconda::gencube
```
<br>

## Email and NCBI API key for E-utilities
When you first run `gencube`, you'll be prompted for your email and NCBI API key, which are saved in the `.gencube_entrez_info` file in your home directory for future use.

All gencube key subcommands use `NCBI's Entrez Utilities (E-Utilities)`, requiring an email. Without an NCBI API key, you can make 3 requests per second; with an NCBI API key, this limit increases to 10 requests per second. If you submit your NCBI API key, you can perform tasks at more than three times the speed when using the seqmeta subcommand, especially when fetching metadata. If possible, it is recommended to submit your API key.
```plaintext
$ gencube
Email address:
NCBI API key (type 'no' to skip):
```
To update the submitted information, run the following command.
```plaintext
$ gencube info
```
<br>

## Tutorials
`gencube` consists of six main subcommands excluding `info`
```plaintext
$ gencube
usage: gencube [-h] {genome,geneset,annotation,sequence,crossgenome,seqmeta,info} ...

gencube v1.0.0

positional arguments:
  {genome,geneset,annotation,sequence,crossgenome,seqmeta,info}
    genome              Search, download, and modify chromosome labels for genome assemblies
    geneset             Search, download, and modify chromosome labels for genesets (gene annotations)
    annotation          Search, download, and modify chromosome labels for various genome annotations, such as gaps and repeats
    sequence            Search and download sequence data of genesets
    crossgenome         Search and download comparative genomics data, such as homology, and codon or protein alignments
    seqmeta             Search, retrive, and integrate metadata of experimental sequencing data
    info                Resubmit email and NCBI API key for use with NCBI's Entrez Utilities (E-Utilities)

options:
  -h, --help            show this help message and exit
```
<img src="https://github.com/snu-cdrc/gencube/blob/main/figures/data_type.jpg?raw=true" width="700">

---
### 0. The positional argument and options shared among the `genome`, `geneset`, `sequence`, `annotation`, and `crossgenome` subcommand
When using the above five subcommands, it's important to find genome assemblies required for personal research.
Below are the positional argument and options shared by the these subcommands to browse and search for specific genome assemblies.

```plaintext
positional arguments:
  keywords              Taxonomic names to search for genomes
                        You can provide various forms such as species names or accession numbers
                        Examples: homo_sapiens, human, GCF_000001405.40, GCA_000001405.29, GRCh38, hg38

                        Multiple names can be combined and will be merged in the search results
                        To specify multiple names, separate them with spaces

options:
  -h, --help            show this help message and exit
  -v level, --level level
                        Specify the genome assembly level (default: complete,chromosome)
                        complete   : Fully assembled genomes
                        chromosome : Assembled at the chromosome level
                        scaffold   : Assembled into scaffolds, but not to the chromosome level
                        contig     : Contiguous sequences without gaps

  -r, --refseq          Show genomes that have RefSeq accession (GCF_* format)
  -u, --ucsc            Show genomes that have UCSC name
  -l, --latest          Show genomes corresponding to the latest version
```
#### Examples
```bash
# Search using scientific or common name
# It is recommended to use the scientific name for a more precise search.
$ gencube genome homo_sapiens
$ gencube genome canis_lupus_familiaris
$ gencube genome human # less accurate
$ gencube genome dog   # less accurate

# Search using assembly name
$ gencube genome GRCh38
$ gencube genome grch38 # Case sensitivity is not an issue
$ gencube genome GRCm39 GRCm38 # Multiple keywords can be searched simultaneously

# Search using UCSC name
$ gencube genome hg38
$ gencube genome hg38 hg19
$ gencube genome mm39 mm10
$ gencube genome canfam4 canfam5 canfam6

# Search using GenBank (GCA_*) or RefSeq (GCF_*) accession
$ gencube genome GCA_021950905.1
$ gencube genome GCF_000001405.40
$ gencube genome GCF_000001405.40 GCA_021950905.1

# Show searched genomes corresponding to all genome assembly levels
$ gencube genome homo_sapiens --level complete,chromosome # default
$ gencube genome homo_sapiens --level complete,chromosome,scaffold,contig
$ gencube genome homo_sapiens --level scaffold,contig

# Only show genomes that have RefSeq accession and UCSC name, and correspond to the latest version
$ gencube genome homo_sapiens --refseq --ucsc --latest
```
#### Example output displayed in the terminal
🔥 **$ gencube genome GCF_000001405.40 GCA_021950905.1**
```plaintext
# Search assemblies in NCBI database
  Keyword: ['GCF_000001405.40', 'GCA_021950905.1']

  Total 3 genomes are searched.

# Convert JSON to dataframe format.
  Filter options
  Level:   ['Complete', 'Chromosome']
  RefSeq:  False
  UCSC:    False
  Latest:  False

# Check accessibility to GenArk, Ensembl Beta
  UCSC GenArk  : 4167 genomes across 2813 species
  Ensembl Rapid: 2272 genomes across 1522 species

+----+------------------------+---------+------------+------------+------------------+--------+----------+-----------+
|    | Assembly name          |   Taxid | Release    | Level      | NCBI             | UCSC   | GenArk   | Ensembl   |
+====+========================+=========+============+============+==================+========+==========+===========+
|  0 | HG002.mat.cur.20211005 |    9606 | 2022/02/04 | Chromosome | GCA_021951015.1  |        | v        | v         |
+----+------------------------+---------+------------+------------+------------------+--------+----------+-----------+
|  1 | HG002.pat.cur.20211005 |    9606 | 2022/02/04 | Chromosome | GCA_021950905.1  |        | v        | v         |
+----+------------------------+---------+------------+------------+------------------+--------+----------+-----------+
|  2 | GRCh38.p14             |    9606 | 2022/02/03 | Chromosome | GCF_000001405.40 | hg38   | v        | v         |
+----+------------------------+---------+------------+------------+------------------+--------+----------+-----------+
```
<br>

---
### 1. `genome` subcommand
**Search, download, and modify chromosome labels for genome assemblies**<br>
You can download genome data in FASTA format from four different databases (GenBank, RefSeq, GenArk, Ensembl Beta). Each database uses a different soft-masking method, and you can selectively download the data as needed. You can also download unmasked and hard-masked genomes from the Ensembl Beta database.
```plaintext
options:
  -m, --metadata        Save metadata for the searched genomes
  -d, --download        Download "fasta" formatted genome file
  -db types, --database types
                        Database where genome file is downloaded (default: refseq)
                        Default is from the RefSeq database
                        If not available, download from the GenBank database
                        genbank : by NCBI GenBank
                        refseq  : by NCBI RefSeq
                        genark  : by UCSC GenArk
                        ensembl : by Ensembl Beta
  -c type, --chr_style type
                        Chromosome label style used in the download file (default: ensembl)
                        ensembl : 1, 2, X, MT & unknowns (GenBank IDs)
                        gencode : chr1, chr2, chrX, chrM & unknowns (GenBank IDs)
                        ucsc    : chr1, chr2, chrX, chrM & unknowns (UCSC-specific IDs)
                                  !! Limited use if UCSC IDs are not issued
                        raw     : Uses raw file labels without modification
                                 - NCBI GenBank: CM_* or other-form IDs
                                 - NCBI RefSeq : NC_*, NW_* or other-form IDs
                                 - GenArk      : GenBank or RefSeq IDs
                                 - Ensembl     : Ensembl IDs
  -mk type, --masking type
                        Masking type for output data (default: soft)
                        soft : soft-masked
                        hard : hard-masked
                        none : unmasked
  -cl 1-9, --compresslevel 1-9
                        Compression level for output data (default: 6)
                        Lower numbers are faster but have lower compression

  --recursive           Download files regardless of their presence only if integrity check is not possible
```
#### Examples
```bash
# Download the full information metadata of searched genomes
$ gencube genome homo_sapiens --metadata
$ gencube genome canis_lupus_familiaris --metadata

# Download genome file under the default conditions (RefSeq or GenBank)
$ gencube genome GCF_011100685.1 --download
$ gencube genome GCF_011100685.1 -d

# Download genome file from a specific database
$ gencube genome GCF_011100685.1 --download --database ensembl
# Download multiple genomes from various databases
$ gencube genome GCF_011100685.1 --download --database refseq,genark,ensembl

# Change the chromosome labels
$ gencube genome GCF_011100685.1 --download --chr_style ensembl # default
$ gencube genome GCF_011100685.1 --download --chr_style gencode
$ gencube genome GCF_011100685.1 --download --chr_style ucsc

# Change the masking type
$ gencube genome GCF_011100685.1 --download --masking soft # default
$ gencube genome GCF_011100685.1 --download --masking hard
$ gencube genome GCF_011100685.1 --download --masking none

# Set the compression level of the file to 2.
$ gencube genome GCF_011100685.1 --download --compresslevel 6 # default
$ gencube genome GCF_011100685.1 --download --compresslevel 1
```
#### Example output displayed in the terminal
🔥 **$ gencube genome GCF_011100685.1 --download --chr_style gencode**
```plaintext
# Search assemblies in NCBI database
  Keyword: ['GCF_011100685.1']

  Total 1 genomes is searched

# Filter genomes based on the following criteria
  Level:   ['Complete', 'Chromosome']
  RefSeq:  False
  UCSC:    False
  Latest:  False

# Check accessibility to GenArk, Ensembl Beta
  UCSC GenArk  : 5396 genomes across 3527 species
  Ensembl Beta : 2938 genomes across 2042 species

+----+-----------------+---------+------------+------------+-----------------+---------+----------+-----------+
|    | Assembly name   |   Taxid | Release    | Level      | NCBI            | UCSC    | GenArk   | Ensembl   |
+====+=================+=========+============+============+=================+=========+==========+===========+
|  0 | UU_Cfam_GSD_1.0 |    9615 | 2020/03/10 | Chromosome | GCF_011100685.1 | canFam4 | v        | v         |
+----+-----------------+---------+------------+------------+-----------------+---------+----------+-----------+

[GCA_011100685.1 / GCF_011100685.1 / UU_Cfam_GSD_1.0]
- refseq
  Assembly Report was already downloaded
  Canis_lupus_familiaris-UU_Cfam_GSD_1.0-refseq.sm.fa.gz: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████| 742M/742M [05:10<00:00, 2.50MB/s]

# Change chromosome names and masking method: gencode-style & soft-masked
[GCA_011100685.1 / GCF_011100685.1 / UU_Cfam_GSD_1.0]
  Downloaded genome: ['refseq.sm']
  - refseq.sm
    Modify chromosome names
    Write compressed fasta file (compresslevel: 6)
    Processing time: 436 seconds

  !! If the file appears to have any problems, please delete it and retry the process
```
<br>

---
### 2. `geneset` subcommand
**Search, download, and modify chromosome labels for genesets (gene annotations)**
```plaintext
options:
  -m, --metadata        Save metadata for the searched genesets
  -d types, --download types
                        Type of gene set
                        refseq_gtf    : RefSeq gene set (GTF format)
                        refseq_gff    : RefSeq gene set (GFF)
                        gnomon        : RefSeq Gnomon gene prediction (GFF)
                        cross         : RefSeq Cross-species alignments (GFF)
                        same          : RefSeq Same-species alignments (GFF)
                        augustus      : GenArk Augustus gene prediction (GFF)
                        xenoref       : GenArk XenoRefGene (GFF)
                        genark_ref    : GenArk RefSeq gene models (GFF)
                        ensembl_gtf   : Ensembl Beta gene set (GTF)
                        ensembl_gff   : Ensembl Beta gene set (GFF)
                        toga_gtf      : Zoonomia TOGA gene set (GTF)
                        toga_bed      : Zoonomia TOGA gene set (BED)
                        toga_pseudo   : Zoonomia TOGA processed pseudogenes (BED)
  -c type, --chr_style type
                        Chromosome label style used in the download file (default: ensembl)
                        ensembl : 1, 2, X, MT & unknowns (GenBank IDs)
                        gencode : chr1, chr2, chrX, chrM & unknowns (GenBank IDs)
                        ucsc    : chr1, chr2, chrX, chrM & unknowns (UCSC-specific IDs)
                                  !! Limited use if UCSC IDs are not issued
                        raw     : Uses raw file labels without modification
                                 - NCBI GenBank: CM_* or other-form IDs
                                 - NCBI RefSeq : NC_*, NW_* or other-form IDs
                                 - GenArk      : GenBank or RefSeq IDs
                                 - Ensembl     : Ensembl IDs

  --recursive           Download files regardless of their presence only if integrity check is not possible
```

#### Examples
```bash
# search usable and accessible data
$ gencube geneset canis_lupus_familiaris
$ gencube geneset GCF_011100685.1

# Download the full information metadata of searched genesets
$ gencube geneset canis_lupus_familiaris --metadata
$ gencube geneset GCF_011100685.1 --metadata

# Download specific geneset file from a database
$ gencube geneset GCF_011100685.1 --download refseq_gtf
$ gencube geneset GCF_011100685.1 --download augustus
$ gencube geneset GCF_011100685.1 --download toga_gtf

# Download multiple genesets from various databases
$ gencube geneset GCF_011100685.1 --download refseq_gtf,augustus,toga_gtf
```
#### Example output displayed in the terminal
🔥 **$ gencube geneset canis_lupus_familiaris**
```plaintext
# Search assemblies in NCBI database
  Keyword: ['canis_lupus_familiaris']

  Total 34 genomes are searched

# Filter genomes based on the following criteria
  Level:   ['Complete', 'Chromosome']
  RefSeq:  False
  UCSC:    False
  Latest:  False

# Check accessibility to GenArk, Ensembl Beta and Zoonomia server
  UCSC GenArk  : 5396 genomes across 3527 species
  Ensembl Beta : 2938 genomes across 2042 species
  Zoonomia TOGA: 1077 genomes across 985 species

+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|    | Assembly name       |   Taxid | Release    | Level      | NCBI            | UCSC    | GenArk   | Ensembl   | Zoonomia   |
+====+=====================+=========+============+============+=================+=========+==========+===========+============+
|  0 | LK_Cfam_Beagle_1.1  |    9615 | 2024/10/30 | Chromosome | GCA_044048985.1 |         |          |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  1 | ASM4364393v1        |    9615 | 2024/10/24 | Chromosome | GCA_043643935.1 |         |          |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  2 | UNSW_CanFamBas_1.2  |    9615 | 2021/02/23 | Chromosome | GCA_013276365.2 |         | v        |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  3 | UMICH_Zoey_3.1      |    9615 | 2019/05/30 | Chromosome | GCF_005444595.1 | canFam5 |          | v         | v          |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  4 | ROS_Cfam_1.0        |    9615 | 2020/09/03 | Chromosome | GCF_014441545.1 |         | v        | v         | v          |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  5 | Yella_v2            |    9615 | 2023/09/07 | Chromosome | GCA_031165255.1 |         | v        |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  6 | CA611_1.0           |    9615 | 2023/08/31 | Chromosome | GCA_031010295.1 |         | v        |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  7 | OD_1.0              |    9615 | 2023/08/31 | Chromosome | GCA_031010635.1 |         | v        |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  8 | BD_1.0              |    9615 | 2023/08/31 | Chromosome | GCA_031010765.1 |         | v        |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
|  9 | Lak_Megaphage_Wal-2 |    9615 | 2024/06/07 | Chromosome | GCA_964164975.1 |         |          |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
| 10 | Lak_Megaphage_Wal-1 |    9615 | 2024/06/07 | Chromosome | GCA_964162815.1 |         |          |           |            |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
| 11 | UU_Cfam_GSD_1.0     |    9615 | 2020/03/10 | Chromosome | GCF_011100685.1 | canFam4 | v        | v         | v          |
+----+---------------------+---------+------------+------------+-----------------+---------+----------+-----------+------------+
...

# Check accessible data in databases
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|    | Assembly name       | UCSC    | RefSeq                        | GenArk                | Ensembl           | Zoonomia     |
+====+=====================+=========+===============================+=======================+===================+==============+
|  0 | LK_Cfam_Beagle_1.1  |         |                               |                       |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  1 | ASM4364393v1        |         |                               |                       |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  2 | UNSW_CanFamBas_1.2  |         |                               | augustus, xenoref      |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  3 | UMICH_Zoey_3.1      | canFam5 | gtf, gff, gnomon, cross, same |                       | ensembl: gtf, gff | human, mouse |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  4 | ROS_Cfam_1.0        |         | gtf, gff, gnomon, cross, same | augustus, xenoref, ref | ensembl: gtf, gff | human        |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  5 | Yella_v2            |         |                               | augustus, xenoref      |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  6 | CA611_1.0           |         |                               | augustus, xenoref      |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  7 | OD_1.0              |         |                               | augustus, xenoref      |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  8 | BD_1.0              |         |                               | augustus, xenoref      |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
|  9 | Lak_Megaphage_Wal-2 |         |                               |                       |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
| 10 | Lak_Megaphage_Wal-1 |         |                               |                       |                   |              |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
| 11 | UU_Cfam_GSD_1.0     | canFam4 | gtf, gff, gnomon, cross, same | augustus, xenoref, ref | ensembl: gtf, gff | human, mouse |
+----+---------------------+---------+-------------------------------+-----------------------+-------------------+--------------+
...
```
<br>

---
### 3. `annotation` subcommand
**Search, download, and modify chromosome labels for various genome annotations, such as gaps and repeats**
```plaintext
options:
  -m, --metadata        Save metadata for the searched annotations
  -d types, --download types
                        Download annotation file.
                        gap : Genomic gaps - AGP defined (bigBed format)
                        sr   : Simple tandem repeats by TRF (bigBed)
                        td   : Tandem duplications (bigBed)
                        wm   : Genomic intervals masked by WindowMasker + SDust (bigBed)
                        rmsk : Repeated elements annotated by RepeatMasker (bigBed)
                        cpg  : CpG Islands - Islands < 300 bases are light green (bigBed)
                        gc   : GC percent in 5-Base window (bigWig)
  -c type, --chr_style type
                        Chromosome label style used in the download file (default: ensembl)
                        ensembl : 1, 2, X, MT & unknowns (GenBank IDs)
                        gencode : chr1, chr2, chrX, chrM & unknowns (GenBank IDs)
                        ucsc    : chr1, chr2, chrX, chrM & unknowns (UCSC-specific IDs)
                                  !! Limited use if UCSC IDs are not issued
                        raw     : Uses raw file labels without modification
                                 - NCBI GenBank: CM_* or other-form IDs
                                 - NCBI RefSeq : NC_*, NW_* or other-form IDs
                                 - GenArk      : GenBank or RefSeq IDs
                                 - Ensembl     : Ensembl IDs

  --recursive           Download files regardless of their presence only if integrity check is not possible
```

#### Examples
```bash
# search usable and accessible data
$ gencube annotation canis_lupus_familiaris
$ gencube annotation GCF_011100685.1

# Download the full information metadata of searched annotations
$ gencube annotation canis_lupus_familiaris --metadata
$ gencube annotation GCF_011100685.1 --metadata

# Download specific annotation file from a database
$ gencube annotation GCF_011100685.1 --download sr
$ gencube annotation GCF_011100685.1 --download td
$ gencube annotation GCF_011100685.1 --download rmsk

# Download multiple annotations from various databases
$ gencube annotation GCF_011100685.1 --download sr,td,rmsk,gc
```
#### Example output displayed in the terminal
🔥 **$ gencube annotation canis_lupus_familiaris**
```plaintext
# Check accessible data in databases
+----+---------------------+---------+------------------------------------+
|    | Assembly name       | UCSC    | GenArk                             |
+====+=====================+=========+====================================+
|  0 | LK_Cfam_Beagle_1.1  |         |                                    |
+----+---------------------+---------+------------------------------------+
|  1 | ASM4364393v1        |         |                                    |
+----+---------------------+---------+------------------------------------+
|  2 | UNSW_CanFamBas_1.2  |         | gap, cpgisland, gc, rm, sr, td, wm |
+----+---------------------+---------+------------------------------------+
|  3 | UMICH_Zoey_3.1      | canFam5 |                                    |
+----+---------------------+---------+------------------------------------+
|  4 | ROS_Cfam_1.0        |         | gap, cpgisland, gc, rm, sr, td, wm |
+----+---------------------+---------+------------------------------------+
|  5 | Yella_v2            |         | gap, cpgisland, gc, rm, sr, td, wm |
+----+---------------------+---------+------------------------------------+
|  6 | CA611_1.0           |         | gap, cpgisland, gc, rm, sr, td, wm |
+----+---------------------+---------+------------------------------------+
|  7 | OD_1.0              |         | gap, cpgisland, gc, rm, sr, td, wm |
+----+---------------------+---------+------------------------------------+
|  8 | BD_1.0              |         | gap, cpgisland, gc, rm, sr, td, wm |
+----+---------------------+---------+------------------------------------+
|  9 | Lak_Megaphage_Wal-2 |         |                                    |
+----+---------------------+---------+------------------------------------+
| 10 | Lak_Megaphage_Wal-1 |         |                                    |
+----+---------------------+---------+------------------------------------+
| 11 | UU_Cfam_GSD_1.0     | canFam4 | gap, cpgisland, gc, rm, sr, td, wm |
+----+---------------------+---------+------------------------------------+
...
```
<br>

---
### 4. `sequence` subcommand
**Search and download sequence data of genesets**
```plaintext
options:
  -m, --metadata        Save metadata for the searched sequence data
  -d types, --download types
                        Download "fasta" formatted sequence file
                        1. Nucleotide sequences:
                           refseq_rna         : Accessioned RNA sequences annotated on the genome assembly
                           refseq_rna_genomic : RNA features based on the genome sequence
                           refseq_cds_genomic : CDS features based on the genome sequence
                           refseq_pseudo      : Pseudogene and other gene regions without transcribed RNA or translated protein products
                           ensembl_cdna       : Ensembl Beta cDNA sequences of transcripts
                           ensembl_repeat     : Ensembl repeat modeler sequences
                        2. Protein sequences:
                           refseq_pep         : Accessioned protein sequences annotated on the genome assembly
                           refseq_pep_cds     : CDS features translated into protein sequences
                           ensembl_pep        : Ensembl Beta protein sequences

  --recursive           Download files regardless of their presence only if integrity check is not possible
```

#### Examples
```bash
# search usable and accessible data
$ gencube sequence canis_lupus_familiaris
$ gencube sequence GCF_011100685.1

# Download the full information metadata of searched sequence data
$ gencube sequence canis_lupus_familiaris --metadata
$ gencube sequence GCF_011100685.1 --metadata

# Download specific sequence file from a database
$ gencube sequence GCF_011100685.1 --download refseq_rna
$ gencube sequence GCF_011100685.1 --download ensembl_cdna
$ gencube sequence GCF_011100685.1 --download ensembl_pep

# Download multiple sequence data from various databases
$ gencube sequence GCF_011100685.1 --download refseq_rna,ensembl_cdna,refseq_pep,ensembl_pep
```
#### Example output displayed in the terminal
🔥 **$ gencube sequence canis_lupus_familiaris**
```plaintext
# Check accessible data in databases
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|    | Assembly name       | UCSC    | RefSeq                                              | Ensembl            |
+====+=====================+=========+=====================================================+====================+
|  0 | LK_Cfam_Beagle_1.1  |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  1 | ASM4364393v1        |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  2 | UNSW_CanFamBas_1.2  |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  3 | UMICH_Zoey_3.1      | canFam5 | rna, rna_genomic, cds_genomic, pseudo, pep, pep_cds | ensembl: cdna, pep |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  4 | ROS_Cfam_1.0        |         | rna, rna_genomic, cds_genomic, pseudo, pep, pep_cds | ensembl: cdna, pep |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  5 | Yella_v2            |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  6 | CA611_1.0           |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  7 | OD_1.0              |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  8 | BD_1.0              |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
|  9 | Lak_Megaphage_Wal-2 |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
| 10 | Lak_Megaphage_Wal-1 |         |                                                     |                    |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
| 11 | UU_Cfam_GSD_1.0     | canFam4 | rna, rna_genomic, cds_genomic, pseudo, pep, pep_cds | ensembl: cdna, pep |
+----+---------------------+---------+-----------------------------------------------------+--------------------+
...
```
<br>

---
### 5. `crossgenome` subcommand
**Search and download comparative genomics data, such as homology, and codon or protein alignment**
```plaintext
options:
  -m, --metadata        Save metadata for the searched comparative genomics data
  -d types, --download types
                        ensembl_homology   : Homology data from Ensembl Beta,
                                             detailing gene orthology relationships across species
                        toga_homology      : Homology data from TOGA, providing predictions of
                                             orthologous genes based on genome alignments
                        toga_align_codon   : Codon alignment data from TOGA, showing aligned codon
                                             sequences between reference and query species
                        toga_align_protein : Protein alignment data from TOGA, detailing aligned
                                             protein sequences between reference and query species
                        toga_inact_mut     : List of inactivating mutations from TOGA, identifying
                                             mutations that disrupt gene function

  --recursive           Download files regardless of their presence only if integrity check is not possible
```

#### Examples
```bash
# search usable and accessible data
$ gencube crossgenome canis_lupus_familiaris
$ gencube crossgenome GCF_011100685.1

# Download the full information metadata of searched comparative genomics data
$ gencube crossgenome canis_lupus_familiaris --metadata
$ gencube crossgenome GCF_011100685.1 --metadata

# Download specific crossgenome file from a database
$ gencube crossgenome GCF_011100685.1 --download toga_homology
$ gencube crossgenome GCF_011100685.1 --download toga_align_codon

# Download multiple crossgenome data
$ gencube crossgenome GCF_011100685.1 --download toga_homology,toga_align_codon
```
#### Example output displayed in the terminal
🔥 **$ gencube crossgenome canis_lupus_familiaris**
```plaintext
# Check accessible data in databases
+----+---------------------+---------+-----------+--------------+
|    | Assembly name       | UCSC    | Ensembl   | Zoonomia     |
+====+=====================+=========+===========+==============+
|  0 | LK_Cfam_Beagle_1.1  |         |           |              |
+----+---------------------+---------+-----------+--------------+
|  1 | ASM4364393v1        |         |           |              |
+----+---------------------+---------+-----------+--------------+
|  2 | UNSW_CanFamBas_1.2  |         |           |              |
+----+---------------------+---------+-----------+--------------+
|  3 | UMICH_Zoey_3.1      | canFam5 | ensembl   | human, mouse |
+----+---------------------+---------+-----------+--------------+
|  4 | ROS_Cfam_1.0        |         | ensembl   | human        |
+----+---------------------+---------+-----------+--------------+
|  5 | Yella_v2            |         |           |              |
+----+---------------------+---------+-----------+--------------+
|  6 | CA611_1.0           |         |           |              |
+----+---------------------+---------+-----------+--------------+
|  7 | OD_1.0              |         |           |              |
+----+---------------------+---------+-----------+--------------+
|  8 | BD_1.0              |         |           |              |
+----+---------------------+---------+-----------+--------------+
|  9 | Lak_Megaphage_Wal-2 |         |           |              |
+----+---------------------+---------+-----------+--------------+
| 10 | Lak_Megaphage_Wal-1 |         |           |              |
+----+---------------------+---------+-----------+--------------+
| 11 | UU_Cfam_GSD_1.0     | canFam4 | ensembl   | human, mouse |
+----+---------------------+---------+-----------+--------------+
...
```
<br>

---
### 6. `seqmeta` subcommand
**Searches for high-throughput sequencing data corresponding to user-specified keywords, retrieves the related sample metadata, and integrates it into experiment-level and study-level tables.**
![seqmeta_scheme](https://github.com/snu-cdrc/gencube/blob/main/figures/seqmeta_scheme.jpg?raw=true)
```plaintext
$ gencube seqmeta
usage: gencube seqmeta [-h] [-o string] [-st string] [-sr string] [-pl string] [-sl string] [-fi string] [-pr string] [-ly string] [-ac string] [-bp string] [-bs string]
                       [-as string] [-ti string] [-at string] [-pd range] [-md range] [-rl range] [-mb string] [-tw string] [-ex keywords] [-d] [-m]
                       [keywords ...]

Search, retrive, and integrate metadata of experimental sequencing data

positional arguments:
  keywords              Keywords to search for sequencing-based experimental data. You can provide various forms
                        Examples: liver, k562, cancer, breast_cancer, etc

                        Multiple keywords can be combined
                        Keywords separated by commas will combine their results
                        Keywords separated by spaces will intersect their results
                        Example: liver,lung cancer,tumor

options:
  -h, --help            show this help message and exit
  -o string, --organism string
                        Scientific name or common name (as found in the NCBI Taxonomy Browser)
                        Example: homo_sapiens or human
  -st string, --strategy string
                        Sequencing strategy:
                        wgs, wga, wxs, targeted_capture, synthetic_long_read, gbs, rad_seq, tn_seq, clone_end, amplicon
                        clone, rna_seq, mrna_seq, ncrna_seq, ribo_seq, rip_seq, mirna_seq, ssrna_seq, est, fl_cdna, atac_seq
                        dnase_hypersensitivity, faire_seq, chip_seq, chip, mre_seq, bisulfite_seq, mbd_seq, medip_seq, hi_c
                        chia_pet, tethered_chromatin_conformation_capture
  -sr string, --source string
                        Source of the biological data:
                        genomic, genomic_single_cell, transcriptomic, transcriptomic_single_cell, metagenomic
                        metatranscriptomic, synthetic, viral_rna, other
  -pl string, --platform string
                        Name of the sequencing platform:
                        abi_solid, bgiseq, capillary, complete_genomics, dnbseq, element, genapsys, genemind, helicos
                        illumina, ion_torrent, ls454, oxford_nanopore, pacbio_smrt, tapestri, ultima, vela_diagnostics
  -sl string, --selection string
                        Library selection methodology:
                        5_methylcytidine_antibody, cage, cdna, cdna_oligo_dt, cdna_randompriming, chip, chip_seq, dnase
                        hmpr, hybrid_selection, inverse_rrna, mbd2_protein_methyl_cpg_binding_domain, mda, mf, mnase, msll
                        oligo_dt, other, padlock_probes_capture_method, pcr, polya, race, random, random_pcr
                        reduced_representation, repeat_fractionation, restriction_digest, rt_pcr, size_fractionation
                        unspecified
  -fi string, --filter string
                        Option to find SRA records that are cross-referenced with other NCBI databases
                        (PubMed, PubMed Central (PMC), Nucleotide, Assembly, and others):
                        sra_all, sra_assembly, sra_bioproject, sra_bioproject_all, sra_biosample, sra_biosample_all, sra_gap
                        sra_gap_all, sra_gds, sra_genome, sra_nuccore, sra_nuccore_alignment, sra_nuccore_wgs, sra_omim
                        sra_pmc, sra_public, sra_pubmed, sra_taxonomy
  -pr string, --properties string
                        Option to narrow search results by controlled-vocabulary library's annotations:
                        aligned_data, cloud_gs, cloud_s3, location_gs_us, location_s3_us_east1, location_s3_us_east_1
                        location_s3_us_west_2, filetype_10x_genomics_bam_file, filetype_ab1, filetype_activ_sars2_vcf
                        filetype_archive, filetype_archive/gzip, filetype_assembled_contigs
                        filetype_assembly/realign_summary, filetype_assembly_of_unidentified_reads, filetype_bai
                        filetype_bam, filetype_bam_header, filetype_basemodification, filetype_complete_genomics
                        filetype_crai, filetype_cram, filetype_fast5, filetype_fasta, filetype_fastq
                        filetype_geo_feature_count, filetype_helicos, filetype_illumina_native, filetype_nanopore
                        filetype_pacbio_base_modification_report, filetype_pacbio_metadata, filetype_pacbio_native
                        filetype_realign_to_de_novo_assembly, filetype_reference_fasta, filetype_run, filetype_run_realign
                        filetype_run_zq, filetype_sff, filetype_solid_native, filetype_source, filetype_sra_lite
                        filetype_sra_normalized, filetype_srf, filetype_tar_archive_of_complete_genomics_tree, filetype_tenx
                        filetype_vcf, filetype_vcf_index, filetype_vdbcache, filetype_vdbcache_zq, filetype_wgmlst_sig
                        filetype_wgmlst_signature, has_data, instrument_454_gs, instrument_454_gs_20, instrument_454_gs_flx
                        instrument_454_gs_flx_titanium, instrument_454_gs_junior, instrument_ab_310_genetic_analyzer
                        instrument_ab_3130_genetic_analyzer, instrument_ab_3130xl_genetic_analyzer
                        instrument_ab_3500_genetic_analyzer, instrument_ab_3500xl_genetic_analyzer
                        instrument_ab_3730_genetic_analyzer, instrument_ab_3730xl_genetic_analyzer
                        instrument_ab_5500_genetic_analyzer, instrument_ab_5500xl_genetic_analyzer
                        instrument_ab_5500xl_w_genetic_analysis_system, instrument_ab_solid_3_plus_system
                        instrument_ab_solid_4_system, instrument_ab_solid_4hq_system, instrument_ab_solid_pi_system
                        instrument_ab_solid_system, instrument_ab_solid_system_2_0, instrument_ab_solid_system_3_0
                        instrument_bgiseq_50, instrument_bgiseq_500, instrument_complete_genomics, instrument_dnbseq_g400
                        instrument_dnbseq_g400_fast, instrument_dnbseq_g50, instrument_dnbseq_t7, instrument_element_aviti
                        instrument_fastaseq_300, instrument_genexus, instrument_genolab_m, instrument_gridion
                        instrument_gs111, instrument_helicos_heliscope, instrument_hiseq_x_five, instrument_hiseq_x_ten
                        instrument_illumina_genome_analyzer, instrument_illumina_genome_analyzer_ii
                        instrument_illumina_genome_analyzer_iix, instrument_illumina_hiscansq
                        instrument_illumina_hiseq_1000, instrument_illumina_hiseq_1500, instrument_illumina_hiseq_2000
                        instrument_illumina_hiseq_2500, instrument_illumina_hiseq_3000, instrument_illumina_hiseq_4000
                        instrument_illumina_hiseq_x, instrument_illumina_hiseq_x_ten, instrument_illumina_iseq_100
                        instrument_illumina_miniseq, instrument_illumina_miseq, instrument_illumina_novaseq_6000
                        instrument_illumina_novaseq_x_plus, instrument_ion_genesudio_s5, instrument_ion_genesudio_s5_plus
                        instrument_ion_genesudio_s5_prime, instrument_ion_s5, instrument_ion_s5_xl
                        instrument_ion_torrent_genexus, instrument_ion_torrent_pgm, instrument_ion_torrent_proton
                        instrument_ion_torrent_s5, instrument_ion_torrent_s5_xl, instrument_mgiseq_2000rs, instrument_minion
                        instrument_nextseq_1000, instrument_nextseq_2000, instrument_nextseq_500, instrument_nextseq_550
                        instrument_onso, instrument_pacbio_rs, instrument_pacbio_rs_ii, instrument_promethion
                        instrument_revio, instrument_sentosa_sq301, instrument_sequel, instrument_sequel_ii
                        instrument_sequel_iie, instrument_tapestri, instrument_ug_100, instrument_unspecified
                        study_type_cancer_genomics, study_type_epigenetics, study_type_exome_sequencing
                        study_type_metagenomics, study_type_other, study_type_pooled_clone_sequencing
                        study_type_population_genomics, study_type_synthetic_genomics, study_type_transcriptome_analysis
                        study_type_transcriptome_sequencing, study_type_whole_genome_sequencing
  -ly string, --layout string
                        Library layout of the sequencing data:
                        paired, single
  -ac string, --access string
                        Data accessibility:
                        public, controlled
  -bp string, --bioproject string
                        BioProject accession in the form of PRJNA#, PRJEB#, or PRJDB#
  -bs string, --biosample string
                        BioSample accession in the form of SAMN#, SAMEA#, or SAMD#
  -as string, --accession string
                        SRA/ENA/DDBJ accession
                        Study with accessions in the form of SRP#, ERP#, or DRP#
                        Sample with accessions in the form of SRS#, ERS#, or DRS#
                        Experiment with accessions in the form of SRX#, ERX#, or DRX#
                        Run with accessions in the form of SRR#, ERR#, or DRR#
  -ti string, --title string
                        Descriptive name of the dataset
  -at string, --author string
                        Researcher or group that submitted the data
                         Example: SON_KH
  -pd range, --publication range
                        Publication Date
                        YYYY.MM.DD : YYYY.MM.DD format
                        Example: 2016, 2016.07, 2016.07.01, 2016.07:2023.02
  -md range, --modification range
                        Modification Date
                        YYYY.MM.DD : YYYY.MM.DD format
                        Example: 2016, 2016.07, 2016.07.01, 2016.07:2023.02
  -rl range, --readlength range
                        Length of the sequencing readsExample: 100 or 100:500
  -mb string, --mbases string
                        Number of mega bases in the SRA Runs
  -tw string, --textword string
                        General search term for finding datasets by specific words in metadata
  -ex keywords, --exclude keywords
                        Exclude the results for the keywords used in this option
                        Example: cell_line,normal,crispr

  -d, --detail          Show the number of searched results for each option and keyword

  -m, --metadata        Save integrated metadata

  -u, --url             Save the file, including the URL address of the raw data (.fastq)

```

#### Examples
```bash
# Search for specific sequencing data for a specific species
$ gencube seqmeta --organism homo_sapiens --strategy rna_seq
$ gencube seqmeta --organism homo_sapiens --strategy chip_seq
$ gencube seqmeta --organism human --strategy chip_seq

# Search for cancer data for specific tissues
$ gencube seqmeta --organism human --strategy chip_seq liver,lung   # liver OR lung
$ gencube seqmeta --organism human --strategy chip_seq cancer,tumor # cancer OR tumor
$ gencube seqmeta --organism human --strategy chip_seq liver,lung cancer,tumor # (liver OR lung) AND (cancer OR tumor)

# Exclude results containing specific keywords
$ gencube seqmeta --organism human --strategy chip_seq liver,lung cancer,tumor --exclude cell_line,crispr # cell_line OR crispr

# Use wild card (*) to search for a broader range of results
$ gencube seqmeta --organism human --strategy chip_seq liver,lung cancer*,tumor* # (liver OR lung) AND (cancer* OR tumor*)

# Use ^ for phrase (not word) search
$ gencube seqmeta --organism human --strategy chip_seq liver,lung cancer,tumor --exclude cell_line^,crispr # "cell line" OR crispr

# Search using accession
$ gencube seqmeta PRJNA838583
$ gencube seqmeta SRP375422
# or specifically
$ gencube seqmeta --bioproject PRJNA838583
$ gencube seqmeta --accession SRP375422

# Search using a custom query
$ gencube seqmeta '(((human[Organism]) AND ("chip seq"[Strategy])) AND ((liver OR lung) AND (cancer OR tumor)))'

# Output the number of search results for each option and keyword
$ gencube seqmeta --organism human --strategy chip_seq --exclude cell_line,crispr liver,lung cancer,tumor --detail

# Save the integrated metadata
$ gencube seqmeta --organism human --strategy chip_seq --exclude cell_line,crispr liver,lung cancer,tumor --metadata
```
#### Using the --url option, you can obtain additional files that include the URLs of the raw files.
🔥 **$ gencube seqmeta -o dog -st chip_seq --metadata --url**
```plaintext
Experiment	Type	URL
SRX23379089	paired	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR277/078/SRR27712678/SRR27712678_1.fastq.gz
SRX23379089	paired	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR277/078/SRR27712678/SRR27712678_2.fastq.gz
SRX23379093	paired	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR277/074/SRR27712674/SRR27712674_1.fastq.gz
SRX23379093	paired	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR277/074/SRR27712674/SRR27712674_2.fastq.gz
SRX1726522	single	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR344/009/SRR3440079/SRR3440079.fastq.gz
SRX1559581	single	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR314/007/SRR3143427/SRR3143427.fastq.gz
SRX1559583	single	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR314/009/SRR3143429/SRR3143429.fastq.gz
...
```
<br>

---
### Structure of output files and directories
🗂️ **Displayed results**
```plaintext
$ tree
├── gencube_output
│   ├── Canis_lupus_familiaris-CanFam3.1_chr_name.txt
│   ├── Canis_lupus_familiaris-Dog10K_Boxer_Tasha_chr_name.txt
│   ├── genome
│   │   ├── Canis_lupus_familiaris-CanFam3.1-refseq.sm.ens-id.fa.gz
│   │   └── Meta_genome_canfam3_complete-chromosome_241125_103030.txt
│   ├── geneset
│   │   ├── Canis_lupus_familiaris-CanFam3.1-refseq.ens-id.gtf.gz
│   │   └── Canis_lupus_familiaris-Dog10K_Boxer_Tasha-ensembl_ensembl.ens-id.gtf.gz
│   ├── annotation
│   │   ├── Canis_lupus_familiaris-CanFam3.1-genark.chrom.sizes.ens-id.txt
│   │   ├── Canis_lupus_familiaris-CanFam3.1-genark.gc5Base.ens-id.bw
│   │   └── Canis_lupus_familiaris-CanFam3.1-genark.simpleRepeat.ens-id.bb
│   └── sequence
│   │   ├── Canis_lupus_familiaris-CanFam3.1-refseq.protein.faa.gz
│   │   └── Canis_lupus_familiaris-CanFam3.1-refseq.rna.fna.gz
│   ├── crossgenome
│   │   ├── Canis_lupus_familiaris-UU_Cfam_GSD_1.0-ensembl_ensembl_homology.tsv.gz
│   │   ├── Canis_lupus_familiaris-UU_Cfam_GSD_1.0-toga_human_homology.tsv.gz
│   │   └── Canis_lupus_familiaris-UU_Cfam_GSD_1.0-toga_mouse_homology.tsv.gz
│   └── seqmeta
│       ├── Meta_seq_dog_241125_112002_experiment_n279.txt
│       ├── Meta_seq_dog_241125_112002_experiment_n279_urls.txt
│       └── Meta_seq_dog_241125_112002_study_n14.txt
└── gencube_raw_download
    ├── bedGraphToBigWig
    ├── bedToBigBed
    ├── bigBedToBed
    ├── bigWigToBedGraph
    ├── Canis_lupus_familiaris-CanFam3.1_assembly_report.txt
    ├── Canis_lupus_familiaris-CanFam3.1-genark.chrom.sizes.txt
    ├── Canis_lupus_familiaris-CanFam3.1-genark.gc5Base.bw
    ├── Canis_lupus_familiaris-CanFam3.1-genark.simpleRepeat.bb
    ├── Canis_lupus_familiaris-CanFam3.1-refseq.gtf.gz
    ├── Canis_lupus_familiaris-CanFam3.1-refseq.sm.fa.gz
    ├── Canis_lupus_familiaris-Dog10K_Boxer_Tasha_assembly_report.txt
    └── Canis_lupus_familiaris-Dog10K_Boxer_Tasha-ensembl_ensembl.gtf.gz
```
🔥 **Commands used to generate the above results**
```bash
# Genome
$ gencube genome canfam3 --donwload
$ gencube genome canfam3 --metadata
# Gene set
$ gencube geneset canfam3 -d refseq_gtf
$ gencube geneset canfam3 -d ensembl_gtf
# Annotation
$ gencube annotation canfam3 -d sr
$ gencube annotation canfam3 -d gc
# Sequence
$ gencube sequence canfam3 -d refseq_rna
$ gencube sequence canfam3 -d refseq_pep
# Cross-genome
$ gencube crossgenome GCF_011100685.1 -d ensembl_homology
$ gencube crossgenome GCF_011100685.1 -d toga_homology
# Seqmeta
$ gencube seqmeta -o dog -st chip_seq --metadata --url
```
<br>

---
### Credits
This package was created with [`Cookiecutter`](https://github.com/audreyr/cookiecutter) and the [`audreyr/cookiecutter-pypackage`](https://github.com/audreyr/cookiecutter-pypackage) project template.
