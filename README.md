Gencube
=======
<!-- 1. GitHub Version Badge: -->
<!-- 2. PyPI Version Badge: -->
<!-- 3. Supported Python Versions Badge: -->
<!-- 6. PyPI Downloads Badge: -->
<!-- 8. License Badge: -->
![github version](https://img.shields.io/badge/Version-1.0.0-informational)
[![pypi version](https://img.shields.io/pypi/v/gencube)](https://pypi.org/project/gencube/)
![python versions](https://img.shields.io/pypi/pyversions/gencube)
[![pypi downloads](https://img.shields.io/pypi/dm/gencube)](https://pypi.org/project/gencube/)
[![license](https://img.shields.io/pypi/l/gencube)](LICENSE)

<!-- 4. GitHub Actions CI Status Badge:
![status](https://github.com/keun-hong/gencube/workflows/CI/badge.svg) -->
<!-- 5. Codecov Badge:
[![codecov](https://codecov.io/gh/keun-hong/gencube/branch/master/graph/badge.svg)](https://codecov.io/gh/keun-hong/gencube) -->
<!-- 7. Documentation Badge:
[![docs](https://readthedocs.org/projects/gencube/badge/?version=latest)](https://gencube.readthedocs.io/en/latest/?badge=latest) -->

### Efficient retrieval, download, and unification of genomic data from leading biodiversity databases

[**Keun Hong Son**](https://keun-hong.github.io/)<sup>1,2,3</sup>, and [**Je-Yoel Cho**](https://vetbio.snu.ac.kr/)<sup>1,2,3</sup>

<sup>1</sup> Department of Biochemistry, College of Veterinary Medicine, Seoul National University, Seoul, Korea<br>
<sup>2</sup> Comparative Medicine and Disease Research Center (CDRC), Science Research Center (SRC), Seoul National University, Seoul, Korea<br>
<sup>3</sup> BK21 PLUS Program for Creative Veterinary Science Research and Research Institute for Veterinary Science, Seoul National University, Seoul, Korea<be>
### Manuscript
[**bioRxiv**]() (uploaded: 2024.07.01)
<!-- Bioinformatics (accepted. 2024.09) -->

---
`gencube` enables researchers to search for, download, and unify genome assemblies and diverse types of annotations, and retrieve metadata for sequencing-based experimental data suitable for specific requirements.

![gencube_overview](https://github.com/snu-cdrc/gencube/blob/main/figures/gencube_overview.jpg?raw=true)

### Databases accessed from gencube
- [GenBank](https://www.ncbi.nlm.nih.gov/genbank/): NCBI GenBank Nucleotide Sequence Database
- [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/): NCBI Reference Sequence Database
- [GenArk](https://hgdownload.soe.ucsc.edu/hubs/): UCSC Genome Archive
- [Ensembl Rapid Release](https://rapid.ensembl.org/index.html): Ensembl genome browser that provides frequent updates for newly sequenced species
- [Zoonomia TOGA](https://zoonomiaproject.org/the-data/): Tool to infer Orthologs from Genome Alignments
- [INSDC](https://www.insdc.org/): International Nucleotide Sequence Database Collaboration
- [SRA](https://www.ncbi.nlm.nih.gov/sra): NCBI Sequence Read Archive
- [ENA](https://www.ebi.ac.uk/ena/browser/home): EMBL-EBI European Nucleotide Archive
- [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html): DNA Data Bank of Japan

### Detailed information of each database
- [GenBank & RefSeq README.txt](https://ftp.ncbi.nlm.nih.gov/genomes/all/README.txt) - `genome`, `geneset`, `sequence`
- [UCSC GenArk paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03057-x) - `genome`, `geneset`, `annotation`
- [Ensembl Rapid Release Help & Docs](https://rapid.ensembl.org/info/index.html) & [Ensembl 2023 paper](https://academic.oup.com/nar/article/51/D1/D933/6786199?login=false) - `genome`, `geneset`, `sequence`, `crossgenome`
- [Zoonomia TOGA README.txt](https://genome.senckenberg.de/download/TOGA/README.txt) & [Paper](https://www.science.org/doi/10.1126/science.abn3107) - `geneset`, `crossgenome`
- [Search in SRA Entrez](https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/), [Entrez Help](https://www.ncbi.nlm.nih.gov/books/NBK3837/) & [SRA Advanced Search Builder](https://www.ncbi.nlm.nih.gov/sra/advanced) - `seqmeta`


## Installation
The latest release can be installed with
```bash
$ pip install gencube
```
Alternative
```bash
$ conda install -c bioconda gencube
```
## Email and NCBI API key for E-utilities
When you first run `gencube`, you'll be prompted for your email and NCBI API key, which are saved in the `.gencube_entrez_info` file in your home directory for future use. To update these details, edit this file.

All gencube subcommands use `E-utilities`, requiring an email. Without an API key, you can make 3 requests per second; with an API key, this limit increases to 10 requests per second. If you hit usage limits, provide an API key.
```plaintext
$ gencube
Email address: 
NCBI API key (type 'no' to skip): 
```

## Tutorials
`gencube` consists of six subcommands
```plaintext
$ gencube
usage: gencube [-h] {genome,geneset,sequence,annotation,crossgenome,seqmeta} ...

gencube v1.0.0

positional arguments:
  {genome,geneset,annotation,sequence,crossgenome,seqmeta}
    genome              Search, download, and modify chromosome labels for genome assemblies
    geneset             Search, download, and modify chromosome labels for genesets (gene annotations)
    annotation          Search, download, and modify chromosome labels for various genome annotations, such as gaps and repeats
    sequence            Search and download sequence data of genesets
    crossgenome         Search and download comparative genomics data, such as homology, and codon or protein alignments
    seqmeta             Search, retrive, and integrate metadata of experimental sequencing data

options:
  -h, --help            show this help message and exit
```
---
### The positional arguments and options shared among the `genome`, `geneset`, `sequence`, `annotation`, and `crossgenome` subcommand
When using the above five subcommands, it's important to find genome assemblies required for personal research.
Below are the positional arguments and options shared by the these subcommands to browse and search for specific genome assemblies.

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
$ gencube genome homo_sapiens canis_lupus_familiaris
$ gencube genome human dog

# Search using assembly name
$ gencube genome T2T-CHM13v2.0 GRCh38

# Search using UCSC name
$ gencube genome hg38 hg19

# Search using GenBank (GCF_*) or RefSeq (GCA_*) accession
$ gencube genome GCF_000001405.40 GCA_021950905.1

# Show searched genomes corresponding to all genome assembly levels
$ gencube genome homo_sapiens --level complete,chromosome,scaffold,contig

# Only show genomes that have RefSeq accession and UCSC name, and correspond to the latest version
$ gencube genome homo_sapiens --refseq --ucsc --latest
```
#### Example output displayed in the terminal
```plaintext
$ gencube genome GCF_000001405.40 GCA_021950905.1

# Search assemblies in NCBI database
  Keyword: ['GCF_000001405.40', 'GCA_021950905.1']

  Total 3 genomes are searched.

# Convert JSON to dataframe format.
  Filter options
  Level:   ['Complete', 'Chromosome']
  RefSeq:  False
  UCSC:    False
  Latest:  False

# Check accessibility to GenArk, Ensembl Rapid Release
  UCSC GenArk  : 4167 genomes across 2813 species
  Ensembl Rapid: 2272 genomes across 1522 species

+----+------------------------+---------+------------+------------------+--------+----------+-----------+
|    | Assembly name          |   Taxid | Release    | NCBI             | UCSC   | GenArk   | Ensembl   |
+====+========================+=========+============+==================+========+==========+===========+
|  0 | HG002.mat.cur.20211005 |    9606 | 2022/02/04 | GCA_021951015.1  |        | v        | v         |
+----+------------------------+---------+------------+------------------+--------+----------+-----------+
|  1 | HG002.pat.cur.20211005 |    9606 | 2022/02/04 | GCA_021950905.1  |        | v        | v         |
+----+------------------------+---------+------------+------------------+--------+----------+-----------+
|  2 | GRCh38.p14             |    9606 | 2022/02/03 | GCF_000001405.40 | hg38   | v        |           |
+----+------------------------+---------+------------+------------------+--------+----------+-----------+
```

---
### `genome`: Search, download, and modify chromosome labels for genome assemblies
You can download genome data in FASTA format from four different databases (GenBank, RefSeq, GenArk, Ensembl Rapid Release). Each database uses a different soft-masking method, and you can selectively download the data as needed. You can also download unmasked and hard-masked genomes from the Ensembl Rapid Release database.
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
                        ensembl : by Ensembl Rapid Release
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
# Download genome files under the default conditions (RefSeq or GenBank)
$ gencube genome GCF_011100685.1 --download
# Download multiple genomes from various databases
$ gencube genome GCF_011100685.1 --download --database refseq,genark,ensembl
# Change the chromosome labels to the GENCODE style and set the compression level of the file to 9.
$ gencube genome GCF_011100685.1 --download --chr_style gencode --compresslevel 9
```
---

### `geneset`: Search, download, and modify chromosome labels for genesets (gene annotations)
```plaintext
options:
  -d types, --download types
                        Type of gene set
                        refseq_gtf    : RefSeq gene set (GTF format)
                        refseq_gff    : RefSeq gene set (GFF)
                        gnomon        : RefSeq Gnomon gene prediction (GFF)
                        cross         : RefSeq Cross-species alignments (GFF)
                        same          : RefSeq Same-species alignments (GFF)
                        agustus       : GenArk Augustus gene prediction (GFF)
                        xenoref       : GenArk XenoRefGene (GFF)
                        genark_ref    : GenArk RefSeq gene models (GFF)
                        ensembl_gtf   : Ensembl Rapid Release gene set (GTF)
                        ensembl_gff   : Ensembl Rapid Release gene set (GFF)
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
gencube geneset GCF_011100685.1

# Download multiple genesets from various databases
$ gencube geneset GCF_011100685.1 --download refseq_gtf,agustus,toga_gtf
```
---

### `annotation`: Search, download, and modify chromosome labels for various genome annotations, such as gaps and repeats
```plaintext
options:
  -d types, --download types
                        Download annotation file.
                        gaps : Genomic gaps - AGP defined (bigBed format)
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
gencube annotation GCF_011100685.1

# Download multiple annotations
gencube annotation GCF_011100685.1 --download sr,td,rmsk,gc
```
---

### `sequence`: Search and download sequence data of genesets
```plaintext
options:
  -d types, --download types
                        Download "fasta" formatted sequence file
                        1. Nucleotide sequences:
                           refseq_rna         : Accessioned RNA sequences annotated on the genome assembly
                           refseq_rna_genomic : RNA features based on the genome sequence
                           refseq_cds_genomic : CDS features based on the genome sequence
                           refseq_pseudo      : Pseudogene and other gene regions without transcribed RNA or translated protein products
                           ensembl_cdna       : Ensembl Rapid Release cDNA sequences of transcripts
                           ensembl_cds        : Ensembl Rapid Release coding sequences (CDS)
                           ensembl_repeat     : Ensembl repeat modeler sequences
                        2. Protein sequences:
                           refseq_pep         : Accessioned protein sequences annotated on the genome assembly
                           refseq_pep_cds     : CDS features translated into protein sequences
                           ensembl_pep        : Ensembl Rapid Release protein sequences
                         
  --recursive           Download files regardless of their presence only if integrity check is not possible
```

#### Examples
```bash
# search usable and accessible data
gencube sequence GCF_011100685.1

# Download multiple genesets from various databases
$ gencube sequence GCF_011100685.1 --download refseq_rna,ensembl_cdna,refseq_pep,ensembl_pep
```
---



### `crossgenome`: Search and download comparative genomics data, such as homology, and codon or protein alignment
```plaintext
options:
  -d types, --download types
                        ensembl_homology   : Homology data from Ensembl Rapid Release,
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
gencube crossgenome GCF_011100685.1

# Download multiple crossgenome data
$ gencube sequence GCF_011100685.1 --download toga_homology,toga_align_codon
```
---

### `seqmeta`: Search, fetch, and integrate metadata of experimental sequencing data
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
```

#### Examples
```bash
# Search for specific sequencing data for a specific species
$ gencube seqmeta --organism human --strategy chip,chip_seq
$ gencube seqmeta --organism homo_sapiens --strategy chip,chip_seq
$ gencube seqmeta --organism homo_sapiens --strategy rna_seq

# Search for cancer data for specific tissues
$ gencube seqmeta --organism human --strategy chip,chip_seq liver,lung cancer,tumor

# Exclude results containing specific keywords
$ gencube seqmeta --organism human --strategy chip,chip_seq --exclude cell_line,crispr liver,lung cancer,tumor

# Use wild card (*) to search for a broader range of results
$ gencube seqmeta --organism human --strategy chip,chip_seq liver,lung cancer*,tumor*

# Use ^ for phrase (not word) search
$ gencube seqmeta --organism human --strategy chip,chip_seq --exclude cell_line^,crispr liver,lung cancer,tumor

# Search using accession
python -m gencube seqmeta PRJNA838583
python -m gencube seqmeta SRP375422
(or specifically)
python -m gencube seqmeta --bioproject PRJNA838583
python -m gencube seqmeta --accession SRP375422

# Search using custom query
python -m gencube seqmeta '(((human[Organism]) AND ("chip"[Strategy] OR "chip seq"[Strategy])) AND ((liver OR lung) AND (cancer OR tumor)))'

# Output the number of search results for each option and keyword
$ gencube seqmeta --organism human --strategy chip,chip_seq --exclude cell_line,crispr liver,lung cancer,tumor --detail

# Save the integrated metadata
$ gencube seqmeta --organism human --strategy chip,chip_seq --exclude cell_line,crispr liver,lung cancer,tumor --metadata
```
---

#### Credits
This package was created with [`Cookiecutter`](https://github.com/audreyr/cookiecutter) and the [`audreyr/cookiecutter-pypackage`](https://github.com/audreyr/cookiecutter-pypackage) project template.
