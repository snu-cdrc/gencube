"""
DIC_ORGANISM = {
    'human': 'homo sapiens',
    'mouse': 'mus musculus',
    'dog': 'canis lupus familiaris',
    'dingo': 'canis lupus dingo',
    'wolf': 'canis lupus',
    # Pet animals
    'cat': 'felis catus',
    # Farm animals
    'pig': 'sus scrofa',
    'pig_domestic': 'sus scrofa domesticus',
    'cow': 'bos taurus',
    'dairy_cow': 'Bos indicus',
    'chicken': 'gallus gallus',
    'horse': 'equus caballus',
    # Farm plants
    'rice': 'oryza sativa',
    'wheat': 'triticum aestivum',
    # Peto's paradox (Long-lived & cancer-free)
    'elephant': 'loxodonta africana', # 아프리카 코끼리
    'whale': 'balaenoptera musculus', # 흰수염고래
    'naked_mole_rat': 'heterocephalus glaber',
    'blind_mole_rat': 'spalax ehrenbergi',
    # Primates
    'gorilla': 'gorilla gorilla',
    'rhesus_monkey': 'macaca mulatta', # 리서스 원숭이
    'cynomolgus_monkey': 'macaca fascicularis', # 시노몰거스 원숭이
    'baboon': 'papio', # 바부인
    'chimpanzee': 'pan troglodytes', # 침팬지
    'marmoset': 'callithrix jacchus',  # 마모셋
    'macaque': 'macaca',
    'capuchin_monkey': 'cebus capucinus',  # Used in behavioral studies and neuroscience
    'squirrel_monkey': 'saimiri sciureus',  # Used in neurobiology, behavioral biology, and pharmacology
    'bonobo': 'pan paniscus',  # Close relation to chimpanzees, used in behavioral and social studies, genetics
    # Experimental models
    'yeast': 'saccharomyces cerevisiae',
    'fruit_fly': 'drosophila melanogaster',
    'nematode': 'caenorhabditis elegans',
    'zebrafish': 'danio rerio',
    'african clawed frog': 'xenopus laevis',
    'rat': 'rattus norvegicus',
    'guinea pig': 'cavia porcellus',
    'rabbit': 'oryctolagus cuniculus',
    # Etc
    'opossum': 'didelphis virginiana'
    }   

#######################################################
# NCBI data 
#######################################################
## 1. Fasta format ------------------------------------
v       _genomic.fna.gz
        _rna.fna.gz
        _rna_from_genomic.fna.gz        
x       _cds_from_genomic.fna.gz
        _pseudo_without_product.fna.gz

        _translated_cds.faa.gz
?       _protein.faa.gz

## 2. Gtf, Gff format ---------------------------------
v       _genomic.gff.gz
v       _genomic.gtf.gz

    Gnomon_models directory 
v        _gnomon_model.gff.gz
    
    Evidence_alignments directory
v        _cross_species_tx_alns.gff.gz
v        _same_species_tx_alns.gff.gz

## 3. RepeatMasker ------------------------------------
v        _rm.out.gz
v        _rm.run

## 4. Others ------------------------------------------
v       _genomic_gaps.txt.gz

        _feature_count.txt.gz
        _feature_table.txt.gz
        _assembly_report.txt
        _assembly_stats.txt
        assembly_status.txt
                
        annotation_hashes.txt
        _gene_ontology.gaf.gz

        _genomic.gbff.gz
        _rna.gbff.gz
        _protein.gpff.gz

    Gnomon_models directory 
            _gnomon_protein.faa.gz
            _gnomon_rna.fna.gz     


File formats and content:
## 1. Fasta format ------------------------------------
*_genomic.fna.gz file
    FASTA format of the genomic sequence(s) in the assembly. Repetitive 
    sequences in eukaryotes are masked to lower-case (see below).
    The FASTA title is formatted as sequence accession.version plus 
    description. The genomic.fna.gz file includes all top-level sequences in
    the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds,
    unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds
    that are part of the chromosomes are not included because they are
    redundant with the chromosome sequences; sequences for these placed 
    scaffolds are provided under the assembly_structure directory.
    
    >NC_006583.4 Canis lupus familiaris isolate Tasha breed boxer chromosome 1, alternate assembly Dog10K_Boxer_Tasha, whole genome shotgun sequence
    ttttttttttttttttttttctcgtttttttttttttttttctttttttctttcttttttttttccgttcTGTGAGCCAC
    TGGGCGGACCTCCCGTGTATCTGAGTGATACCGCGATTTCCAACGCTTTCATCAATCAGAGTAGCCCGTTCCTACCCATG
    TCCGAATGGGTAGCAGACATttgcctcccatgcatacacttctgagagttgagcttctGGCCTGTACCTCCCATCCCGCC
    GCTGGTACCCTTGGCTTCtcaaaggcctaggctcgctgtttGACGGAAGTGCTTGGAGAGAGAATTGGAATcaccggcac
    
*_cds_from_genomic.fna.gz
    FASTA format of the nucleotide sequences corresponding to all CDS 
    features annotated on the assembly, based on the genome sequence. See 
    the "Description of files" section below for details of the file format.
    
    >lcl|NC_006583.4_cds_XP_003638833.2_1 [gene=ENPP1] [db_xref=GeneID:100856150] [protein=ectonucleotide pyrophosphatase/phosphodiesterase family member 1] [protein_id=XP_003638833.2] [location=complement(join(1248058..1248228,1250132..1250294,1251863..1251995,1252828..1252908,1257503..1257632,1259193..1259347,1261248..1261299,1262129..1262298,1263921..1264008,1264736..1264805,1266156..1266283,1267437..1267468,1269237..1269368,1271589..1271697,1273918..1273990,1274301..1274366,1275268..1275377,1276443..1276562,1277567..1277646,1279545..1279642,1281587..1281647,1282525..1282650,1283470..1283586,1285357..1285429,1317509..1317592))] [gbkey=CDS]
    ATGGACCTGGGGGAGGAGCCGCTGGAGAAGGCGGCCCGCGCCCGCCCGGCCAAGGACCCCAACACCTACAAGGTGCTCTC
    GCTGGTATTGTCAGTCTGTGTGTTAACAACAATTCTTGGTTGTATATTTGGTTTGAAACCAAGCTGTGCCAAAGAAGTTA

*_rna.fna.gz file
    FASTA format of accessioned RNA products annotated on the genome 
    assembly; Provided for RefSeq assemblies as relevant (Note, RNA and mRNA 
    products are not instantiated as a separate accessioned record in GenBank
    but are provided for some RefSeq genomes, most notably the eukaryotes.)
    The FASTA title is provided as sequence accession.version plus 
    description.

    >NM_001002930.1 Canis lupus familiaris neuropeptide Y receptor Y1 (NPY1R), mRNA
    ATGAATTCAACATCATTTTCCCAGGTTGAAAACCATTCAATCTTCTGTAATTTTTCAGAGAATTCCCAGTTTTTGGCTTT
    TGAAAGTGATGATTGTCACCTGCCCTTGGCCATGATATTTACATTAGCTCTTGCTTACGGAGCTGTAATAATTCTTGGGG
    TCACTGGAAACCTGGCCTTGATCATGATCATCTTGAAACAAAAGGAGATGAGAAATGTTACCAATATCCTGATTGTGAAC
    CTTTCCTTCTCAGACTTGCTTGTTGCCATCATGTGTCTTCCCTTCACATTTGTCTACACCTTAATGGACCACTGGGTTTT
        
*_rna_from_genomic.fna.gz
    FASTA format of the nucleotide sequences corresponding to all RNA 
    features annotated on the assembly, based on the genome sequence. See 
    the "Description of files" section below for details of the file format.
    
    >lcl|NC_006583.4_ncrna_XR_005316961.1_1 [gene=LOC119870367] [db_xref=GeneID:119870367] [product=uncharacterized LOC119870367] [ncRNA_class=lncRNA] [transcript_id=XR_005316961.1] [location=join(651385..652415,695111..695247,696316..696387,696906..698305,788676..789323)] [gbkey=ncRNA]
    CGCCTGAGCGTCCCGCAGGGGGCACCTGAGCCCAGGCCCAACCCAGCAGGCATATATGTGCAGTAGGCCTCCTCTGCCTG
    CACTATCAGGAGCCAGCCTGTGCTGGGGAGGCCTGAGGAGGACGTGGGGACGTGTGCACCCTGATGAGCCACAGCTTCAG
    GAAATTGTGGGCTCGCATCCACACTGAAAGGAGCCCCTGGGCCTCCTCAGGACCCACGGGGAGGCTTCAGGAAGATGTGG
    GCCTGTATCCACGCTGACAGGAGCCCCTGGGCCACCTCAGGACATAAGTCGAGGCTTCTGGAAGATGGGGGCCGCATTCA
        
*_pseudo_without_product.fna.gz
    FASTA format of the genomic sequence corresponding to pseudogene and 
    other gene regions which do not have any associated transcribed RNA 
    products or translated protein products. It includes annotated gene 
    regions that require rearrangement to provide the final product, e.g.
    immunoglobulin segments. These sequences are not assigned accession 
    numbers, and are derived directly from the assembled genomic sequences.
    The FASTA title has a local sequence identifier, the Gene ID and gene 
    name.

    >lcl|NC_006583.4:LOC100684726 GeneID:100684726 (LOC100684726) on NC_006583.4 [Canis lupus familiaris]
    GAGAGATTCATTCGGCAGCTGCTGGAAGACACGGTCAAGGTAGCGGCAGTGTCTGGCCTGATGCGGAGACCCCCGGAACA
    GGTGTCGGGATTGCTGAGGAGGCGCTTTCACCGGACGGCGCTGGCGGCGCTGCAGGTGACAGTTCGTGAGGCTCTAAACC
    AAGGTACGGATGAAGAGCTGGAAAGAGATGAGAAGGTATTTATGCTTGGGGAAGAAGTTGCCCAGTATGATGGTGCATAT
    AAGGTTAGTCGAGGCCTGTGGAAGAAATATGGCGATAAGAGAATCATAGATACTCCCATATCTGAGATGGGCTTTGCTGG

*_translated_cds.faa.gz
    FASTA sequences of individual CDS features annotated on the genomic 
    records, conceptually translated into protein sequence. The sequence 
    corresponds to the translation of the nucleotide sequence provided in the
    *_cds_from_genomic.fna.gz file.
    
    >lcl|NC_006583.4_prot_XP_003638833.2_1 [gene=ENPP1] [db_xref=GeneID:100856150] [protein=ectonucleotide pyrophosphatase/phosphodiesterase family member 1] [protein_id=XP_003638833.2] [location=complement(join(1248058..1248228,1250132..1250294,1251863..1251995,1252828..1252908,1257503..1257632,1259193..1259347,1261248..1261299,1262129..1262298,1263921..1264008,1264736..1264805,1266156..1266283,1267437..1267468,1269237..1269368,1271589..1271697,1273918..1273990,1274301..1274366,1275268..1275377,1276443..1276562,1277567..1277646,1279545..1279642,1281587..1281647,1282525..1282650,1283470..1283586,1285357..1285429,1317509..1317592))] [gbkey=CDS]
    MDLGEEPLEKAARARPAKDPNTYKVLSLVLSVCVLTTILGCIFGLKPSCAKEVKSCKGRCFERTFGNCRCDVACVDLGNC
    CLDYQETCIEPERIWTCSKFRCGEKRLSRSLCSCSDDCRDKGDCCVNYSSVCLGEKSWVEETCESIDEPQCPAGFEMPPT
    LLFSLDGFRAEYLHTWGGLLPVISKLKNCGTYTKNMRPVYPTKTFPNHYSIVTGLYPESHGIIDNKIYDPKMNAFFALKS
    KEKFNPEWYKGEPIWLTTKYQGLKSGTFFWPGSDVEIKGILPDIYKMYNGSIPFEERILAVLKWLQLPKDERPHFYTLYL

*_protein.faa.gz file
    FASTA format sequences of the accessioned protein products annotated on
    the genome assembly. The FASTA title is formatted as sequence 
    accession.version plus description.

Gnomon_models directory 
*_gnomon_rna.fna.gz
    FASTA format sequences of Gnomon transcript models annotated on the 
    genome assembly. The FASTA title is the Gnomon identifier for the 
    transcript (>gnl|GNOMON|XXX.m)
 
*_gnomon_protein.faa.gz
    FASTA format sequences of Gnomon protein models annotated on the genome
    assembly. The FASTA title is the Gnomon identifier for the protein model
    (>gnl|GNOMON|XXX.p) 

    >NP_001002930.1 neuropeptide Y receptor type 1 [Canis lupus familiaris]
    MNSTSFSQVENHSIFCNFSENSQFLAFESDDCHLPLAMIFTLALAYGAVIILGVTGNLALIMIILKQKEMRNVTNILIVN
    LSFSDLLVAIMCLPFTFVYTLMDHWVFGEAMCKLNPFVQCVSITVSIFSLVLIAVERHQLIINPRGWRPNNRHAYVGIAV
    IWVLAVVSSLPFLIYQVLTDEPFQNVTLDAFKDKYVCFDKFPSDSHRLSYTTLLLMLQYFGPLCFIFICYFKIYIRLKRR
    NNMMDKMRDNKYRSSETKRINIMLLSIVVAFAVCWLPLTIFNTVFDWNHQIIATCNHNLLFLLCHLTAMISTCVNPIFYG
    

## 2. Gtf, Gff format ---------------------------------

*_genomic.gff.gz file
    Annotation of the genomic sequence(s) in Generic Feature Format Version 3
    (GFF3). Sequence identifiers are provided as accession.version.
    Additional information about NCBI's GFF files is available at 
    https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.

    ##gff-version 3
    #!gff-spec-version 1.21
    #!processor NCBI annotwriter
    #!genome-build Dog10K_Boxer_Tasha
    #!genome-build-accession NCBI_Assembly:GCF_000002285.5
    
*_genomic.gtf.gz file
    Annotation of the genomic sequence(s) in Gene Transfer Format Version 2.2
    (GTF2.2). Sequence identifiers are provided as accession.version.

    #gtf-version 2.2
    #!genome-build Dog10K_Boxer_Tasha
    #!genome-build-accession NCBI_Assembly:GCF_000002285.5
    #!annotation-source NCBI Canis lupus familiaris Annotation Release 106
    NC_006583.4     Gnomon  gene    29337   30631   .       +       .       gene_id "LOC100684726"; transcript_id ""; db_xref "GeneID:100684726"; gbkey "Gene"; gene "LOC100684726"; gene_biotype "pseudogene"; pseudo "true"; 

Gnomon_models directory 
*_gnomon_model.gff.gz
    Gnomon annotation of the genomic sequence(s) in Generic Feature Format
    Version 3 (GFF3). Sequence identifiers are provided as accession.version
    for the genomic sequences and Gnomon identifiers for the Gnomon models:
    gene.XXX for genes, GNOMON.XXX.m for transcripts and GNOMON.XXX.p for 
    proteins. These identifiers are NOT universally unique. They are unique
    per annotation release only. Additional information about NCBI's GFF 
    files is available at 
    https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.
    
Evidence_alignments directory
*_cross_species_tx_alns.gff.gz
    Alignments of cDNAs, ESTs and TSAs from other species to the genomic
    sequence(s) in Generic Feature Format Version 3 (GFF3) [not all 
    annotation releases have cross-species alignments]. These alignments may
    have been used as evidence for gene prediction by the annotation 
    pipeline. Sequence identifiers are provided as accession.version. 
    Additional information about NCBI's GFF files is available at 
    https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.

*_same_species_tx_alns.gff.gz
    Alignments of same-species cDNAs, ESTs and TSAs to the genomic 
    sequence(s) in Generic Feature Format Version 3 (GFF3). These alignments
    were used as evidence for gene prediction by the annotation pipeline. 
    Sequence identifiers are provided as accession.version. Additional 
    information about NCBI's GFF files is available at 
    https://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.
    
## 3. RepeatMasker ------------------------------------

*_rm.out.gz file
    RepeatMasker output; 
    Provided for Eukaryotes
    
    SW     perc    perc    perc  query           position in query                   matching       repeat           position in repeat
    score     div.    del.    ins.  sequence        begin      end             (left)   repeat         class/family   begin  end    (left)     ID

    40   10.200   0.000   0.000  NC_006583.4              1         64 (122014004) + (T)n           Simple_repeat       1     64    (0)      1
    40   10.200   0.000   0.000  NC_006583.4              1         64 (122014004) + (T)n           Simple_repeat       1     64    (0)      1
    
*_rm.run file
    Documentation of the RepeatMasker version, parameters, and library; 
    Provided for Eukaryotes 
    
## 4. Others ------------------------------------------

*_genomic_gaps.txt.gz
    Tab-delimited text file reporting the coordinates of all gaps in the 
    top-level genomic sequences. The gaps reported include gaps specified in
    the AGP files, gaps annotated on the component sequences, and any other 
    run of 10 or more Ns in the sequences. See the "Description of files" 
    section below for details of the file format.

    # accession.version     start   stop    gap_length      gap_type        linkage_evidence
    NC_006583.4     13769   13793   25      within_scaffold paired-ends
    NC_006583.4     33049   33148   100     within_scaffold paired-ends
    NC_006583.4     58510   58609   100     within_scaffold paired-ends
    NC_006583.4     79254   79278   25      within_scaffold paired-ends

*_feature_count.txt.gz
    Tab-delimited text file reporting counts of gene, RNA, CDS, and similar
    features, based on data reported in the *_feature_table.txt.gz file.
    See the "Description of files" section below for details of the file 
    format.
    
    # Feature       Class   Full Assembly   Assembly-unit accession Assembly-unit name      Unique Ids      Placements
    CDS     with_protein    GCF_000002285.5         all     62347   62347
    CDS     without_protein GCF_000002285.5         all     na      90
    C_region                GCF_000002285.5         all     na      16
    V_segment               GCF_000002285.5         all     na      75
        
*_feature_table.txt.gz
    Tab-delimited text file reporting locations and attributes for a subset 
    of annotated features. Included feature types are: gene, CDS, RNA (all 
    types), operon, C/V/N/S_region, and V/D/J_segment. Replaces the .ptt & 
    .rnt format files that were provided in the old genomes FTP directories.
    See the "Description of files" section below for details of the file 
    format.
    
    # feature       class   assembly        assembly_unit   seq_type        chromosome      genomic_accession       start   end     strand  product_accession       non-redundant_refseq    related_accession       name    symbol GeneID   locus_tag       feature_interval_length product_length  attributes
    gene    pseudogene      GCF_000002285.5 Primary Assembly        chromosome      1       NC_006583.4     29337   30631   +                                       LOC100684726    100684726               1295            pseudo
    gene    lncRNA  GCF_000002285.5 Primary Assembly        chromosome      1       NC_006583.4     651385  789323  +                                       LOC119870367    119870367               137939
    ncRNA   lncRNA  GCF_000002285.5 Primary Assembly        chromosome      1       NC_006583.4     651385  789323  +       XR_005316961.1                  uncharacterized LOC119870367    LOC119870367    119870367              3288     3288
    gene    lncRNA  GCF_000002285.5 Primary Assembly        chromosome      1       NC_006583.4     722713  747463  +                                       LOC119870369    119870369               24751

*_assembly_report.txt file
    Tab-delimited text file reporting the name, role and sequence 
    accession.version for objects in the assembly. The file header contains 
    meta-data for the assembly including: assembly name, assembly 
    accession.version, scientific name of the organism and its taxonomy ID, 
    assembly submitter, and sequence release date.
    
    # Assembly name:  Dog10K_Boxer_Tasha
    # Organism name:  Canis lupus familiaris (dog)
    # Infraspecific name:  breed=boxer
    # Isolate:  Tasha
    # Sex:  female

*_assembly_stats.txt file
    Tab-delimited text file reporting statistics for the assembly including: 
    total length, ungapped length, contig & scaffold counts, contig-N50, 
    scaffold-L50, scaffold-N50, scaffold-N75, and scaffold-N90

    # Assembly Statistics Report
    # Assembly name:  Dog10K_Boxer_Tasha
    # Organism name:  Canis lupus familiaris (dog)
    # Infraspecific name:  breed=boxer
    # Isolate:  Tasha

assembly_status.txt
    A text file reporting the current status of the version of the assembly
    for which data is provided. Any assembly anomalies are also reported.

    status=latest


    
*_ani_contam_ranges.tsv file
    Tab-delimited text file reporting potentially contaminated regions in
    the assembly identified based on Average Nucleotide Identity (ANI)
    calculations against a reference set of assemblies from type materials
    and provides actionable recommendations for addressing identified 
    contamination. This file is provided only for prokaryotic assemblies
    and only when potential contamination is detected.
*_ani_report.txt file
    Tab-delimited text file reporting Average Nucleotide Identity (ANI)
    based evaluation of the taxonomic identity of the assembly against a
    reference set of assemblies from type materials.
    
*_assembly_regions.txt
    Provided for assemblies that include alternate or patch assembly units. 
    Tab-delimited text file reporting the location of genomic regions and 
    listing the alt/patch scaffolds placed within those regions.
    
*_assembly_structure directory
    This directory will only be present if the assembly has internal 
    structure. When present, it will contain AGP files that define how 
    component sequences are organized into scaffolds and/or chromosomes. 
    Other files define how scaffolds and chromosomes are organized into 
    non-nuclear and other assembly-units, and how any alternate or patch 
    scaffolds are placed relative to the chromosomes. Refer to the README.txt
    file in the assembly_structure directory for additional information.

*_normalized_gene_expression_counts.txt.gz
    Tab-delimited text file with normalized counts of RNA-seq reads mapped to 
    each gene. See "Description of files" section below for details of the 
    file format.
*_gene_expression_counts.txt.gz
    Tab-delimited text file with counts of RNA-seq reads mapped to each gene.
    See "Description of files" section below for details of the file format.
*_gene_ontology.gaf.gz
    Gene Ontology (GO) annotation of the annotated genes in GO Annotation 
    File (GAF) format. Additional information about the GAF format is 
    available at 
    http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/ 

    !gaf-version: 2.2
    !generated-by: NCBI
    !date-generated: 2023-11-11
    !
    !interProScan version: 5.64-96.0
    
*_genomic.gbff.gz file
    GenBank flat file format of the genomic sequence(s) in the assembly. This
    file includes both the genomic sequence and the CONTIG description (for 
    CON records), hence, it replaces both the .gbk & .gbs format files that 
    were provided in the old genomes FTP directories.
    
    LOCUS       NC_006583          122014068 bp    DNA     linear   CON 08-JAN-2021
    DEFINITION  Canis lupus familiaris isolate Tasha breed boxer chromosome 1,
                alternate assembly Dog10K_Boxer_Tasha, whole genome shotgun
                sequence.
    ACCESSION   NC_006583

*_rna.gbff.gz file
    GenBank flat file format of RNA products annotated on the genome 
    assembly; Provided for RefSeq assemblies as relevant

    LOCUS       NM_001002930            1149 bp    mRNA    linear   MAM 05-DEC-2020
    DEFINITION  Canis lupus familiaris neuropeptide Y receptor Y1 (NPY1R), mRNA.
    ACCESSION   NM_001002930
    VERSION     NM_001002930.1
    KEYWORDS    RefSeq.
    
*_protein.gpff.gz file
    GenPept format of the accessioned protein products annotated on the 
    genome assembly

    LOCUS       NP_001002930             382 aa            linear   MAM 05-DEC-2020
    DEFINITION  neuropeptide Y receptor type 1 [Canis lupus familiaris].
    ACCESSION   NP_001002930
    VERSION     NP_001002930.1
    DBSOURCE    REFSEQ: accession NM_001002930.1
    
annotation_hashes.txt
    Tab-delimited text file reporting hash values for different aspects
    of the annotation data. See the "Description of files" section below 
    for details of the file format.

    # Assembly accession    Descriptors hash        Descriptors last changed        Features hash   Features last changed   Locations hash  Locations last change   Protein names hash      Protein names last changed
    GCF_000002285.5 AD02EF5A517E29D662E47AC80F6AC064        2021/01/09 09:05:00     430AC2FC31A11C8625AC3994B706B6F7        2021/01/09 09:05:00     39A44E1A63D662505565159FAD66DDC2        2021/01/09 09:05:00     25084D25D718F848ADDCAF0A1927DE7F        2021/01/09 09:05:00


*_rnaseq_alignment_summary.txt
    Tab-delimited text file containing counts of alignments that were either
    assigned to a gene or skipped for a specific reason. See "Description of
    files" section below for details of the file format.
*_rnaseq_runs.txt
    Tab-delimited text file containing information about RNA-seq runs used 
    for gene expression analyses (See *_gene_expression_counts.txt file and 
    *.bw files within "RNASeq_coverage_graphs" directory). 
    
*_wgsmaster.gbff.gz
    GenBank flat file format of the WGS master for the assembly (present only
    if a WGS master record exists for the sequences in the assembly).
    
	   
Additional directories and files provided for organisms annotated by the NCBI 
Eukaryotic Genome Annotation Pipeline:

RefSeq_transcripts_alignments directory
*_knownrefseq_alns.bam
    Alignments of the annotated Known RefSeq transcripts (identified with 
    accessions prefixed with NM_ and NR_) to the genome in BAM format [not 
    all annotation releases have Known RefSeq transcripts]. For more 
    information about the BAM format see: 
    https://samtools.github.io/hts-specs/SAMv1.pdf
*_knownrefseq_alns.bam.bai
    Index of the BAM alignments of the annotated Known RefSeq transcripts 
    to the genome. [not all annotation releases have Known RefSeq 
    transcripts]
*_modelrefseq_alns.bam
    Alignments of the annotated Model RefSeq transcripts (identified with 
    accessions prefixed with XM_ and XR_) to the genome in BAM format. For 
    more information about the BAM format see:
    https://samtools.github.io/hts-specs/SAMv1.pdf
*_modelrefseq_alns.bam.bai
    Index of the BAM alignments of the annotated Model RefSeq transcripts to
    the genome.

Annotation_comparison directory
    This directory is only provided for re-annotations of the same species.
*_compare_prev.txt.gz
    Matching genes and transcripts in the current and previous annotation 
    releases binned by type of difference (column 1 for genes and column 14 
    for transcripts), in tabular format.

RNASeq_coverage_graphs directory
*_graph.bw
    RNA-seq read coverage graphs in UCSC BigWig file format
    (https://genome.ucsc.edu/goldenPath/help/bigWig.html)
    
md5checksums.txt file
    file checksums are provided for all data files in the directory


#######################################################
# GENEARK data 
#######################################################
Basic: 
“GC Percent,” “CpG Islands,” a “Simple Repeats” track (Tandem Repeats Finder)
RepeatMasker track

Gene models: 
AUGUSTUS de novo predictor
“Xeno RefGene”
“NCBI RefSeq”
Ensembl Rapid Release
TOGA (Tool to infer Orthologs from Genome Alignments) gene models.



 
"""