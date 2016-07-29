# Transcriptomics-Assembly
De-novo Transcriptomics Assembly workflow for four Dictyostelium species (e.g.- Dictyostelium discoideum, Polysphondylium pallidum, Dictyostelium Lacteum and Dictyostelium Fasciculatum). This is the standard assembly workflow that should ideally work on any organism.Before normalizing first concatenate all RNAseq data across all samples into a single set of inputs to generate a single reference transcriptomics assembly.  Combine all left reads in one file and all right reads in another file. In order to reduce the number, raw reads were normalized using in silico digital normalization implemented in trinity at 50X coverage. The reads were assembled with Trinity using kmer parameter of 25.

### Read Normalization
    trinity/util/normalize_by_kmer_coverage.pl --seqType fq --JM 10G --max_cov 75 --left All_Left_Reads.fq --right All_Right_Reads.fq --pairs_together --PARALLEL_STATS --JELLY_CPU 6
### Read Assembly using Trinity
    alignReads.pl --left reads.ALL.left.fq.normalized_K25_C50_pctSD200.fq --right reads.ALL.right.fq.normalized_K25_C50_pctSD200.fq --seqType fq --target NewAssembly_35631.fasta.clean --aligner bowtie --retain_intermediate_files
### Assembly Statistics
    SAM_nameSorted_to_uniq_count_stats.pl bowtie_out/bowtie_out.nameSorted.sam >AlignmentStatistics
### PASA Assembly
    PASA/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g PN500.fa -t NewAssembly_35631.fasta.clean --TDN tdn.accs --TRANSDECODER --ALT_SPLICE --ALIGNERS blat,gmap
### PASA updation with the existing Annotation
    PASA/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g PN500.fa -P P_Pal_2.gff
    PASA/scripts/Launch_PASA_pipeline.pl -c annotCompare.config  -A -g PN500.fa -t NewAssembly_35631.fasta.clean
### Quality Control Measurement
#### Annotation
#### Blobplot
#### CEGMA
#### Transrate

