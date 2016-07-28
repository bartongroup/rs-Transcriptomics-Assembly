# Transcriptomics-Assembly
De-novo Transcriptomics Assembly workflow for four Dictyostelium species (e.g.- Dictyostelium discoideum, Polysphondylium pallidum, Dictyostelium Lacteum and Dictyostelium Fasciculatum). This is the standard assembly workflow that should ideally work on any organism.

### Read Normalization
### Read Assembly using Trinity
    alignReads.pl --left reads.ALL.left.fq.normalized_K25_C50_pctSD200.fq --right reads.ALL.right.fq.normalized_K25_C50_pctSD200.fq --seqType fq --target NewAssembly_35631.fasta.clean --aligner bowtie --retain_intermediate_files
### Statistics
### PASA Assembly
### PASA updation with the existing Annotation
### Quality Control Measurement
#### Annotation
#### Blobplot
#### CEGMA
#### Transrate

