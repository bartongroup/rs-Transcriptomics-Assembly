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
##### In case the existing genome annotation in Augustus file format
    grep "AUGUSTUS" PN500_augustus_prediction.gff | awk -F "\t" '$3 ~/gene|transcript|exon|CDS/' >P_Pal.gff 
    python Format_Gff.py >P_Pal_1.gff
    python Format_Gff2.py >P_Pal_2.gff
    PASA/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g PN500.fa -P P_Pal_2.gff
    PASA/scripts/Launch_PASA_pipeline.pl -c annotCompare.config  -A -g PN500.fa -t NewAssembly_35631.fasta.clean
### Quality Control Measurement
#### Annotation
    blat PN500.fa ../NewAssembly_35631.fasta.clean -t=dna -q=dna -out=psl Trinity-Genome.psl
    perl blat2gbrowse.pl Trinity-Genome.psl Reema1.gff
    less Reema1.gff | grep "HSP" >NewAssembly.gtf
    perl GTF1.pl NewAssembly.gtf >NewAssembly_Final.gtf
#### Blobplot
    
#### CEGMA
#### Transrate
    transrate --assembly NewAssembly_35631.fasta.clean --left reads.ALL.left.fq.normalized_K25_C50_pctSD200.fq --right reads.ALL.right.fq.normalized_K25_C50_pctSD200.fq --reference PN500_augustus_prediction_test.aa --outfile Reference_Based
    transrate --assembly DNA.fas --left reads.ALL.left.fq.normalized_K25_C50_pctSD200.fq --right reads.ALL.right.fq.normalized_K25_C50_pctSD200.fq --reference PN500_augustus_prediction_test.aa --outfile Reference_Based
    transrate --assembly pasa_Pallidum_Gernot.assemblies.fasta --left reads.ALL.left.fq.normalized_K25_C50_pctSD200.fq --right reads.ALL.right.fq.normalized_K25_C50_pctSD200.fq --reference PN500_augustus_prediction_test.aa --outfile Reference_Based
    transrate --assembly PASA2.fasta --left reads.ALL.left.fq.normalized_K25_C50_pctSD200.fq --right reads.ALL.right.fq.normalized_K25_C50_pctSD200.fq --reference PN500_augustus_prediction_test.aa --outfile Reference_Based

### Plots
#### Assembly Depth Plots
    dicty <- read.table(file="NewAssembly_31259_LengthnCount",sep="\t")
    pal <- read.table(file="Pallidum_Length_Count",sep="\t")
    fas <- read.table(file="Dfas_Length_Count",sep="\t")
    lac <- read.table(file="Dlac_Length_Count",sep="\t")

    t_frame <- data.frame(Species="D.discoideum",Count=dicty$V3)
    p1_frame <- data.frame(Species="P.pallidum",Count=pal$V3)
    p2_frame <- data.frame(Species="D.fasciculatum",Count=fas$V3)
    d_frame <- data.frame(Species="D.lacteum",Count=lac$V3)
    all_frame <- rbind(t_frame,p1_frame,p2_frame,d_frame)
    
    library(ggplot2)
    pdf("All_depth.pdf")
    ggplot(all_frame, aes(x=Species, y=Count,fill=Species)) +geom_boxplot()+ scale_y_continuous(trans=log10_trans())+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())
    dev.off()

#### CEGMA Plots

    a <- read.table(file="Cegma-Dicty1.txt",header=TRUE)
    b <- read.table(file="cegma-Pal.txt",header=TRUE)
    c <- read.table(file="cegma-Fasc.txt",header=TRUE)
    d <- read.table(file="cegma-Lactum.txt",header=TRUE)
    
    require(ggplot2)
    theme_set(theme_bw(14))
    p <- ggplot(data=a,aes(x=Name,y=Score,fill=Type))+geom_bar(stat="identity")+xlab("")+ggtitle("Dictyostelium discoideum")+ylab("Number of CEGs")+theme(text = element_text(size=20))+scale_fill_manual(values = c("skyblue3","skyblue2"))
    p1 <- ggplot(data=b,aes(x=Name,y=Score,fill=Type))+geom_bar(stat="identity")+xlab("")+ggtitle("Polysphondylium pallidum")+ylab("Number of CEGs")+theme(text = element_text(size=20))+ scale_fill_manual(values = c("burlywood3","burlywood2"))
    p2 <- ggplot(data=c,aes(x=Name,y=Score,fill=Type))+geom_bar(stat="identity")+xlab("")+ggtitle("Dictyostelium fasciculatum")+ylab("Number of CEGs")+theme(text = element_text(size=20))+ scale_fill_manual(values = c("indianred3","indianred2"))
    p3 <- ggplot(data=d,aes(x=Name,y=Score,fill=Type))+geom_bar(stat="identity")+xlab("")+ggtitle("Dictyostelium lacteum")+ylab("Number of CEGs")+theme(text = element_text(size=20))+scale_fill_manual(values = c("palegreen3","palegreen2"))
    
    source("multiplot.R")
    pdf("Cegma_staggerd.pdf",width=15,height=15)
    multiplot(p,p2,p1,p3)
    dev.off()

#### Transrate
##### Assembly Score
    pal <- read.table(file="QA_Ref_Coveragre.txt",sep="\t",header=TRUE)
    dp1_frame <- data.frame(name="Ddis_PASA1",n=pal$X,N=pal$DPasa1)
    pp1_frame <- data.frame(name="Ppal_PASA1",n=pal$X,N=pal$PPasa1)
    Fp1_frame <- data.frame(name="DFas_PASA1",n=pal$X,N=pal$FPasa1)
    lp1_frame <- data.frame(name="DLac_PASA1",n=pal$X,N=pal$LPasa1)
    dp2_frame <- data.frame(name="Ddis_PASA2",n=pal$X,N=pal$DPasa2)
    lp2_frame <- data.frame(name="DLac_PASA2",n=pal$X,N=pal$LPasa2)
    pp2_frame <- data.frame(name="Ppal_PASA2",n=pal$X,N=pal$PPasa2)
    Fp2_frame <- data.frame(name="DFas_PASA2",n=pal$X,N=pal$FPasa2)
    
    all_frame <- rbind(dp1_frame,pp1_frame,Fp1_frame,lp1_frame,dp2_frame,pp2_frame,Fp2_frame,lp2_frame)
    library(ggplot2)
    library(scales)
    p <- ggplot(all_frame,aes(x=n,y=N,group=name,shape=name))+geom_point(aes(colour=name),size=4)+geom_line(aes(colour=name, group=name))+scale_color_manual(values=c("#FF3333","#0000CC","#00FF00","#FF00CC","#FF3333","#0000CC","#00FF00","#FF00CC"))+scale_shape_manual(values = c(1,1,1,1,16,16,16,16))+expand_limits(y=0)+xlab("Reference Coverage") + ylab("Proportion of reference transcripts")+ theme_bw() +theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),text = element_text (size=20), legend.position ="top")+geom_hline(yintercept=c(0.22,0.35),linetype="dotted",color=c("#006633","#66CCFF"))
    
    pal <- read.table(file="QA_Score.txt",sep="\t",header=TRUE)
    dp1_frame <- data.frame(name="Ddis_PASA1",n=pal$X,N=pal$DPasa1)
    pp1_frame <- data.frame(name="Ppal_PASA1",n=pal$X,N=pal$PPasa1)
    Fp1_frame <- data.frame(name="DFas_PASA1",n=pal$X,N=pal$FPasa1)
    lp1_frame <- data.frame(name="DLac_PASA1",n=pal$X,N=pal$LPasa1)
    dp2_frame <- data.frame(name="Ddis_PASA2",n=pal$X,N=pal$DPasa2)
    lp2_frame <- data.frame(name="DLac_PASA2",n=pal$X,N=pal$LPasa2)
    pp2_frame <- data.frame(name="Ppal_PASA2",n=pal$X,N=pal$PPasa2)
    Fp2_frame <- data.frame(name="DFas_PASA2",n=pal$X,N=pal$FPasa2)
    
    all_frame <- rbind(dp1_frame,pp1_frame,Fp1_frame,lp1_frame,dp2_frame,pp2_frame,Fp2_frame,lp2_frame)
    p1 <- ggplot(all_frame,aes(x=n,y=N,group=name,shape=name))+geom_point(aes(colour=name),size=4)+geom_line(aes(colour=name, group=name))+scale_color_manual(values= c("#FF3333","#0000CC","#00FF00","#FF00CC","#FF3333","#0000CC","#00FF00","#FF00CC"))+scale_shape_manual(values = c(1,1,1,1,16,16,16,16)) +expand_limits(y=0)+ xlab(" ") +ylab("Score")+theme_bw() +theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),text = element_text(size=20))+geom_hline(yintercept=c(0.22,0.35),linetype="dotted",color=c("#006633","#66CCFF"))
    
    source("../multiplot.R")
    pdf("Assembly_Score.pdf",width=15,height=15)
    multiplot(p,p1)
    dev.off()
    
##### Contig Score

    

