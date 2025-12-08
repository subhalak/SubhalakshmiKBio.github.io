# Bulk-RNAseq analysis - From Raw reads to Count matrix

RNA sequencing (RNA-Seq) is a powerful and widely used technique for studying the transcriptome, the complete set of RNA transcripts produced by the genome under specific conditions or in a specific cell type. 

Unlike traditional microarrays, RNA-Seq provides a high-resolution, quantitative, and unbiased view of gene expression, enabling researchers to identify differentially expressed genes, discover novel transcripts, detect alternative splicing events, and quantify gene expression across the entire genome.

The RNA-Seq data analysis process can be broadly divided into two main stages: 

* The transformation of raw sequencing reads into a count matrix, and the 
* Downstream statistical and functional interpretation of gene expression data.

Altogether, the RNA-Seq analysis workflow is a multi-step process that transforms raw sequencing data into biologically meaningful insights, enabling researchers to understand gene regulation, identify biomarkers, and study disease mechanisms at the transcriptomic level.

# A. The transformation of raw sequencing reads into a count matrix:

1. It begins with quality control of raw sequencing reads in FASTQ format. 

2. Tools like FastQC are commonly used to evaluate sequencing quality, identify adapter contamination, and detect technical biases. If necessary, low-quality reads and adapter sequences are trimmed using programs such as Trimmomatic.

3. The cleaned reads are then aligned to a reference genome or transcriptome using aligners like HISAT2.

4. Following alignment, gene- or transcript-level quantification is performed with tools like featureCounts, resulting in a count matrix that represents the number of reads mapped to each gene across all samples.

# About the Dataset:

This repository provides a complete, step-by-step hands-on/workshop for performing RNA-seq data analysis using publicly available datasets from the Gene Expression Omnibus (GEO). The dataset used in this project, GSE106305, contains

The aim of this tutorial is to demonstrate how to process, analyze, and interpret RNA-seq data in a reproducible manner, covering the entire workflow from raw FASTQ files to the identification of differentially expressed genes and downstream biological insights.

The dataset for this study includes two groups: 

SRR



# B.Downstream statistical and functional interpretation of gene expression data:

1. The downstream analysis of this count data to derive meaningful biological insights. 

2. Differential gene expression (DGE) analysis is typically performed using statistical frameworks such as DESeq2 to identify genes that show significant changes in expression between conditions or experimental groups. 

3. The results are often visualized using principal component analysis (PCA), heatmaps, and volcano plots to explore sample relationships and highlight key regulatory genes. 

4. To interpret the biological relevance of differentially expressed genes, functional enrichment analysis is carried out using gene ontology (GO) and pathway analysis tools such as clusterProfiler or GSEA. 

5. This integrative workflow allows researchers to move from raw sequencing data to a deeper understanding of cellular states, molecular functions, and regulatory mechanisms underlying the biological system under study.

# Differential expression analysis using DeSeq2 in R:

DESeq2 uses a statistical model based on the negative binomial distribution to identify genes that show significant differences in expression between experimental conditions
(for example, Control vs. Disease).

In this workshop, we'll perform Differential Expression Analysis (DEA) using the DESeq2 package in R. The goal is to identify genes that are significantly upregulated or downregulated between the two conditions (Normaxia vs Hypoxia) in 2 different celllines (LNCaP vs PC3) using RNA-Seq count data generated from preprocessing the sequence reads from GSE106305.

All analyses in this step have been conducted in RStudio. Make sure you have both R and RStudio installed. You can download and install them by following the instructions available here: https://posit.co/download/rstudio-desktop/

# Step 1: Install Required Packages in R

Before any analysis, install all the necessary packages that will be used throughout the workflow. R’s Bioconductor repository hosts specialized bioinformatics tools. DESeq2 and other packages are maintained there to ensure reproducibility and compatibility.


