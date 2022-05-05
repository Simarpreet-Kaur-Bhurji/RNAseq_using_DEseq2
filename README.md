# RNAseq analysis with DEseq2 and functional enrichment analysis

Interested in exploring more applications of the RNASeq, Please read here more https://doi.org/10.1093/bib/bbab259

# About the RNA-Seq analysis
The R script performs several steps in RNAseq gene differential expression analysis, including filtering, preprocessing, visualization, clustering, and Enrichment. For the analysis, several R Bioconductor packages are required to be installed (Installation commands are provided in the script. However, users can also refer to the Bioconductor website for detailed instructions). 

# Required data files
You should have a raw count and annotation/metadata file for running this analysis. Raw count files are usually obtained from tools such as featureCount, Rsem etc.

# Bioconductor packages to be installed
 DESeq2
 
 edgeR
 
 biomaRt (Very useful for gene filtering and annotations)
 
 PCAtools (PCA detailed analysis)
 
 ReactomePA (enrichment analysis)
