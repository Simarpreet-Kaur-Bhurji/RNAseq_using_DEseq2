### RNAseq using Deseq2 and Functional enrichment Analysis ####
### Dr. Amarinder Singh Thind and Simarpreet Kaur
### Date : 18-19 April, 2022

##### Install packages, if not done before 

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("biomaRt")
# BiocManager::install('PCAtools')
# BiocManager::install('EnhancedVolcano')

###################### load the raw count matrix #######################

#setwd("/Users/athind/Dropbox/RNAseq_using_DEseq2-april16/") #Path_to_working_directory

rawcount<-read.table ("RawCount_input.csv",header=TRUE,  sep=",",  row.names=1)

## Replace NAs by zero and 
rawcount <- round(rawcount) 
rawcount[is.na(rawcount)] <- 0


###################### Data annotation  #################################

anno <-read.table ("Annotation_of_samples_12_Samples_ALL.csv",header=TRUE,  sep=",", row.names = 1) ##In this case we have 3 coulmns (a) sample (b) Condition (c) batch
#rownames(anno) <- anno$sample  ##add rownames as sample name (if not already), because pca function check rownames of anno == col of data matrix

##############################################################
############ PCA plot for pre DE investigation ##############


library(PCAtools)

anno <- anno[match(colnames(rawcount), anno$Sample),] ## reordering anno row with colnmaes of rawcount
lograwcount <- as.matrix(log2(rawcount +1))  ## log transformation of rawcount for PCA plot 

top1000.order <- head(order(matrixStats::rowVars(lograwcount), decreasing = TRUE), 1000)
p <- PCAtools::pca(mat = lograwcount[top1000.order,], metadata = anno, removeVar = 0.01)

biplot(p,lab = paste0(p$metadata$Sample),
        colby = 'Batch',  #Sample #Batch #Condition #sex
        hline = 0, vline = 0,
        legendPosition = 'right',
        encircle = T )


##############################################################
################# Lets check combat normalization ############
############## SVA #####################
  
#BiocManager::install("sva")
  
library('sva')
rawcount <- as.matrix(rawcount)
adjusted_counts <- ComBat_seq(rawcount, batch=anno$Batch, group=anno$Condition) ##In ComBat-seq, user may specify biological covariates, whose signals will be preserved in the adjusted data. I
  
nor_set <- as.matrix(log2(adjusted_counts+1)) ## log transformation of adjusted count
top1000.order <- head(order(matrixStats::rowVars(nor_set), decreasing = TRUE), 1000)
pp <- PCAtools::pca(mat =nor_set[top1000.order,] , metadata = anno, removeVar = 0.01)
biplot(pp,
      lab = paste0(p$metadata$Sample),
      #colby = 'Batch',   #Batch_log', #Condition
      colby = 'Condition',
      hline = 0, vline = 0,
      legendPosition = 'right',encircle = T)

   

  ##### Do we suppose to remove any defaulty sample  #########
   
   ### subset raw and conditional data for defined pairs
   ##### Removing sample number 7 ########## 
   
anno <- anno[!(anno$Sample == 'sample_7' | anno$Sample == 'sample_8'),]
rawcount <- as.data.frame(rawcount)
rawcount <- rawcount[,names(rawcount) %in% anno$Sample]
    
### Go back to PCA plot and check what happned 
### perform combat normalization again after removal of sample
rawcount <- as.matrix(rawcount)
adjusted_counts <- ComBat_seq(rawcount, batch=anno$Batch, group=anno$Condition) ##In ComBat-seq, user may specify biological covariates, whose signals will be preserved in the adjusted data. I
  

# Define conditions (for contrast) that you want to compare if you have more than one #control #case
# This is pair-wise comparison, so only consider one pair at one time

firstC<-"Condition_A"       #case1 #case2 #case3 etc          
SecondC <-"Condition_B"     
p.threshold <- 0.05   ##define threshold for filtering

############################### Create DESeq2 datasets #############################
library(DESeq2)

##dds <- DESeqDataSetFromMatrix(countData = rawcount, colData = anno, design = ~Condition )   ##rawcount
#dds <- DESeqDataSetFromMatrix(countData = rawcount, colData = anno, design =  ~Batch+Condition )  ###USE this one if you have extra col in anno data with Batch info
dds = DESeq2::DESeqDataSetFromMatrix(countData = adjusted_counts, colData = anno, design = ~ Condition)  ##https://github.com/zhangyuqing/ComBat-seq/issues/7

##When considering batch effects in group design, it takes into account the mean differences across batch, 
##not necessarily the variance differences. ComBat-Seq is designed to address both mean and variance batch effects.
###In theory, no, you do not need to include batch as a covariate any more. However, you can always try both and evaluate the results.


#View(counts(dds))

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)  ## extract normalization count after executing Deseq2 for visualization purpose
vst <- vst(dds, blind=TRUE)  ### Transform counts for data visualization #options (1) vst (2) rld
plotPCA(vst, intgroup="Condition")  ### Plot PCA 

## Run DESEQ2
dds <- DESeq(dds)

##ensure your data is a good fit for the DESeq2 model
plotDispEsts(dds)

################# contrast based  comparison ##########################

#In case of multiple comparisons ## we need to change the contrast for every comparision
contrast<- c("Condition",firstC,SecondC)

res <- results(dds, contrast=contrast)  ## extract result dataframe 
View(as.data.frame(res))

### Valcono plot
library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')   ## Default cut-off for log2FC is >|2| and for P value is 10e-6. USE  pCutoff = 10e-6, FCcutoff = 2.0 


res$threshold <- as.logical(res$padj < p.threshold)  #Threshold defined earlier

nam <- paste('down_in',firstC, sep = '_')
#res$nam <- as.logical(res$log2FoldChange < 0)
res[, nam] <- as.logical(res$log2FoldChange < 0)

genes.deseq <- row.names(res)[which(res$threshold)]   ### list of gene with Padjust < defined threshold
genes_deseq2_sig <- res[which(res$threshold),]


########### Plots normalized count of top 20 genes ## sorted based on padjust and filter by |logFC| >=1

res$gene <- row.names(res)
View(as.data.frame(res))

# Order results by padj values

#library(dplyr)
library(tidyverse)

top20 <- res %>% 
  as.data.frame %>%
  arrange(padj) %>% 	#Arrange rows by padj values
  filter(abs(log2FoldChange) >=1) %>%   #filter based on logFC
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes

top20_norm <- as.data.frame(normalized_counts[rownames(normalized_counts) %in% top20,])

top20_norm_v2 <- top20_norm ## will use later for heatmap

top20_norm <- (top20_norm+1) ## in later step to remove infinity bias due to log

top20_norm$gene <-  row.names(top20_norm)  
top20_norm <- top20_norm %>% 
  pivot_longer(!gene, names_to = "samplename", values_to = "normalized_counts") # Gathering the columns to have normalized counts to a single column

# Create tibbles including row names
mov10_meta <- anno %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

top20_norm <- inner_join(mov10_meta, top20_norm)

################3
## plot using ggplot2

ggplot(top20_norm) +
  geom_point(aes(x = gene, y = normalized_counts, color = Condition)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log 10 CPM Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes with abs(logFC) =>1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))

##################
library(RColorBrewer)
### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap
library(pheatmap) 
pheatmap(top20_norm_v2 , 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         annotation_col = anno[,1:2], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)



file <- paste('Deseq2_',firstC,'_v_',SecondC,'_results_significant_padj',p.threshold,'.csv',sep = '') 
all_results <- paste('Deseq2_',firstC,'_v_',SecondC,'_all_results.csv',sep = '')

write.table(genes_deseq2_sig,all_results,sep = ",")  ## no LogFC threshold


######################  Filter for coding genes (In case want to filter non-coding Genes) ########################
library("biomaRt")

#new_config <- httr::config(ssl_verifypeer = FALSE) ############For certificate error
#httr::set_config(new_config, override = FALSE)     ############For certificate error

### define the mart for h_sapiens

#ensembl_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")  ## either this or following line
ensembl_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")

#all_genes <- getBM(attributes = c( "hgnc_symbol","ensembl_gene_id","ensembl_gene_id_version"),  mart =ensembl_mart)  ## etract df of verious types of ID

#############################################################3

### Add EntrezID column to results dataframe for easier downstream processing ####
genes_deseq2_sig <- as.data.frame(genes_deseq2_sig)
genes_deseq2_sig$hgnc_symbol = row.names(genes_deseq2_sig)  ## significant gene table from previous DE analysis
row.names(genes_deseq2_sig) <- NULL 


genes.deseq.entrezid <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = genes_deseq2_sig$hgnc_symbol, mart = ensembl_mart)
#genes.deseq.entrezid = as.data.frame(genes.deseq.entrezid) ## if not

merged <- merge(genes_deseq2_sig, genes.deseq.entrezid, by.x= "hgnc_symbol", by.y="hgnc_symbol")

##### You may want to filter genes based on LOGFC threshold

merged <- merged[(merged$log2FoldChange >=1 | merged$log2FoldChange <= -1),]
print(merged)

######### Rank all genes based on their fold change #########

#BiocManager::install("clusterProfiler", force = TRUE)
#BiocManager::install("pathview", force = TRUE)
#BiocManager::install("enrichplot", force = TRUE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"

#BiocManager::install(organism, character.only = TRUE, force = TRUE)

library(organism, character.only = TRUE)
keytypes(org.Hs.eg.db)
#We will take the log2FoldChange value from previously saved significant results file
#Deseq2_case1_v_Control_results_significant.csv

 
df <- read.csv("Deseq2_case1_v_Control_results_significant_padj0.05.csv") 
#df <- as.data.frame(merged) 
#print(df)
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
print(original_gene_list)

# name the vector
names(original_gene_list) <- df$entrezgene_id

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

print(gene_list)

#### Gene Set Enrichment #####

library(stats)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENTREZID", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

# require(DOSE)
# dotplot(gse, showCategory=10, split=".sign", orderBy = "X")

### We can see that in our dataset not a single value is enriched at a pvalue cut-off of 0.05.

### Lets exlore other functions with a sample dataset and see what analysis we can do with
## a list of differentially expressed genes
###### geneList dataset of DOSE package #######


##GO Enrichment Analysis of a gene set. 
##Given a vector of genes, enrichGO function will return the 
##enrichment GO categories after FDR control.


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(ggnewscale)
library(DOSE)

data(geneList)
View(geneList)
gene <- names(geneList)[abs(geneList) > 2]
ego <- enrichGO(gene  = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

##### Visualization of enrichGO ######
d <- godata('org.Hs.eg.db', ont="BP")
ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
emapplot(ego2)
emapplot_cluster(ego2)

### Try GO with all different ont methods parameter 
## BP = Biological Processes, CC= Cellular component, MF = Molecular functions

###In the following example, we selected fold change above 1 as the differential genes 
##and analyzing their disease association.

#### enrich DO #####

library(ggupset)


gene = names(geneList)[abs(geneList) > 1.5]
head(gene)
X = enrichDO(gene,ont = "DO", 
              pvalueCutoff=0.05, 
              pAdjustMethod = "BH", 
              universe = names(geneList),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

#The readable is a logical parameter, 
#indicates whether the entrezgene IDs will mapping to gene symbols or not
head(X)

#setReadable function helps to convert entrezgene IDs to gene symbols
X <- setReadable(X, 'org.Hs.eg.db')
head(X)

## Visualization of enrichDO results ##
barplot(X, showCategory=10)
dotplot(X)

#gene may belong to multiple annotation categories, 
#we developed cnetplot function to extract the complex association between genes and diseases
cnetplot(X, categorySize="pvalue", foldChange=geneList)

#upsetplot is an alternative to 
#cnetplot for visualizing the complex association between genes and diseases.

upsetplot(X)

## Enrichment Map ##
###Enrichment map organizes enriched terms into a network with edges 
##connecting overlapping gene sets. 
##In this way, mutually overlapping gene sets are tend to cluster together, 
##making it easy to identify functional modules.

X2 <- pairwise_termsim(X)
emapplot(X2, showCategory = 10, layout = "star")


#Network of Cancer Gene (NCG)3 is a manually curated repository of cancer genes. 
#NCG release 5.0 (Aug. 2015) collects 1,571 cancer genes from 175 published studies. 
#DOSE supports analyzing gene list and 
#determine whether they are enriched in genes known to be mutated in a given cancer type.

######## enrichNCG function ####

gene2 <- names(geneList)[abs(geneList) < 3]
ncg <- enrichNCG(gene2)
head(ncg)

##The enrichment analysis of disease-gene associations is supported by the enrichDGN function 
##to determine whether the genes have associations with any known diseases 
#### gene disease association #####

dgn <- enrichDGN(gene)
head(dgn)

####### Gene set enrichment analysis GSEA Plot ######
gsecc <- gseGO(geneList=geneList, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)
head(summary(gsecc))
gseaplot(gsecc, geneSetID="GO:0000779")


###### KEGG Enrichment Analysis #######

library(clusterProfiler)

## KEGG pathway over-representation analysis
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

## KEGG module over-representation analysis
#KEGG Module is a collection of manually defined function units. In some situation,
#KEGG Modules have a more straightforward interpretation

mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)                   

## KEGG module gene set enrichment analysis ##

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 1)
head(mkk2)


####### Visualize enriched KEGG pathways #########

## To view the KEGG pathway,use the browseKEGG function,
#which will open a web browser and highlight enriched genes.

browseKEGG(kk, 'hsa04110')

###use the pathview() function from the pathview to visualize enriched KEGG 
##pathways identified by the clusterProfiler package

library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
