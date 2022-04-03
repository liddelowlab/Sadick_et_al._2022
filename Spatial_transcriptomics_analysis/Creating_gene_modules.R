# This script provides the code necessary for the generation of astrocyte cluster gene signature modules based on our single-nucleus RNA-seq results from querying the Hasel et al. (2021) and Maynard et al. (2021) spatial transcriptomics datasets  

# load required packages
library(dplyr)
library(readxl)
library(Seurat)
library(biomaRt)

### get package requirements file
pkgs <- loadedNamespaces()
for (i in 1:length(pkgs)) {
  cat(print(paste(pkgs[i], getNamespaceVersion(pkgs[i]), sep = " ")), file = "Creating_gene_modules_packages.txt", append = TRUE)
  cat("\n", file = "Creating_gene_modules_packages.txt", append = TRUE)
}

# load preprocessed Hasel et al Seurat objects
cnt.A <- readRDS("file_path/hasel_cntA_spatial_seurat_object.rds")
cnt.B <- readRDS("file_path/hasel_cntB_spatial_seurat_object.rds")
cnt.C <- readRDS("file_path/hasel_cntC_spatial_seurat_object.rds")

lps.A <- readRDS("file_path/hasel_lpsA_spatial_seurat_object.rds")
lps.B <- readRDS("file_path/hasel_lpsB_spatial_seurat_object.rds")
lps.C <- readRDS("file_path/hasel_lpsC_spatial_seurat_object.rds")

# retrieve marts 
human <- useEnsembl(biomart = "ensembl", 
                    dataset = "hsapiens_gene_ensembl", 
                    mirror = "useast")
mouse <- useEnsembl(biomart = "ensembl", 
                    dataset = "mmusculus_gene_ensembl", 
                    mirror = "uswest")

# function for converting human gene symbol to mouse gene symbol; retrieves one-to-one, one-to-many, and many-to-one orthologs
## please note: biomaRt ortholog annotations may change over time; orthologous gene modules generated at the time of our analysis can found in Table S7
convertHumanGeneList <- function(x){
  require(biomaRt)
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T, )
  return(genesV2)
}

# get vector of gene names from all Seurat objects
all_genes <- unique(c(rownames(cnt.A), rownames(cnt.B), rownames(cnt.C), rownames(lps.A), rownames(lps.B), rownames(lps.C)))

# create marker genes lists: (using excel sheet with data from Table S6)
# cluster 0
c0 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 1)
# filter out only genes with average log2FC greater than or equal to 0.30
c0.thresh <- c0 %>% filter(avg_log2FC >= 0.3)
c0.genes.df <- convertHumanGeneList(c0.thresh$gene)
colnames(c0.genes.df) <- c("gene", "mus_gene")
c0.genes.df <- merge(c0.thresh, c0.genes.df, by = "gene", all.x = TRUE)
c0.genes.df <- c0.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c0.genes.df)) 
c0.genes.df.o <- c0.genes.df
c0.genes.df <- c0.genes.df[!(duplicated(c0.genes.df$gene) | duplicated(c0.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c0.genes.df <- c0.genes.df[!(duplicated(c0.genes.df$mus_gene) | duplicated(c0.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c0.genes <- c0.genes.df$mus_gene
c0.genes <- c0.genes[!is.na(c0.genes)]
c0.mod.genes <- c()
for (x in 1:length(c0.genes)){
  if (c0.genes[x] %in% all_genes){
    c0.mod.genes <- c(c0.mod.genes, c0.genes[x])
  }
}
c0.mod.genes 


# cluster 1
c1 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 2)
# filter out only genes with average log2FC greater than or equal to 0.30
c1.thresh <- c1 %>% filter(avg_log2FC >= 0.3)
c1.genes.df <- convertHumanGeneList(c1.thresh$gene)
colnames(c1.genes.df) <- c("gene", "mus_gene")
c1.genes.df <- merge(c1.thresh, c1.genes.df, by = "gene", all.x = TRUE)
c1.genes.df <- c1.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c1.genes.df)) # 
c1.genes.df.o <- c1.genes.df
c1.genes.df <- c1.genes.df[!(duplicated(c1.genes.df$gene) | duplicated(c1.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c1.genes.df <- c1.genes.df[!(duplicated(c1.genes.df$mus_gene) | duplicated(c1.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c1.genes <- c1.genes.df$mus_gene
c1.genes <- c1.genes[!is.na(c1.genes)]
c1.mod.genes <- c()
for (x in 1:length(c1.genes)){
  if (c1.genes[x] %in% all_genes){
    c1.mod.genes <- c(c1.mod.genes, c1.genes[x])
  }
}
c1.mod.genes

### 

# cluster 2
c2 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 3)
# filter out only genes with average log2FC greater than or equal to 0.30
c2.thresh <- c2 %>% filter(avg_log2FC >= 0.3)
c2.genes.df <- convertHumanGeneList(c2.thresh$gene)
colnames(c2.genes.df) <- c("gene", "mus_gene")
c2.genes.df <- merge(c2.thresh, c2.genes.df, by = "gene", all.x = TRUE)
c2.genes.df <- c2.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c2.genes.df)) 
c2.genes.df.o <- c2.genes.df
c2.genes.df <- c2.genes.df[!(duplicated(c2.genes.df$gene) | duplicated(c2.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c2.genes.df <- c2.genes.df[!(duplicated(c2.genes.df$mus_gene) | duplicated(c2.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c2.genes <- c2.genes.df$mus_gene
c2.genes <- c2.genes[!is.na(c2.genes)]
c2.mod.genes <- c()
for (x in 1:length(c2.genes)){
  if (c2.genes[x] %in% all_genes){
    c2.mod.genes <- c(c2.mod.genes, c2.genes[x])
  }
}
c2.mod.genes 

# cluster 3
c3 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 4)
# filter out only genes with average log2FC greater than or equal to 0.30
c3.thresh <- c3 %>% filter(avg_log2FC >= 0.3)
c3.genes.df <- convertHumanGeneList(c3.thresh$gene)
colnames(c3.genes.df) <- c("gene", "mus_gene")
c3.genes.df <- merge(c3.thresh, c3.genes.df, by = "gene", all.x = TRUE)
c3.genes.df <- c3.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c3.genes.df)) 
c3.genes.df.o <- c3.genes.df
c3.genes.df <- c3.genes.df[!(duplicated(c3.genes.df$gene) | duplicated(c3.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c3.genes.df <- c3.genes.df[!(duplicated(c3.genes.df$mus_gene) | duplicated(c3.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c3.genes <- c3.genes.df$mus_gene
c3.genes <- c3.genes[!is.na(c3.genes)]
c3.mod.genes <- c()
for (x in 1:length(c3.genes)){
  if (c3.genes[x] %in% all_genes){
    c3.mod.genes <- c(c3.mod.genes, c3.genes[x])
  }
}
c3.mod.genes 


# cluster 4
c4 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 5)
# filter out only genes with average log2FC greater than or equal to 0.30
c4.thresh <- c4 %>% filter(avg_log2FC >= 0.3)
c4.genes.df <- convertHumanGeneList(c4.thresh$gene)
colnames(c4.genes.df) <- c("gene", "mus_gene")
c4.genes.df <- merge(c4.thresh, c4.genes.df, by = "gene", all.x = TRUE)
c4.genes.df <- c4.genes.df[c("gene", "mus_gene")]
c4.genes.df.o <- c4.genes.df
colSums(!is.na(c4.genes.df)) 
c4.genes.df <- c4.genes.df[!(duplicated(c4.genes.df$gene) | duplicated(c4.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c4.genes.df <- c4.genes.df[!(duplicated(c4.genes.df$mus_gene) | duplicated(c4.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c4.genes <- c4.genes.df$mus_gene
c4.genes <- c4.genes[!is.na(c4.genes)]
c4.mod.genes <- c()
for (x in 1:length(c4.genes)){
  if (c4.genes[x] %in% all_genes){
    c4.mod.genes <- c(c4.mod.genes, c4.genes[x])
  }
}
c4.mod.genes 


# cluster 5
c5 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 6)
# filter out only genes with average log2FC greater than or equal to 0.30
c5.thresh <- c5 %>% filter(avg_log2FC >= 0.3)
c5.genes.df <- convertHumanGeneList(c5.thresh$gene)
colnames(c5.genes.df) <- c("gene", "mus_gene")
c5.genes.df <- merge(c5.thresh, c5.genes.df, by = "gene", all.x = TRUE)
c5.genes.df <- c5.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c5.genes.df)) 
c5.genes.df.o <- c5.genes.df
c5.genes.df <- c5.genes.df[!(duplicated(c5.genes.df$gene) | duplicated(c5.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c5.genes.df <- c5.genes.df[!(duplicated(c5.genes.df$mus_gene) | duplicated(c5.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c5.genes <- c5.genes.df$mus_gene
c5.genes <- c5.genes[!is.na(c5.genes)]
c5.mod.genes <- c()
for (x in 1:length(c5.genes)){
  if (c5.genes[x] %in% all_genes){
    c5.mod.genes <- c(c5.mod.genes, c5.genes[x])
  }
}
c5.mod.genes 


# cluster 6
c6 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 7)
# filter out only genes with average log2FC greater than or equal to 0.30
c6.thresh <- c6 %>% filter(avg_log2FC >= 0.3)
c6.genes.df <- convertHumanGeneList(c6.thresh$gene)
colnames(c6.genes.df) <- c("gene", "mus_gene")
c6.genes.df <- merge(c6.thresh, c6.genes.df, by = "gene", all.x = TRUE)
c6.genes.df <- c6.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c6.genes.df)) 
c6.genes.df.o <- c6.genes.df
c6.genes.df <- c6.genes.df[!(duplicated(c6.genes.df$gene) | duplicated(c6.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c6.genes.df <- c6.genes.df[!(duplicated(c6.genes.df$mus_gene) | duplicated(c6.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c6.genes <- c6.genes.df$mus_gene
c6.genes <- c6.genes[!is.na(c6.genes)]
c6.mod.genes <- c()
for (x in 1:length(c6.genes)){
  if (c6.genes[x] %in% all_genes){
    c6.mod.genes <- c(c6.mod.genes, c6.genes[x])
  }
}
c6.mod.genes 

# cluster 7
c7 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 8)
c7.thresh <- c7 %>% filter(avg_log2FC >= 0.3)
c7.genes.df <- convertHumanGeneList(c7.thresh$gene)
colnames(c7.genes.df) <- c("gene", "mus_gene")
c7.genes.df <- merge(c7.thresh, c7.genes.df, by = "gene", all.x = TRUE)
c7.genes.df <- c7.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c7.genes.df))
c7.genes.df.o <- c7.genes.df
c7.genes.df <- c7.genes.df[!(duplicated(c7.genes.df$gene) | duplicated(c7.genes.df$gene, fromLast = TRUE)), ] 
c7.genes.df <- c7.genes.df[!(duplicated(c7.genes.df$mus_gene) | duplicated(c7.genes.df$mus_gene, fromLast = TRUE)), ] 
c7.genes <- c7.genes.df$mus_gene
c7.genes <- c7.genes[!is.na(c7.genes)]
c7.mod.genes <- c()
for (x in 1:length(c7.genes)){
  if (c7.genes[x] %in% all_genes){
    c7.mod.genes <- c(c7.mod.genes, c7.genes[x])
  }
}
c7.mod.genes 


# cluster 8
c8 = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 9)
# filter out only genes with average log2FC greater than or equal to 0.30
c8.thresh <- c8 %>% filter(avg_log2FC >= 0.3)
c8.genes.df <- convertHumanGeneList(c8.thresh$gene)
colnames(c8.genes.df) <- c("gene", "mus_gene")
c8.genes.df <- merge(c8.thresh, c8.genes.df, by = "gene", all.x = TRUE)
c8.genes.df <- c8.genes.df[c("gene", "mus_gene")]
colSums(!is.na(c8.genes.df)) 
c8.genes.df.o <- c8.genes.df
c8.genes.df <- c8.genes.df[!(duplicated(c8.genes.df$gene) | duplicated(c8.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
c8.genes.df <- c8.genes.df[!(duplicated(c8.genes.df$mus_gene) | duplicated(c8.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
c8.genes <- c8.genes.df$mus_gene
c8.genes <- c8.genes[!is.na(c8.genes)]
c8.mod.genes <- c()
for (x in 1:length(c8.genes)){
  if (c8.genes[x] %in% all_genes){
    c8.mod.genes <- c(c8.mod.genes, c8.genes[x])
  }
}
c8.mod.genes 

## AD associated up-regulated genes per cluster: 
ds_genes = read_excel("file_path/sadick_cluster_degs.xlsx", sheet = 15)
colnames(ds_genes) <- c("gene", "pvalue", "log2FC", "cluster")

# setting log2FC threshold of 0.5
clus.0 <- ds_genes %>% filter(cluster == 0) %>% filter(log2FC > 0.5)
clus.1 <- ds_genes %>% filter(cluster == 1) %>% filter(log2FC > 0.5)
clus.2 <- ds_genes %>% filter(cluster == 2) %>% filter(log2FC > 0.5)
clus.3 <- ds_genes %>% filter(cluster == 3) %>% filter(log2FC > 0.5)
clus.4 <- ds_genes %>% filter(cluster == 4) %>% filter(log2FC > 0.5)
clus.5 <- ds_genes %>% filter(cluster == 5) %>% filter(log2FC > 0.5)
clus.6 <- ds_genes %>% filter(cluster == 6) %>% filter(log2FC > 0.5)
clus.7 <- ds_genes %>% filter(cluster == 7) %>% filter(log2FC > 0.5)
clus.8 <- ds_genes %>% filter(cluster == 8) %>% filter(log2FC > 0.5)

# cluster 0
clus.0.genes.df <- convertHumanGeneList(clus.0$gene)
colnames(clus.0.genes.df) <- c("gene", "mus_gene")
clus.0.genes.df <- merge(clus.0, clus.0.genes.df, by = "gene", all.x = TRUE)
clus.0.genes.df.o <- clus.0.genes.df
clus.0.genes.df <- clus.0.genes.df[!(duplicated(clus.0.genes.df$gene) | duplicated(clus.0.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.0.genes.df <- clus.0.genes.df[!(duplicated(clus.0.genes.df$mus_gene) | duplicated(clus.0.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.0.genes <- clus.0.genes.df$mus_gene
clus.0.genes <- clus.0.genes[!is.na(clus.0.genes)]
clus.0.ad.genes <- c()
for (x in 1:length(clus.0.genes)){
  if (clus.0.genes[x] %in% all_genes){
    clus.0.ad.genes <- c(clus.0.ad.genes, clus.0.genes[x])
  }
}
clus.0.ad.genes 

# cluster 1: 
clus.1.genes.df <- convertHumanGeneList(clus.1$gene)
colnames(clus.1.genes.df) <- c("gene", "mus_gene")
clus.1.genes.df <- merge(clus.1, clus.1.genes.df, by = "gene", all.x = TRUE)
clus.1.genes.df.o <- clus.1.genes.df
clus.1.genes.df <- clus.1.genes.df[!(duplicated(clus.1.genes.df$gene) | duplicated(clus.1.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.1.genes.df <- clus.1.genes.df[!(duplicated(clus.1.genes.df$mus_gene) | duplicated(clus.1.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.1.genes <- clus.1.genes.df$mus_gene
clus.1.genes <- clus.1.genes[!is.na(clus.1.genes)]
clus.1.ad.genes <- c()
for (x in 1:length(clus.1.genes)){
  if (clus.1.genes[x] %in% all_genes){
    clus.1.ad.genes <- c(clus.1.ad.genes, clus.1.genes[x])
  }
}
clus.1.ad.genes 

# cluster 2
clus.2.genes.df <- convertHumanGeneList(clus.2$gene)
colnames(clus.2.genes.df) <- c("gene", "mus_gene")
clus.2.genes.df <- merge(clus.2, clus.2.genes.df, by = "gene", all.x = TRUE)
clus.2.genes.df.o <- clus.2.genes.df
clus.2.genes.df <- clus.2.genes.df[!(duplicated(clus.2.genes.df$gene) | duplicated(clus.2.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.2.genes.df <- clus.2.genes.df[!(duplicated(clus.2.genes.df$mus_gene) | duplicated(clus.2.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.2.genes <- clus.2.genes.df$mus_gene
clus.2.genes <- clus.2.genes[!is.na(clus.2.genes)]
clus.2.ad.genes <- c()
for (x in 1:length(clus.2.genes)){
  if (clus.2.genes[x] %in% all_genes){
    clus.2.ad.genes <- c(clus.2.ad.genes, clus.2.genes[x])
  }
}
clus.2.ad.genes

# cluster 3
clus.3.genes.df <- convertHumanGeneList(clus.3$gene)
colnames(clus.3.genes.df) <- c("gene", "mus_gene")
clus.3.genes.df <- merge(clus.3, clus.3.genes.df, by = "gene", all.x = TRUE)
clus.3.genes.df.o <- clus.3.genes.df
clus.3.genes.df <- clus.3.genes.df[!(duplicated(clus.3.genes.df$gene) | duplicated(clus.3.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.3.genes.df <- clus.3.genes.df[!(duplicated(clus.3.genes.df$mus_gene) | duplicated(clus.3.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.3.genes <- clus.3.genes.df$mus_gene
clus.3.genes <- clus.3.genes[!is.na(clus.3.genes)]
clus.3.ad.genes <- c()
for (x in 1:length(clus.3.genes)){
  if (clus.3.genes[x] %in% all_genes){
    clus.3.ad.genes <- c(clus.3.ad.genes, clus.3.genes[x])
  }
}
clus.3.ad.genes

# cluster 4:
clus.4.genes.df <- convertHumanGeneList(clus.4$gene)
colnames(clus.4.genes.df) <- c("gene", "mus_gene")
clus.4.genes.df <- merge(clus.4, clus.4.genes.df, by = "gene", all.x = TRUE)
clus.4.genes.df.o <- clus.4.genes.df
clus.4.genes.df <- clus.4.genes.df[!(duplicated(clus.4.genes.df$gene) | duplicated(clus.4.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.4.genes.df <- clus.4.genes.df[!(duplicated(clus.4.genes.df$mus_gene) | duplicated(clus.4.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.4.genes <- clus.4.genes.df$mus_gene
clus.4.genes <- clus.4.genes[!is.na(clus.4.genes)]
clus.4.ad.genes <- c()
for (x in 1:length(clus.4.genes)){
  if (clus.4.genes[x] %in% all_genes){
    clus.4.ad.genes <- c(clus.4.ad.genes, clus.4.genes[x])
  }
}
clus.4.ad.genes

# cluster 5:
clus.5.genes.df <- convertHumanGeneList(clus.5$gene)
colnames(clus.5.genes.df) <- c("gene", "mus_gene")
clus.5.genes.df <- merge(clus.5, clus.5.genes.df, by = "gene", all.x = TRUE)
clus.5.genes.df.o <- clus.5.genes.df
clus.5.genes.df <- clus.5.genes.df[!(duplicated(clus.5.genes.df$gene) | duplicated(clus.5.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.5.genes.df <- clus.5.genes.df[!(duplicated(clus.5.genes.df$mus_gene) | duplicated(clus.5.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.5.genes <- clus.5.genes.df$mus_gene
clus.5.genes <- clus.5.genes[!is.na(clus.5.genes)]
clus.5.ad.genes <- c()
for (x in 1:length(clus.5.genes)){
  if (clus.5.genes[x] %in% all_genes){
    clus.5.ad.genes <- c(clus.5.ad.genes, clus.5.genes[x])
  }
}
clus.5.ad.genes


# cluster 6:
clus.6.genes.df <- convertHumanGeneList(clus.6$gene)
colnames(clus.6.genes.df) <- c("gene", "mus_gene")
clus.6.genes.df <- merge(clus.6, clus.6.genes.df, by = "gene", all.x = TRUE)
clus.6.genes.df.o <- clus.6.genes.df
clus.6.genes.df <- clus.6.genes.df[!(duplicated(clus.6.genes.df$gene) | duplicated(clus.6.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.6.genes.df <- clus.6.genes.df[!(duplicated(clus.6.genes.df$mus_gene) | duplicated(clus.6.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.6.genes <- clus.6.genes.df$mus_gene
clus.6.genes <- clus.6.genes[!is.na(clus.6.genes)]
clus.6.ad.genes <- c()
for (x in 1:length(clus.6.genes)){
  if (clus.6.genes[x] %in% all_genes){
    clus.6.ad.genes <- c(clus.6.ad.genes, clus.6.genes[x])
  }
}
clus.6.ad.genes

# cluster 7
clus.7.genes.df <- convertHumanGeneList(clus.7$gene)
colnames(clus.7.genes.df) <- c("gene", "mus_gene")
clus.7.genes.df <- merge(clus.7, clus.7.genes.df, by = "gene", all.x = TRUE)
clus.7.genes.df.o <- clus.7.genes.df
clus.7.genes.df <- clus.7.genes.df[!(duplicated(clus.7.genes.df$gene) | duplicated(clus.7.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.7.genes.df <- clus.7.genes.df[!(duplicated(clus.7.genes.df$mus_gene) | duplicated(clus.7.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.7.genes <- clus.7.genes.df$mus_gene
clus.7.genes <- clus.7.genes[!is.na(clus.7.genes)]
clus.7.ad.genes <- c()
for (x in 1:length(clus.7.genes)){
  if (clus.7.genes[x] %in% all_genes){
    clus.7.ad.genes <- c(clus.7.ad.genes, clus.7.genes[x])
  }
}
clus.7.ad.genes 

# cluster 8
clus.8.genes.df <- convertHumanGeneList(clus.8$gene)
colnames(clus.8.genes.df) <- c("gene", "mus_gene")
clus.8.genes.df <- merge(clus.8, clus.8.genes.df, by = "gene", all.x = TRUE)
clus.8.genes.df.o <- clus.8.genes.df
clus.8.genes.df <- clus.8.genes.df[!(duplicated(clus.8.genes.df$gene) | duplicated(clus.8.genes.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
clus.8.genes.df <- clus.8.genes.df[!(duplicated(clus.8.genes.df$mus_gene) | duplicated(clus.8.genes.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
clus.8.genes <- clus.8.genes.df$mus_gene
clus.8.genes <- clus.8.genes[!is.na(clus.8.genes)]
clus.8.ad.genes <- c()
for (x in 1:length(clus.8.genes)){
  if (clus.8.genes[x] %in% all_genes){
    clus.8.ad.genes <- c(clus.8.ad.genes, clus.8.genes[x])
  }
}
clus.8.ad.genes 

# now all gene modules have been created; saving as rds files
final.modules = list("Astro_0" = c0.mod.genes, "Astro_1" = c1.mod.genes, "Astro2" = c2.mod.genes, "Astro3" = c3.mod.genes, "Astro4" = c4.mod.genes, "Astro5" = c5.mod.genes, "Astro6" = c6.mod.genes,
                     "Astro7" = c7.mod.genes, "Astro8" = c8.mod.genes,
                     "Astro_0_AD" = clus.0.ad.genes, "Astro_1_AD" = clus.1.ad.genes, "Astro2_AD" = clus.2.ad.genes, "Astro3_AD" = clus.3.ad.genes, "Astro4_AD" = clus.4.ad.genes, "Astro5_AD" = clus.5.ad.genes, 
                     "Astro6_AD" = clus.6.ad.genes, "Astro7_AD" = clus.7.ad.genes, "Astro8_AD" = clus.8.ad.genes)

saveRDS(final.modules, "file_path/final_mouse_gene_modules_forhaselsections.rds")

final.human.modules = list("Astro_0" = c0.thresh$gene, "Astro_1" = c1.thresh$gene, "Astro2" = c2.thresh$gene, "Astro3" = c3.thresh$gene, "Astro4" = c4.thresh$gene, "Astro5" = c5.thresh$gene, "Astro6" = c6.thresh$gene,
                           "Astro7" = c7.thresh$gene, "Astro8" = c8.thresh$gene)

saveRDS(final.human.modules, "file_path/final_human_gene_modules_formaynardsections.rds")

