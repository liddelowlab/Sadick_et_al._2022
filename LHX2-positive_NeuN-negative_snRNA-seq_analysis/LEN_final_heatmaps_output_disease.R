#Generate heatmaps for disease-specific differentially expressed genes (DEGs) in astrocytes or oligodendrocytes in LIM Homeobox 2 (LHX2)-positive/NeuN-negative *FINAL* donor snRNA-seq dataset
#Goal: To generate heatmaps highlighting identified DEGs in either astrocytes or oligodendrocyte by cluster and disease state.
#Pipeline prepared by Jessica S. Sadick

#Load libraries
library(ggplot2)
library(gplots)
library(plyr)
library(dplyr)
library(matrixStats)
library(ComplexHeatmap)
library(Seurat)
library(SingleCellExperiment)

#Set working directory
setwd("file_path")

#-----------------------------------------------------------------------------------------------------
#GENERATE HEATMAPS FOR ASTROCYTE DISEASE-SPECIFIC DEGs BY CLUSTER AND DISEASE STATE
#Read in csv files for astrocyte disease-specific DEGs identified per cluster
Cluster0 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster0-0.25lfc.csv")
Cluster1 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster1-0.25lfc.csv")
Cluster2 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster2-0.25lfc.csv")
Cluster3 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster3-0.25lfc.csv")
Cluster4 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster4-0.25lfc.csv")
Cluster5 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster5-0.25lfc.csv")
Cluster6 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster6-0.25lfc.csv")
Cluster7 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster7-0.25lfc.csv")
Cluster8 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_sig_genes_cluster8-0.25lfc.csv")

##Subset astrocyte disease-specific upregulated DEGs per cluster
C0 <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC >= 1)
C1 <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC >= 1)
C2 <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC >= 1)
C3 <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC >= 1)
C4 <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC >= 1)
C5 <- subset(Cluster5, select = -c(X,Pval), Cluster5$logFC >= 1)
C6 <- subset(Cluster6, select = -c(X,Pval), Cluster6$logFC >= 1)
C7 <- subset(Cluster7, select = -c(X,Pval), Cluster7$logFC >= 1)
C8 <- subset(Cluster8, select = -c(X,Pval), Cluster8$logFC >= 1)
##Subset astrocyte disease-specific downregulated DEGs per cluster
C0 <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC <= -1)
C1 <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC <= -1)
C2 <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC <= -1)
C3 <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC <= -1)
C4 <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC <= -1)
C5 <- subset(Cluster5, select = -c(X,Pval), Cluster5$logFC <= -1)
C6 <- subset(Cluster6, select = -c(X,Pval), Cluster6$logFC <= -1)
C7 <- subset(Cluster7, select = -c(X,Pval), Cluster7$logFC <= -1)
C8 <- subset(Cluster8, select = -c(X,Pval), Cluster8$logFC <= -1)

#Identify the top 10 DEGs that are either upregulated or downregulated in each astrocyte cluster
C0_down <- Cluster0[order(Cluster0$logFC),][1:10,]
C0_up <- Cluster0[rev(order(Cluster0$logFC)),][1:10,]
C1_down <- Cluster1[order(Cluster1$logFC),][1:10,]
C1_up <- Cluster1[rev(order(Cluster1$logFC)),][1:10,]
C2_down <- Cluster2[order(Cluster2$logFC),][1:10,]
C2_up <- Cluster2[rev(order(Cluster2$logFC)),][1:10,]
C3_down <- Cluster3[order(Cluster3$logFC),][1:10,]
C3_up <- Cluster3[rev(order(Cluster3$logFC)),][1:10,]
C4_down <- Cluster4[order(Cluster4$logFC),][1:10,]
C4_up <- Cluster4[rev(order(Cluster4$logFC)),][1:10,]
C5_down <- Cluster5[order(Cluster5$logFC),][1:10,]
C5_up <- Cluster5[rev(order(Cluster5$logFC)),][1:10,]
C6_down <- Cluster6[order(Cluster6$logFC),][1:10,]
C6_up <- Cluster6[rev(order(Cluster6$logFC)),][1:10,]
C7_down <- Cluster7[order(Cluster7$logFC),][1:10,]
C7_up <- Cluster7[rev(order(Cluster7$logFC)),][1:10,]
C8_down <- Cluster8[order(Cluster8$logFC),][1:10,]
C8_up <- Cluster8[rev(order(Cluster8$logFC)),][1:10,]
##Subset DEG name and pval for upregulated DEGs by astrocyte cluster
C0_U <- subset(C0_up, select = -c(X,Pval))
C1_U <- subset(C1_up, select = -c(X,Pval))
C2_U <- subset(C2_up, select = -c(X,Pval))
C3_U <- subset(C3_up, select = -c(X,Pval))
C4_U <- subset(C4_up, select = -c(X,Pval))
C5_U <- subset(C5_up, select = -c(X,Pval))
C6_U <- subset(C6_up, select = -c(X,Pval))
C7_U <- subset(C7_up, select = -c(X,Pval))
C8_U <- subset(C8_up, select = -c(X,Pval))
##Subset DEG name and pval for downregulated DEGs by astrocyte cluster
C0_D <- subset(C0_down, select = -c(X,Pval))
C1_D <- subset(C1_down, select = -c(X,Pval))
C2_D <- subset(C2_down, select = -c(X,Pval))
C3_D <- subset(C3_down, select = -c(X,Pval))
C4_D <- subset(C4_down, select = -c(X,Pval))
C5_D <- subset(C5_down, select = -c(X,Pval))
C6_D <- subset(C6_down, select = -c(X,Pval))
C7_D <- subset(C7_down, select = -c(X,Pval))
C8_D <- subset(C8_down, select = -c(X,Pval))

#Merge lists of DEGs into a dataframe
##For astrocyte upregulated DEGs
DEG_df_up <- C0_U %>% full_join(C1_U, by=c("row.names.genes.")) %>%
  full_join(., C2_U, by=c("row.names.genes.")) %>%
  full_join(., C3_U, by=c("row.names.genes.")) %>%
  full_join(., C4_U, by=c("row.names.genes.")) %>%
  full_join(., C5_U, by=c("row.names.genes.")) %>%
  full_join(., C6_U, by=c("row.names.genes.")) %>%
  full_join(., C7_U, by=c("row.names.genes.")) %>%
  full_join(., C8_U, by=c("row.names.genes."))
##Update column names
colnames(DEG_df_up)<- c("gene_id", "Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4",
                        "Cluster5", "Cluster6", "Cluster7", "Cluster8")

##For astrocyte downregulated DEGs
DEG_df_down <- C0_D %>% full_join(C1_D, by=c("row.names.genes.")) %>%
  full_join(., C2_D, by=c("row.names.genes.")) %>%
  full_join(., C3_D, by=c("row.names.genes.")) %>%
  full_join(., C4_D, by=c("row.names.genes.")) %>%
  full_join(., C5_D, by=c("row.names.genes.")) %>%
  full_join(., C6_D, by=c("row.names.genes.")) %>%
  full_join(., C7_D, by=c("row.names.genes.")) %>%
  full_join(., C8_D, by=c("row.names.genes."))
##Update column names
colnames(DEG_df_down)<- c("gene_id", "Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4",
                          "Cluster5", "Cluster6", "Cluster7", "Cluster8")

#Save row names, convert data to matrix, and add row names to matrix
##For astrocyte upregulated DEGs
rnames_up <- DEG_df_up[,1]
##For astrocyte downregulated DEGs
rnames_down <- DEG_df_down[,1]

#Load SeuratObject and convert to SCE
so_astro_r2 <- readRDS(file.path("file_path", "CR4_LEN_final_so_astro_r2_20PC.rds"))
sce <- as.SingleCellExperiment(so_astro_r2, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
mutate_if(is.character, as.factor) %>%
DataFrame(row.names = colnames(sce))
##Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so_astro_r2 <- SetIdent(so_astro_r2, value = "integrated_snn_res.0.3")
so_astro_r2@meta.data$cluster_id <- Idents(so_astro_r2)
sce$cluster_id <- Idents(so_astro_r2)
##Set ident ##in order to evaluate DEGs by cluster AND by disease state
DefaultAssay(so_astro_r2) <- "RNA"
so_astro_r2$cluster.disease_id <- paste(Idents(so_astro_r2), so_astro_r2$disease_id, sep = "_")
so_astro_r2$cluster <- Idents(so_astro_r2)
Idents(so_astro_r2) <- "cluster.disease_id"
designated_levels <- c("0_NS", "0_AD", "1_NS", "1_AD", "2_NS", "2_AD", "3_NS", "3_AD", 
                       "4_NS", "4_AD", "5_NS", "5_AD", "6_NS", "6_AD", "7_NS", "7_AD", 
                       "8_NS", "8_AD")
Idents(so_astro_r2) <- factor(Idents(so_astro_r2), levels= designated_levels)

#Calculate average expression per cluster/disease state
so_astro_r2_average <- as.data.frame(AverageExpression(so_astro_r2, assays = "RNA"))

#Subset DEGs from average expression dataframe
##For astrocyte upregulated DEGs
so_astro_r2_average_DEG_up <- so_astro_r2_average %>% filter(row.names(so_astro_r2_average) %in% rnames_up)
##For astrocyte downregulated DEGs
so_astro_r2_average_DEG_down <- so_astro_r2_average %>% filter(row.names(so_astro_r2_average) %in% rnames_down)

#Generate z-scores from average expression dataframe
##For astrocyte upregulated DEGs
so_astro_r2_average_DEG_zscore_up <- (so_astro_r2_average_DEG_up-rowMeans(so_astro_r2_average_DEG_up))/(rowSds(as.matrix(so_astro_r2_average_DEG_up)))[row(so_astro_r2_average_DEG_up)]
##For astrocyte downregulated DEGs
so_astro_r2_average_DEG_zscore_down <- (so_astro_r2_average_DEG_down-rowMeans(so_astro_r2_average_DEG_down))/(rowSds(as.matrix(so_astro_r2_average_DEG_down)))[row(so_astro_r2_average_DEG_down)]

#Generate heatmap
##For astrocyte upregulated DEGs
Heatmap(so_astro_r2_average_DEG_zscore_up, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))
##For astrocyte downregulated DEGs
Heatmap(so_astro_r2_average_DEG_zscore_down, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))

#-----------------------------------------------------------------------------------------------------
#GENERATE HEATMAPS FOR ASTROCYTE DISEASE-SPECIFIC UPREGULATED DEGs ASSOCIATED WITH SPECIFIC 'GENE ONTOLOGY (GO) DESCRIPTIONS' BY CLUSTER AND DISEASE STATE
#Read in csv with astrocyte disease-specific upregulated gene lists by cluster ##manually curated by cross-comparing GO terms and associated DEGs
DEG_up <- read.csv(file = "./Astro_diseaseUP_features_by_cluster.csv")
##Make lists of DEGs for specific cluster heatmaps
rname_C0 <- DEG_up[,1]
rname_C1 <- DEG_up[,2]
rname_C2 <- DEG_up[,3]
rname_C3 <- DEG_up[,4]
rname_C4 <- DEG_up[,5]
rname_C5 <- DEG_up[,6]
rname_C6 <- DEG_up[,7]
rname_C0146 <- DEG_up[,8]
rname_C46 <- DEG_up[,9]
rname_C134 <- DEG_up[,10]

#Load SeuratObject and convert to SCE
so_astro_r2 <- readRDS(file.path("file_path", "CR4_LEN_final_so_astro_r2_20PC.rds"))
sce <- as.SingleCellExperiment(so_astro_r2, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
mutate_if(is.character, as.factor) %>%
DataFrame(row.names = colnames(sce))
##Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so_astro_r2 <- SetIdent(so_astro_r2, value = "integrated_snn_res.0.3")
so_astro_r2@meta.data$cluster_id <- Idents(so_astro_r2)
sce$cluster_id <- Idents(so_astro_r2)
##Set ident ##in order to evaluate DEGs by cluster AND by disease state
DefaultAssay(so_astro_r2) <- "RNA"
so_astro_r2$cluster.disease_id <- paste(Idents(so_astro_r2), so_astro_r2$disease_id, sep = "_")
so_astro_r2$cluster <- Idents(so_astro_r2)
Idents(so_astro_r2) <- "cluster.disease_id"
designated_levels <- c("0_NS", "0_AD", "1_NS", "1_AD", "2_NS", "2_AD", "3_NS", "3_AD",
"4_NS", "4_AD", "5_NS", "5_AD", "6_NS", "6_AD", "7_NS", "7_AD",
"8_NS", "8_AD")
Idents(so_astro_r2) <- factor(Idents(so_astro_r2), levels= designated_levels)

#Calculate average expression per cluster/disease state
so_astro_r2_average <- as.data.frame(AverageExpression(so_astro_r2, assays = "RNA"))

#Subset DEGs from average expression dataframe for each cluster list of DEGs ##repeat to cover all clusters
so_astro_r2_average_DEG <- so_astro_r2_average %>% filter(row.names(so_astro_r2_average) %in% rname_C0)

#Generate z-scores from average expression dataframe
so_astro_r2_average_DEG_zscore <- (so_astro_r2_average_DEG-rowMeans(so_astro_r2_average_DEG))/(rowSds(as.matrix(so_astro_r2_average_DEG)))[row(so_astro_r2_average_DEG)]

#Generate heatmap
Heatmap(so_astro_r2_average_DEG_zscore, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))

#Highlight specific DEGs of interest using violin plot function
##Astrocyte disease upregulated: Unique C1 == cell death/oxidative stress
C1_up <- subset(x = so_astro_r2, idents = c("1_NS", "1_AD"), invert = FALSE)
VlnPlot(C1_up, features = c("CST3"), pt.size = 0) + theme(aspect.ratio = 1)

##Astrocyte disease upregulated: Unique C5 == lipid localization
C5_up <- subset(x = so_astro_r2, idents = c("5_NS", "5_AD"), invert = FALSE)
VlnPlot(C5_up, features = c("ABCA1"), pt.size = 0) + theme(aspect.ratio = 1)

##Astrocyte disease upregulated: Common C46 == protein folding/unfolding
C46_up <- subset(x = so_astro_r2, idents = c("4_NS", "4_AD", "6_NS", "6_AD"), invert = FALSE)
VlnPlot(C46_up, features = c("DNAJB1"), pt.size = 0) + theme(aspect.ratio = 1)

##Astrocyte disease downregulated: Common C0146 == response to metal ions
C0146_up <- subset(x = so_astro_r2, idents = c("0_NS", "0_AD", "1_NS", "1_AD", "4_NS", "4_AD", "6_NS", "6_AD"), invert = FALSE)
VlnPlot(C0146_up, features = c("JUND"), pt.size = 0) + theme(aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------
#GENERATE HEATMAPS FOR ASTROCYTE DISEASE-SPECIFIC DOWNREGULATED DEGs ASSOCIATED WITH SPECIFIC 'GO DESCRIPTIONS' BY CLUSTER AND DISEASE STATE
#Read in csv with astrocyte disease-specific downregulated gene lists by cluster ##manually curated by cross-comparing GO terms and associated DEGs
DEG_down <- read.csv(file = "./Astro_diseaseDOWN_features_by_cluster.csv")
##Make lists of DEGs for specific cluster heatmaps
rname_C2 <- DEG_down[,1]
rname_C3 <- DEG_down[,2]
rname_C4 <- DEG_down[,3]
rname_C5 <- DEG_down[,4]
rname_C46 <- DEG_down[,5]
rname_C2456 <- DEG_down[,6]
rname_C2357 <- DEG_down[,7]

#Load SeuratObject and convert to SCE
so_astro_r2 <- readRDS(file.path("file_path", "CR4_LEN_final_so_astro_r2_20PC.rds"))
sce <- as.SingleCellExperiment(so_astro_r2, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
mutate_if(is.character, as.factor) %>%
DataFrame(row.names = colnames(sce))
##Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so_astro_r2 <- SetIdent(so_astro_r2, value = "integrated_snn_res.0.3")
so_astro_r2@meta.data$cluster_id <- Idents(so_astro_r2)
sce$cluster_id <- Idents(so_astro_r2)
##Set ident ##in order to evaluate DEGs by cluster AND by disease state
DefaultAssay(so_astro_r2) <- "RNA"
so_astro_r2$cluster.disease_id <- paste(Idents(so_astro_r2), so_astro_r2$disease_id, sep = "_")
so_astro_r2$cluster <- Idents(so_astro_r2)
Idents(so_astro_r2) <- "cluster.disease_id"
designated_levels <- c("0_NS", "0_AD", "1_NS", "1_AD", "2_NS", "2_AD", "3_NS", "3_AD",
"4_NS", "4_AD", "5_NS", "5_AD", "6_NS", "6_AD", "7_NS", "7_AD",
"8_NS", "8_AD")
Idents(so_astro_r2) <- factor(Idents(so_astro_r2), levels= designated_levels)

#Calculate average expression per cluster/disease state
so_astro_r2_average <- as.data.frame(AverageExpression(so_astro_r2, assays = "RNA"))

#Subset DEGs from average expression dataframe for each cluster list of DEGs ##repeat to cover all clusters
so_astro_r2_average_DEG <- so_astro_r2_average %>% filter(row.names(so_astro_r2_average) %in% rname_C2)

#Generate z-scores from average expression dataframe
so_astro_r2_average_DEG_zscore <- (so_astro_r2_average_DEG-rowMeans(so_astro_r2_average_DEG))/(rowSds(as.matrix(so_astro_r2_average_DEG)))[row(so_astro_r2_average_DEG)]

#Generate heatmap
Heatmap(so_astro_r2_average_DEG_zscore, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))

#Highlight specific DEGs of interest using violin plot function
##Astrocyte disease downregulated: Unique C3 == angiogenesis regulation/ BBB maintenance
C3_down <- subset(x = so_astro_r2, idents = c("3_NS", "3_AD"), invert = FALSE)
VlnPlot(C3_down, features = c("PML"), pt.size = 0) + theme(aspect.ratio = 1)

##Astrocyte disease downregulated: Unique C5 == ERK1/2 cascade signaling
C5_down <- subset(x = so_astro_r2, idents = c("5_NS", "5_AD"), invert = FALSE)
VlnPlot(C5_down, features = c("RIPK2"), pt.size = 0) + theme(aspect.ratio = 1)

##Astrocyte disease downregulated: Common C46 == receptor signaling and synaptic/axonal regulation
C46_down <- subset(x = so_astro_r2, idents = c("4_NS", "4_AD", "6_NS", "6_AD"), invert = FALSE)
VlnPlot(C46_down, features = c("SHISA6"), pt.size = 0) + theme(aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------
#GENERATE HEATMAPS FOR OLIGODENDROCYTE DISEASE-SPECIFIC DEGs BY CLUSTER AND DISEASE STATE
#Read in csv files for oligodendrocyte disease-specific DEGs identified per cluster
Cluster0 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster0-0.25lfc.csv")
Cluster1 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster1-0.25lfc.csv")
Cluster2 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster2-0.25lfc.csv")
Cluster3 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster3-0.25lfc.csv")
Cluster4 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster4-0.25lfc.csv")

##Subset oligodendrocyte disease-specific upregulated DEGs per cluster
C0 <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC >= 1)
C1 <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC >= 1)
C2 <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC >= 1)
C3 <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC >= 1)
C4 <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC >= 1)
##Subset oligodendrocyte disease-specific downregulated DEGs per cluster
C0 <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC <= -1)
C1 <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC <= -1)
C2 <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC <= -1)
C3 <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC <= -1)
C4 <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC <= -1)

#Identify the top 10 DEGs that are either upregulated or downregulated in each oligodendrocyte cluster
C0_down <- Cluster0[order(Cluster0$logFC),][1:10,]
C0_up <- Cluster0[rev(order(Cluster0$logFC)),][1:10,]
C1_down <- Cluster1[order(Cluster1$logFC),][1:10,]
C1_up <- Cluster1[rev(order(Cluster1$logFC)),][1:10,]
C2_down <- Cluster2[order(Cluster2$logFC),][1:10,]
C2_up <- Cluster2[rev(order(Cluster2$logFC)),][1:10,]
C3_down <- Cluster3[order(Cluster3$logFC),][1:10,]
C3_up <- Cluster3[rev(order(Cluster3$logFC)),][1:10,]
C4_down <- Cluster4[order(Cluster4$logFC),][1:10,]
C4_up <- Cluster4[rev(order(Cluster4$logFC)),][1:10,]
##Subset DEG name and pval for upregulated DEGs by oligodendrocyte cluster
C0_U <- subset(C0_up, select = -c(X,Pval))
C1_U <- subset(C1_up, select = -c(X,Pval))
C2_U <- subset(C2_up, select = -c(X,Pval))
C3_U <- subset(C3_up, select = -c(X,Pval))
C4_U <- subset(C4_up, select = -c(X,Pval))
##Subset DEG name and pval for downregulated DEGs by oligodendrocyte cluster
C0_D <- subset(C0_down, select = -c(X,Pval))
C1_D <- subset(C1_down, select = -c(X,Pval))
C2_D <- subset(C2_down, select = -c(X,Pval))
C3_D <- subset(C3_down, select = -c(X,Pval))
C4_D <- subset(C4_down, select = -c(X,Pval))

#Merge lists of DEGs into a dataframe
##For oligodendrocyte upregulated DEGs
DEG_df_up <- C0_U %>% full_join(C1_U, by=c("row.names.genes.")) %>%
  full_join(., C2_U, by=c("row.names.genes.")) %>%
  full_join(., C3_U, by=c("row.names.genes.")) %>%
  full_join(., C4_U, by=c("row.names.genes."))
##Update column names
colnames(DEG_df_up)<- c("gene_id", "Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4")

##For oligodendrocyte downregulated DEGs
DEG_df_down <- C0_D %>% full_join(C1_D, by=c("row.names.genes.")) %>%
  full_join(., C2_D, by=c("row.names.genes.")) %>%
  full_join(., C3_D, by=c("row.names.genes.")) %>%
  full_join(., C4_D, by=c("row.names.genes."))
##Update column names
colnames(DEG_df_down)<- c("gene_id", "Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4")

#Save row names, convert data to matrix, and add row names to matrix
##For oligodendrocyte upregulated DEGs
rnames_up <- DEG_df_up[,1]
##For oligodendrocyte downregulated DEGs
rnames_down <- DEG_df_down[,1]

#Load SeuratObject and convert to SCE
so_oligo <- readRDS(file.path("file_path", "CR4_LEN_final_so_oligo_15PC.rds"))
sce <- as.SingleCellExperiment(so_oligo, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
mutate_if(is.character, as.factor) %>%
DataFrame(row.names = colnames(sce))
##Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so_oligo <- SetIdent(so_oligo, value = "integrated_snn_res.0.1")
so_oligo@meta.data$cluster_id <- Idents(so_oligo)
sce$cluster_id <- Idents(so_oligo)
##Set ident ##in order to evaluate DEGs by cluster AND by disease state
DefaultAssay(so_oligo) <- "RNA"
so_oligo$cluster.disease_id <- paste(Idents(so_oligo), so_oligo$disease_id, sep = "_")
so_oligo$cluster <- Idents(so_oligo)
Idents(so_oligo) <- "cluster.disease_id"
designated_levels <- c("0_NS", "0_AD", "1_NS", "1_AD", "2_NS", "2_AD", "3_NS", "3_AD", 
                       "4_NS", "4_AD")
Idents(so_oligo) <- factor(Idents(so_oligo), levels= designated_levels)

#Calculate average expression per cluster/disease state
so_oligo_average <- as.data.frame(AverageExpression(so_oligo, assays = "RNA"))

#Subset DEGs from average expression dataframe
##For oligodendrocyte upregulated DEGs
so_oligo_average_DEG_up <- so_oligo_average %>% filter(row.names(so_oligo_average) %in% rnames_up)
##For oligodendrocyte downregulated DEGs
so_oligo_average_DEG_down <- so_oligo_average %>% filter(row.names(so_oligo_average) %in% rnames_down)

#Generate z-scores from average expression dataframe
##For oligodendrocyte upregulated DEGs
so_oligo_average_DEG_zscore_up <- (so_oligo_average_DEG_up-rowMeans(so_oligo_average_DEG_up))/(rowSds(as.matrix(so_oligo_average_DEG_up)))[row(so_oligo_average_DEG_up)]
##For oligodendrocyte downregulated DEGs
so_oligo_average_DEG_zscore_down <- (so_oligo_average_DEG_down-rowMeans(so_oligo_average_DEG_down))/(rowSds(as.matrix(so_oligo_average_DEG_down)))[row(so_oligo_average_DEG_down)]

#Generate heatmap
##For oligodendrocyte upregulated DEGs
Heatmap(so_oligo_average_DEG_zscore_up, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))
##For oligodendrocyte downregulated DEGs
Heatmap(so_oligo_average_DEG_zscore_down, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))

#-----------------------------------------------------------------------------------------------------
#GENERATE HEATMAPS FOR OLIGODENDROCYTE DISEASE-SPECIFIC UPREGULATED DEGs ASSOCIATED WITH SPECIFIC 'GO DESCRIPTIONS' BY CLUSTER AND DISEASE STATE
#Read in csv with oligodendrocyte disease-specific upregulated gene lists by cluster ##manually curated by cross-comparing GO terms and associated DEGs
DEG_up <- read.csv(file = "./oligo_diseaseUP_features_by_cluster.csv")
##Make lists of DEGs for specific cluster heatmaps
rname_C1 <- DEG_up[,1]
rname_C2_chol <- DEG_up[,2]
rname_C2_syndep <- DEG_up[,3]
rname_C2_gliadev <- DEG_up[,4]
rname_C3 <- DEG_up[,5]
rname_C4 <- DEG_up[,6]
rname_C14 <- DEG_up[,7]
rname_C12 <- DEG_up[,8]
rname_C24 <- DEG_up[,9]
rname_C02 <- DEG_up[,10]

#Load SeuratObject and convert to SCE
so_oligo <- readRDS(file.path("file_path", "CR4_LEN_final_so_oligo_15PC.rds"))
sce <- as.SingleCellExperiment(so_oligo, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
mutate_if(is.character, as.factor) %>%
DataFrame(row.names = colnames(sce))
##Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so_oligo <- SetIdent(so_oligo, value = "integrated_snn_res.0.1")
so_oligo@meta.data$cluster_id <- Idents(so_oligo)
sce$cluster_id <- Idents(so_oligo)
#Set ident ##in order to evaluate DEGs by cluster AND by disease state
DefaultAssay(so_oligo) <- "RNA"
so_oligo$cluster.disease_id <- paste(Idents(so_oligo), so_oligo$disease_id, sep = "_")
so_oligo$cluster <- Idents(so_oligo)
Idents(so_oligo) <- "cluster.disease_id"
designated_levels <- c("0_NS", "0_AD", "1_NS", "1_AD", "2_NS", "2_AD", "3_NS", "3_AD",
"4_NS", "4_AD", "5_NS", "5_AD", "6_NS", "6_AD", "7_NS", "7_AD",
"8_NS", "8_AD")
Idents(so_oligo) <- factor(Idents(so_oligo), levels= designated_levels)

#Calculate average expression per cluster/disease state
so_oligo_average <- as.data.frame(AverageExpression(so_oligo, assays = "RNA"))

#Subset DEGs from average expression dataframe for each cluster list of DEGs ##repeat to cover all clusters
so_oligo_average_DEG <- so_oligo_average %>% filter(row.names(so_oligo_average) %in% rname_C1)

#Generate zscores for average expression dataframe
so_oligo_average_DEG_zscore <- (so_oligo_average_DEG-rowMeans(so_oligo_average_DEG))/(rowSds(as.matrix(so_oligo_average_DEG)))[row(so_oligo_average_DEG)]

#Generate heatmap
Heatmap(so_oligo_average_DEG_zscore, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))

#Highlight specific DEGs of interest using violin plot function
##Oligodendrocyte disease upregulated: Unique C1 == axonogenesis/synapse organization
C1_up <- subset(x = so_oligo, idents = c("1_NS", "1_AD"), invert = FALSE)
VlnPlot(C1_up, features = c("DCC"), pt.size = 0) + theme(aspect.ratio = 1)

##Oligodendrocyte disease upregulated: Unique C2 == cholesterol metabolism
C2_up <- subset(x = so_oligo, idents = c("2_NS", "2_AD"), invert = FALSE)
VlnPlot(C2_up, features = c("LDLR"), pt.size = 0) + theme(aspect.ratio = 1)

##Oligodendrocyte disease upregulated: Common C02 == regulation of GTPase-mediated signaling
C02_up <- subset(x = so_oligo, idents = c("0_NS", "0_AD", "2_NS", "2_AD"), invert = FALSE)
VlnPlot(C02_up, features = c("RELN"), pt.size = 0) + theme(aspect.ratio = 1)

#-----------------------------------------------------------------------------------------------------
#GENERATE HEATMAPS FOR OLIGODENDROCYTE DISEASE-SPECIFIC DOWNREGULATED DEGs ASSOCIATED WITH SPECIFIC 'GO DESCRIPTIONS' BY CLUSTER AND DISEASE STATE
#Read in csv with oligodendrocyte disease-specific downregulated gene lists by cluster ##manually curated by cross-comparing GO terms and associated DEGs
DEG_down <- read.csv(file = "./Oligo_diseaseDOWN_features_by_cluster.csv")
##Make lists of DEGs for specific cluster heatmaps
rname_C0_synapse <- DEG_down[,1]
rname_C0_metab <- DEG_down[,2]
rname_C0_glialdiff <- DEG_down[,3]
rname_C3 <- DEG_down[,4]
rname_C4 <- DEG_down[,5]
rname_C02_synapse <- DEG_down[,6]
rname_C02_AAsyn <- DEG_down[,7]
rname_C03 <- DEG_down[,8]
rname_C04 <- DEG_down[,9]
rname_C0234 <- DEG_down[,10]

#Load SeuratObject and convert to SCE
so_oligo <- readRDS(file.path("file_path", "CR4_LEN_final_so_oligo_15PC.rds"))
sce <- as.SingleCellExperiment(so_oligo, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
mutate_if(is.character, as.factor) %>%
DataFrame(row.names = colnames(sce))
##Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)
so_oligo <- SetIdent(so_oligo, value = "integrated_snn_res.0.1")
so_oligo@meta.data$cluster_id <- Idents(so_oligo)
sce$cluster_id <- Idents(so_oligo)
#Set ident ##in order to evaluate DEGs by cluster AND by disease state
DefaultAssay(so_oligo) <- "RNA"
so_oligo$cluster.disease_id <- paste(Idents(so_oligo), so_oligo$disease_id, sep = "_")
so_oligo$cluster <- Idents(so_oligo)
Idents(so_oligo) <- "cluster.disease_id"
designated_levels <- c("0_NS", "0_AD", "1_NS", "1_AD", "2_NS", "2_AD", "3_NS", "3_AD",
"4_NS", "4_AD", "5_NS", "5_AD", "6_NS", "6_AD", "7_NS", "7_AD",
"8_NS", "8_AD")
Idents(so_oligo) <- factor(Idents(so_oligo), levels= designated_levels)

#Calculate average expression per cluster/disease state
so_oligo_average <- as.data.frame(AverageExpression(so_oligo, assays = "RNA"))

#Subset DEGs from average expression dataframe for each cluster list of DEGs ##repeat to cover all clusters
so_oligo_average_DEG <- so_oligo_average %>% filter(row.names(so_oligo_average) %in% rname_C0_synapse)

#Generate zscores for average expression dataframe
so_oligo_average_DEG_zscore <- (so_oligo_average_DEG-rowMeans(so_oligo_average_DEG))/(rowSds(as.matrix(so_oligo_average_DEG)))[row(so_oligo_average_DEG)]

#Generate heatmap
Heatmap(so_oligo_average_DEG_zscore, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"))

#Highlight specific DEGs of interest using violin plot function
##Oligodendrocyte disease downregulated: Unique C0 == synapse and metabolism
C0_down <- subset(x = so_oligo, idents = c("0_NS", "0_AD"), invert = FALSE)
VlnPlot(C0_down, features = c("CDH1"), pt.size = 0) + theme(aspect.ratio = 1)

##Oligodendrocyte disease downregulated: Common C0 and C2 == AA synthesis
C02_down <- subset(x = so_oligo, idents = c("0_NS", "0_AD", "2_NS", "2_AD"), invert = FALSE)
VlnPlot(C02_down, features = c("FOLH1"), pt.size = 0) + theme(aspect.ratio = 1)
