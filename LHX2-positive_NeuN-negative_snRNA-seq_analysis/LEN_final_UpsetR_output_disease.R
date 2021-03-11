#Generate UpSetR plots and identify disease-specific gene ontology(GO) terms for astrocyte or oligodendrocyte clusters in LIM Homeobox 2 (LHX2)-positive/NeuN-negative *FINAL* donor snRNA-seq dataset
#Goal: To generate UpSetR plots and to identify disease-specific GO terms in unique to single/multiple astrocyte or oligodendrocyte cluster(s).
#Pipeline prepared by Jessica S. Sadick

#Load library
library(dplyr)
library(UpSetR)
library(ggplot2)
library(plyr)
library(gridExtra)
library(grid)

#Set working directory
setwd("file_path")

#-----------------------------------------------------------------------------------------------------
#GENERATE UpSetR PLOTS FOR ASTROCYTE DISEASE-SPECIFIC DIFFERENTIALLY EXPRESSED GENES (DEGs) BY CLUSTER
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
C0_up <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC >= 0.25)
C1_up <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC >= 0.25)
C2_up <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC >= 0.25)
C3_up <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC >= 0.25)
C4_up <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC >= 0.25)
C5_up <- subset(Cluster5, select = -c(X,Pval), Cluster5$logFC >= 0.25)
C6_up <- subset(Cluster6, select = -c(X,Pval), Cluster6$logFC >= 0.25)
C7_up <- subset(Cluster7, select = -c(X,Pval), Cluster7$logFC >= 0.25)
C8_up <- subset(Cluster8, select = -c(X,Pval), Cluster8$logFC >= 0.25)
##Subset astrocyte disease-specific downregulated DEGs per cluster
C0_down <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC <= -0.25)
C1_down <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC <= -0.25)
C2_down <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC <= -0.25)
C3_down <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC <= -0.25)
C4_down <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC <= -0.25)
C5_down <- subset(Cluster5, select = -c(X,Pval), Cluster5$logFC <= -0.25)
C6_down <- subset(Cluster6, select = -c(X,Pval), Cluster6$logFC <= -0.25)
C7_down <- subset(Cluster7, select = -c(X,Pval), Cluster7$logFC <= -0.25)
C8_down <- subset(Cluster8, select = -c(X,Pval), Cluster8$logFC <= -0.25)

#Subset row names as a vector
##For astrocyte disease-specific upregulated DEGs
C0_U <- as.vector(C0_up$row.names.genes.)
C1_U <- as.vector(C1_up$row.names.genes.)
C2_U <- as.vector(C2_up$row.names.genes.)
C3_U <- as.vector(C3_up$row.names.genes.)
C4_U <- as.vector(C4_up$row.names.genes.)
C5_U <- as.vector(C5_up$row.names.genes.)
C6_U <- as.vector(C6_up$row.names.genes.)
C7_U <- as.vector(C7_up$row.names.genes.)
C8_U <- as.vector(C8_up$row.names.genes.)
##For astrocyte disease-specific downregulated DEGs
C0_D <- as.vector(C0_down$row.names.genes.)
C1_D <- as.vector(C1_down$row.names.genes.)
C2_D <- as.vector(C2_down$row.names.genes.)
C3_D <- as.vector(C3_down$row.names.genes.)
C4_D <- as.vector(C4_down$row.names.genes.)
C5_D <- as.vector(C5_down$row.names.genes.)
C6_D <- as.vector(C6_down$row.names.genes.)
C7_D <- as.vector(C7_down$row.names.genes.)
C8_D <- as.vector(C8_down$row.names.genes.)

#Make a set list
##For astrocyte disease-specific upregulated DEGs
upregulated_read_sets = list(Cluster_0 = C0_U, Cluster_1 = C1_U, Cluster_2 = C2_U, Cluster_3 = C3_U,
Cluster_4 = C4_U, Cluster_5 = C5_U, Cluster_6 = C6_U, Cluster_7 = C7_U, Cluster_8 = C8_U)
##For astrocyte disease-specific downregulated DEGs
downregulated_read_sets = list(Cluster_0 = C0_D, Cluster_1 = C1_D, Cluster_2 = C2_D, Cluster_3 = C3_D,
Cluster_4 = C4_D, Cluster_5 = C5_D, Cluster_6 = C6_D, Cluster_7 = C7_D, Cluster_8 = C8_D)

#Visualize with UpSetR
##For astrocyte disease-specific upregulated DEGs
upset(fromList(upregulated_read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5", "Cluster_6", "Cluster_7", "Cluster_8"),order.by = "freq", keep.order = TRUE)
##For astrocyte disease-specific downregulated DEGs
upset(fromList(downregulated_read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5", "Cluster_6", "Cluster_7", "Cluster_8"),order.by = "freq", keep.order = TRUE)

#---------------------------------------------------------------------------------------------------
#GENERATE UpSetR PLOTS FOR ASTROCYTES DISEASE-SPECIFIC UPREGULATED GO TERMS BY CLUSTER
#Read in csv files for astrocyte disease-specific upregulated GO terms identified per cluster
Cluster0 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster0_upregulated_pathways.csv")
Cluster1 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster1_upregulated_pathways.csv")
Cluster2 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster2_upregulated_pathways.csv")
Cluster3 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster3_upregulated_pathways.csv")
Cluster4 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster4_upregulated_pathways.csv")
Cluster5 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster5_upregulated_pathways.csv")
Cluster6 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster6_upregulated_pathways.csv")
Cluster7 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster7_upregulated_pathways.csv")
Cluster8 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster8_upregulated_pathways.csv")

#Pull GO column as a vector
C0 <- as.vector(Cluster0$ID)
C1 <- as.vector(Cluster1$ID)
C2 <- as.vector(Cluster2$ID)
C3 <- as.vector(Cluster3$ID)
C4 <- as.vector(Cluster4$ID)
C5 <- as.vector(Cluster5$ID)
C6 <- as.vector(Cluster6$ID)
C7 <- as.vector(Cluster7$ID)
C8 <- as.vector(Cluster8$ID)

#Make a set list
read_sets = list(Cluster_0 = C0, Cluster_1 = C1, Cluster_2 = C2, Cluster_3 = C3,
Cluster_4 = C4, Cluster_5 = C5, Cluster_6 = C6, Cluster_7 = C7,
Cluster_8 = C8)

#Visualize with UpSetR
upset(fromList(read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5", "Cluster_6", "Cluster_7", "Cluster_8"),order.by = "freq", keep.order = TRUE)

#In order to identify overlapping GO terms
#ref: https://github.com/hms-dbmi/UpSetR/issues/85
fromList <- function (input) {
    #Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
}

upreg_pathways <- fromList(read_sets)

#Identify unique upregulated pathways per astrocyte cluster(s) by using logical operators
##Astrocyte cluster0
unique_C0 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 0 
                          & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                          & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                          & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster1
unique_C1 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                             & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster2
unique_C2 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 1 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                             & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster3
unique_C3 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                             & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                             & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster4
unique_C4 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                             & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster5
unique_C5 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 1
                             & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster6
unique_C6 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                             & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster46
unique_C46 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                             & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster16
unique_C16 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                              & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster146
unique_C146 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                              & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster13
unique_C13 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                              & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster14
unique_C14 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                              & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster016
unique_C016 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                              & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster06
unique_C06 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 0 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                               & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster134
unique_C134 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                               & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster01
unique_C01 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                               & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster0136
unique_C0136 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                              & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster34
unique_C34 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                                & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                                & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                                & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster346
unique_C346 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                              & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster45
unique_C45 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                               & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 1
                               & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster1346
unique_C1346 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                              & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster015
unique_C015 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                              & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 1
                              & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster0146
unique_C0146 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                               & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster013
unique_C013 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                               & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster01346
unique_C01346 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                               & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster136
unique_C136 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                                 & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                                 & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                                 & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster04
unique_C04 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 0 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                               & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster345
unique_C345 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                              & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 1
                              & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster3456
unique_C3456 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                               & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 1
                               & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster045
unique_C045 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 0 
                                & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                                & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 1
                                & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster03
unique_C03 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 0 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                               & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster036
unique_C036 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 0 
                              & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                              & upreg_pathways$Cluster_4 == 0 & upreg_pathways$Cluster_5 == 0
                              & upreg_pathways$Cluster_6 == 1), ]
##Astrocyte cluster0134
unique_C0134 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                               & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                               & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 0
                               & upreg_pathways$Cluster_6 == 0), ]
##Astrocyte cluster013456
unique_C013456 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 1 
                                & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                                & upreg_pathways$Cluster_4 == 1 & upreg_pathways$Cluster_5 == 1
                                & upreg_pathways$Cluster_6 == 1), ]

#Make a dataframe using rownames
C0 <- as.data.frame(row.names(unique_C0))
colnames(C0)<- c("GO")
C1 <- as.data.frame(row.names(unique_C1))
colnames(C1)<- c("GO")
C2 <- as.data.frame(row.names(unique_C2))
colnames(C2)<- c("GO")
C3 <- as.data.frame(row.names(unique_C3))
colnames(C3)<- c("GO")
C4 <- as.data.frame(row.names(unique_C4))
colnames(C4)<- c("GO")
C5 <- as.data.frame(row.names(unique_C5))
colnames(C5)<- c("GO")
C6 <- as.data.frame(row.names(unique_C6))
colnames(C6)<- c("GO")
C46 <- as.data.frame(row.names(unique_C46))
colnames(C46)<- c("GO")
C16 <- as.data.frame(row.names(unique_C16))
colnames(C16)<- c("GO")
C146 <- as.data.frame(row.names(unique_C146))
colnames(C146)<- c("GO")
C13 <- as.data.frame(row.names(unique_C13))
colnames(C13)<- c("GO")
C14 <- as.data.frame(row.names(unique_C14))
colnames(C14)<- c("GO")
C016 <- as.data.frame(row.names(unique_C016))
colnames(C016)<- c("GO")
C06 <- as.data.frame(row.names(unique_C06))
colnames(C06)<- c("GO")
C134 <- as.data.frame(row.names(unique_C134))
colnames(C134)<- c("GO")
C01 <- as.data.frame(row.names(unique_C01))
colnames(C01)<- c("GO")
C0136 <- as.data.frame(row.names(unique_C0136))
colnames(C0136)<- c("GO")
C34 <- as.data.frame(row.names(unique_C34))
colnames(C34)<- c("GO")
C346 <- as.data.frame(row.names(unique_C346))
colnames(C346)<- c("GO")
C45 <- as.data.frame(row.names(unique_C45))
colnames(C45)<- c("GO")
C1346 <- as.data.frame(row.names(unique_C1346))
colnames(C1346)<- c("GO")
C015 <- as.data.frame(row.names(unique_C015))
colnames(C015)<- c("GO")
C0146 <- as.data.frame(row.names(unique_C0146))
colnames(C0146)<- c("GO")
C013 <- as.data.frame(row.names(unique_C013))
colnames(C013)<- c("GO")
C01346 <- as.data.frame(row.names(unique_C01346))
colnames(C01346)<- c("GO")
C136 <- as.data.frame(row.names(unique_C136))
colnames(C136)<- c("GO")
C04 <- as.data.frame(row.names(unique_C04))
colnames(C04)<- c("GO")
C345 <- as.data.frame(row.names(unique_C345))
colnames(C345)<- c("GO")
C3456 <- as.data.frame(row.names(unique_C3456))
colnames(C3456)<- c("GO")
C045 <- as.data.frame(row.names(unique_C045))
colnames(C045)<- c("GO")
C03 <- as.data.frame(row.names(unique_C03))
colnames(C03)<- c("GO")
C036 <- as.data.frame(row.names(unique_C036))
colnames(C036)<- c("GO")
C0134 <- as.data.frame(row.names(unique_C0134))
colnames(C0134)<- c("GO")
C013456 <- as.data.frame(row.names(unique_C013456))
colnames(C013456)<- c("GO")

#Read in unique list of astrocyte disease-specific upregulated pathways
total_upreg_pathways <- read.csv(file = "./Astro_disease-specific_upregulated_pathways.csv")
colnames(total_upreg_pathways)<- c("GO", "Description")

#Subset pathways by astrocyte cluster
C0_subset <- inner_join(total_upreg_pathways, C0 %>% select(GO))
C1_subset <- inner_join(total_upreg_pathways, C1 %>% select(GO))
C2_subset <- inner_join(total_upreg_pathways, C2 %>% select(GO))
C3_subset <- inner_join(total_upreg_pathways, C3 %>% select(GO))
C4_subset <- inner_join(total_upreg_pathways, C4 %>% select(GO))
C5_subset <- inner_join(total_upreg_pathways, C5 %>% select(GO))
C6_subset <- inner_join(total_upreg_pathways, C6 %>% select(GO))
C46_subset <- inner_join(total_upreg_pathways, C46 %>% select(GO))
C16_subset <- inner_join(total_upreg_pathways, C16 %>% select(GO))
C146_subset <- inner_join(total_upreg_pathways, C146 %>% select(GO))
C13_subset <- inner_join(total_upreg_pathways, C13 %>% select(GO))
C14_subset <- inner_join(total_upreg_pathways, C14 %>% select(GO))
C016_subset <- inner_join(total_upreg_pathways, C016 %>% select(GO))
C06_subset <- inner_join(total_upreg_pathways, C06 %>% select(GO))
C134_subset <- inner_join(total_upreg_pathways, C134 %>% select(GO))
C01_subset <- inner_join(total_upreg_pathways, C01 %>% select(GO))
C0136_subset <- inner_join(total_upreg_pathways, C0136 %>% select(GO))
C34_subset <- inner_join(total_upreg_pathways, C34 %>% select(GO))
C346_subset <- inner_join(total_upreg_pathways, C346 %>% select(GO))
C45_subset <- inner_join(total_upreg_pathways, C45 %>% select(GO))
C1346_subset <- inner_join(total_upreg_pathways, C1346 %>% select(GO))
C015_subset <- inner_join(total_upreg_pathways, C015 %>% select(GO))
C0146_subset <- inner_join(total_upreg_pathways, C0146 %>% select(GO))
C013_subset <- inner_join(total_upreg_pathways, C013 %>% select(GO))
C01346_subset <- inner_join(total_upreg_pathways, C01346 %>% select(GO))
C136_subset <- inner_join(total_upreg_pathways, C136 %>% select(GO))
C04_subset <- inner_join(total_upreg_pathways, C04 %>% select(GO))
C345_subset <- inner_join(total_upreg_pathways, C345 %>% select(GO))
C3456_subset <- inner_join(total_upreg_pathways, C3456 %>% select(GO))
C045_subset <- inner_join(total_upreg_pathways, C045 %>% select(GO))
C03_subset <- inner_join(total_upreg_pathways, C03 %>% select(GO))
C036_subset <- inner_join(total_upreg_pathways, C036 %>% select(GO))
C0134_subset <- inner_join(total_upreg_pathways, C0134 %>% select(GO))
C013456_subset <- inner_join(total_upreg_pathways, C013456 %>% select(GO))

#Add astrocyte cluster(s) identifier to each subset
C0_subset$cluster <- "0"
C1_subset$cluster <- "1"
C2_subset$cluster <- "2"
C3_subset$cluster <- "3"
C4_subset$cluster <- "4"
C5_subset$cluster <- "5"
C6_subset$cluster <- "6"
C46_subset$cluster <- "4--6"
C16_subset$cluster <- "1--6"
C146_subset$cluster <- "1--4--6"
C13_subset$cluster <- "1--3"
C14_subset$cluster <- "1--4"
C016_subset$cluster <- "0--1--6"
C06_subset$cluster <- "0--6"
C134_subset$cluster <- "1--3--4"
C01_subset$cluster <- "0--1"
C0136_subset$cluster <- "0--1--3--6"
C34_subset$cluster <- "3--4"
C346_subset$cluster <- "3--4--6"
C45_subset$cluster <- "4--5"
C1346_subset$cluster <- "1--3--4--6"
C015_subset$cluster <- "0--1--5"
C0146_subset$cluster <- "0--1--4--6"
C013_subset$cluster <- "0--1--3"
C01346_subset$cluster <- "0--1--3--4--6"
C136_subset$cluster <- "1--3--6"
C04_subset$cluster <- "0--4"
C345_subset$cluster <- "3--4--5"
C3456_subset$cluster <- "3--4--5--6"
C045_subset$cluster <- "0--4--5"
C03_subset$cluster <- "0--3"
C036_subset$cluster <- "0--3--6"
C0134_subset$cluster <- "0--1--3--4"
C013456_subset$cluster <- "0--1--3--4--5--6"

#Join subsetted dataframes
unique <- rbind(C0_subset, C1_subset, C2_subset, C3_subset, C4_subset, C5_subset, C6_subset,
                C46_subset, C16_subset, C146_subset, C13_subset, C14_subset, C016_subset,
                C06_subset, C134_subset, C01_subset, C0136_subset, C34_subset, C346_subset, 
                C45_subset, C1346_subset, C015_subset, C0146_subset, C013_subset, 
                C01346_subset, C136_subset, C04_subset, C345_subset, C3456_subset, 
                C045_subset, C03_subset, C036_subset, C0134_subset, C013456_subset)
write.csv(unique, "R:/liddes01labspace/Jess/GRCh38_premRNA_CR4_objects/Astro_unique_disease-specific_upregulated_pathways_by_cluster.csv")

#---------------------------------------------------------------------------------------------------
#GENERATE UpSetR PLOTS FOR ASTROCYTES DISEASE-SPECIFIC DOWNREGULATED GO TERMS BY CLUSTER
#Read in csv files for astrocyte disease-specific downregulated GO terms identified per cluster
Cluster0 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster0_downregulated_pathways.csv")
Cluster1 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster1_downregulated_pathways.csv")
Cluster2 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster2_downregulated_pathways.csv")
Cluster3 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster3_downregulated_pathways.csv")
Cluster4 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster4_downregulated_pathways.csv")
Cluster5 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster5_downregulated_pathways.csv")
Cluster6 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster6_downregulated_pathways.csv")
Cluster7 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster7_downregulated_pathways.csv")
Cluster8 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_astro_r2_20PC_zinb_cluster8_downregulated_pathways.csv")

#Pull GO column as a vector
C0 <- as.vector(Cluster0$ID)
C1 <- as.vector(Cluster1$ID)
C2 <- as.vector(Cluster2$ID)
C3 <- as.vector(Cluster3$ID)
C4 <- as.vector(Cluster4$ID)
C5 <- as.vector(Cluster5$ID)
C6 <- as.vector(Cluster6$ID)
C7 <- as.vector(Cluster7$ID)
C8 <- as.vector(Cluster8$ID)

#Make a set list
read_sets = list(Cluster_1 = C1, Cluster_2 = C2, Cluster_3 = C3,
Cluster_4 = C4, Cluster_5 = C5, Cluster_6 = C6, Cluster_7 = C7,
Cluster_8 = C8, Cluster_0 = C0)

#Visualize with UpSetR
upset(fromList(read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5", "Cluster_6", "Cluster_7", "Cluster_8"),order.by = "freq", keep.order = TRUE)

#In order to identify overlapping GO terms
#ref: https://github.com/hms-dbmi/UpSetR/issues/85
fromList <- function (input) {
    # Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
}

downreg_pathways <- fromList(read_sets)

#Identify unique downregulated pathways per astrocyte cluster(s) by using logical operators
##Astrocyte cluster2
unique_C2 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                               & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                               & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                               & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster3
unique_C3 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                               & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                               & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                               & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster4
unique_C4 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                               & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                               & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                               & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster5
unique_C5 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                               & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                               & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                               & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster6
unique_C6 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                               & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                               & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                               & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster7
unique_C7 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                               & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                               & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                               & downreg_pathways$Cluster_7 == 1 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster8
unique_C8 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                               & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                               & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                               & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 1), ]
##Astrocyte cluster13
unique_C13 <- downreg_pathways[(downreg_pathways$Cluster_1 == 1 & downreg_pathways$Cluster_2 == 0 
                               & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                               & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                               & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster23
unique_C23 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster24
unique_C24 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster25
unique_C25 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster27
unique_C27 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 1 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster28
unique_C28 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 1), ]
##Astrocyte cluster34
unique_C34 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster35
unique_C35 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster36
unique_C36 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster45
unique_C45 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster46
unique_C46 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster56
unique_C56 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster123
unique_C123 <- downreg_pathways[(downreg_pathways$Cluster_1 == 1 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster234
unique_C234 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster235
unique_C235 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster236
unique_C236 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster245
unique_C245 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster246
unique_C246 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster248
unique_C248 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 1), ]
##Astrocyte cluster256
unique_C256 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster258
unique_C258 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 1), ]
##Astrocyte cluster346
unique_C346 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster456
unique_C456 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 0 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster1236
unique_C1236 <- downreg_pathways[(downreg_pathways$Cluster_1 == 1 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster2345
unique_C2345 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster2346
unique_C2346 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 0 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster2356
unique_C2356 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster2357
unique_C2357 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 0 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 0
                                & downreg_pathways$Cluster_7 == 1 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster2456
unique_C2456 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 0 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster23456
unique_C23456 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 0), ]
##Astrocyte cluster234568
unique_C234568 <- downreg_pathways[(downreg_pathways$Cluster_1 == 0 & downreg_pathways$Cluster_2 == 1 
                                & downreg_pathways$Cluster_3 == 1 & downreg_pathways$Cluster_4 == 1 
                                & downreg_pathways$Cluster_5 == 1 & downreg_pathways$Cluster_6 == 1
                                & downreg_pathways$Cluster_7 == 0 & downreg_pathways$Cluster_8 == 1), ]

#Make a dataframe using rownames
C2 <- as.data.frame(row.names(unique_C2))
colnames(C2)<- c("GO")
C3 <- as.data.frame(row.names(unique_C3))
colnames(C3)<- c("GO")
C4 <- as.data.frame(row.names(unique_C4))
colnames(C4)<- c("GO")
C5 <- as.data.frame(row.names(unique_C5))
colnames(C5)<- c("GO")
C6 <- as.data.frame(row.names(unique_C6))
colnames(C6)<- c("GO")
C7 <- as.data.frame(row.names(unique_C7))
colnames(C7)<- c("GO")
C8 <- as.data.frame(row.names(unique_C8))
colnames(C8)<- c("GO")
C13 <- as.data.frame(row.names(unique_C13))
colnames(C13)<- c("GO")
C23 <- as.data.frame(row.names(unique_C23))
colnames(C23)<- c("GO")
C24 <- as.data.frame(row.names(unique_C24))
colnames(C24)<- c("GO")
C25 <- as.data.frame(row.names(unique_C25))
colnames(C25)<- c("GO")
C27 <- as.data.frame(row.names(unique_C27))
colnames(C27)<- c("GO")
C28 <- as.data.frame(row.names(unique_C28))
colnames(C28)<- c("GO")
C34 <- as.data.frame(row.names(unique_C34))
colnames(C34)<- c("GO")
C35 <- as.data.frame(row.names(unique_C35))
colnames(C35)<- c("GO")
C36 <- as.data.frame(row.names(unique_C36))
colnames(C36)<- c("GO")
C45 <- as.data.frame(row.names(unique_C45))
colnames(C45)<- c("GO")
C46 <- as.data.frame(row.names(unique_C46))
colnames(C46)<- c("GO")
C56 <- as.data.frame(row.names(unique_C56))
colnames(C56)<- c("GO")
C123 <- as.data.frame(row.names(unique_C123))
colnames(C123)<- c("GO")
C234 <- as.data.frame(row.names(unique_C234))
colnames(C234)<- c("GO")
C235 <- as.data.frame(row.names(unique_C235))
colnames(C235)<- c("GO")
C236 <- as.data.frame(row.names(unique_C236))
colnames(C236)<- c("GO")
C245 <- as.data.frame(row.names(unique_C245))
colnames(C245)<- c("GO")
C246 <- as.data.frame(row.names(unique_C246))
colnames(C246)<- c("GO")
C248 <- as.data.frame(row.names(unique_C248))
colnames(C248)<- c("GO")
C256 <- as.data.frame(row.names(unique_C256))
colnames(C256)<- c("GO")
C258 <- as.data.frame(row.names(unique_C258))
colnames(C258)<- c("GO")
C346 <- as.data.frame(row.names(unique_C346))
colnames(C346)<- c("GO")
C456 <- as.data.frame(row.names(unique_C456))
colnames(C456)<- c("GO")
C1236 <- as.data.frame(row.names(unique_C1236))
colnames(C1236)<- c("GO")
C2345 <- as.data.frame(row.names(unique_C2345))
colnames(C2345)<- c("GO")
C2346 <- as.data.frame(row.names(unique_C2346))
colnames(C2346)<- c("GO")
C2356 <- as.data.frame(row.names(unique_C2356))
colnames(C2356)<- c("GO")
C2357 <- as.data.frame(row.names(unique_C2357))
colnames(C2357)<- c("GO")
C2456 <- as.data.frame(row.names(unique_C2456))
colnames(C2456)<- c("GO")
C23456 <- as.data.frame(row.names(unique_C23456))
colnames(C23456)<- c("GO")
C234568 <- as.data.frame(row.names(unique_C234568))
colnames(C234568)<- c("GO")

#Read in unique list of astrocyte disease-specific downregulated pathways
total_downreg_pathways <- read.csv(file = "./Astro_disease-specific_downregulated_pathways.csv")
colnames(total_downreg_pathways)<- c("GO", "Description")

#Subset total pathways by astrocyte cluster
C2_subset <- inner_join(total_downreg_pathways, C2 %>% select(GO))
C3_subset <- inner_join(total_downreg_pathways, C3 %>% select(GO))
C4_subset <- inner_join(total_downreg_pathways, C4 %>% select(GO))
C5_subset <- inner_join(total_downreg_pathways, C5 %>% select(GO))
C6_subset <- inner_join(total_downreg_pathways, C6 %>% select(GO))
C7_subset <- inner_join(total_downreg_pathways, C7 %>% select(GO))
C8_subset <- inner_join(total_downreg_pathways, C8 %>% select(GO))
C13_subset <- inner_join(total_downreg_pathways, C13 %>% select(GO))
C23_subset <- inner_join(total_downreg_pathways, C23 %>% select(GO))
C24_subset <- inner_join(total_downreg_pathways, C24 %>% select(GO))
C25_subset <- inner_join(total_downreg_pathways, C25 %>% select(GO))
C27_subset <- inner_join(total_downreg_pathways, C27 %>% select(GO))
C28_subset <- inner_join(total_downreg_pathways, C28 %>% select(GO))
C34_subset <- inner_join(total_downreg_pathways, C34 %>% select(GO))
C35_subset <- inner_join(total_downreg_pathways, C35 %>% select(GO))
C36_subset <- inner_join(total_downreg_pathways, C36 %>% select(GO))
C45_subset <- inner_join(total_downreg_pathways, C45 %>% select(GO))
C46_subset <- inner_join(total_downreg_pathways, C46 %>% select(GO))
C56_subset <- inner_join(total_downreg_pathways, C56 %>% select(GO))
C123_subset <- inner_join(total_downreg_pathways, C123 %>% select(GO))
C234_subset <- inner_join(total_downreg_pathways, C234 %>% select(GO))
C235_subset <- inner_join(total_downreg_pathways, C235 %>% select(GO))
C236_subset <- inner_join(total_downreg_pathways, C236 %>% select(GO))
C245_subset <- inner_join(total_downreg_pathways, C245 %>% select(GO))
C246_subset <- inner_join(total_downreg_pathways, C246 %>% select(GO))
C248_subset <- inner_join(total_downreg_pathways, C248 %>% select(GO))
C256_subset <- inner_join(total_downreg_pathways, C256 %>% select(GO))
C258_subset <- inner_join(total_downreg_pathways, C258 %>% select(GO))
C346_subset <- inner_join(total_downreg_pathways, C346 %>% select(GO))
C456_subset <- inner_join(total_downreg_pathways, C456 %>% select(GO))
C1236_subset <- inner_join(total_downreg_pathways, C1236 %>% select(GO))
C2345_subset <- inner_join(total_downreg_pathways, C2345 %>% select(GO))
C2346_subset <- inner_join(total_downreg_pathways, C2346 %>% select(GO))
C2356_subset <- inner_join(total_downreg_pathways, C2356 %>% select(GO))
C2357_subset <- inner_join(total_downreg_pathways, C2357 %>% select(GO))
C2456_subset <- inner_join(total_downreg_pathways, C2456 %>% select(GO))
C23456_subset <- inner_join(total_downreg_pathways, C23456 %>% select(GO))
C234568_subset <- inner_join(total_downreg_pathways, C234568 %>% select(GO))

#Add astrocyte cluster(s) identifier to each subset
C2_subset$cluster <- "2"
C3_subset$cluster <- "3"
C4_subset$cluster <- "4"
C5_subset$cluster <- "5"
C6_subset$cluster <- "6"
C7_subset$cluster <- "7"
C8_subset$cluster <- "8"
C13_subset$cluster <- "1--3"
C23_subset$cluster <- "2--3"
C24_subset$cluster <- "2--4"
C25_subset$cluster <- "2--5"
C27_subset$cluster <- "2--7"
C28_subset$cluster <- "2--8"
C34_subset$cluster <- "3--4"
C35_subset$cluster <- "3--5"
C36_subset$cluster <- "3--6"
C45_subset$cluster <- "4--5"
C46_subset$cluster <- "4--6"
C56_subset$cluster <- "5--6"
C123_subset$cluster <- "1--2--3"
C234_subset$cluster <- "2--3--4"
C235_subset$cluster <- "2--3--5"
C236_subset$cluster <- "2--3--6"
C245_subset$cluster <- "2--4--5"
C246_subset$cluster <- "2--4--6"
C248_subset$cluster <- "2--4--8"
C256_subset$cluster <- "2--5--6"
C258_subset$cluster <- "2--5--8"
C346_subset$cluster <- "3--4--6"
C456_subset$cluster <- "4--5--6"
C1236_subset$cluster <- "1--2--3--6"
C2345_subset$cluster <- "2--3--4--5"
C2346_subset$cluster <- "2--3--4--6"
C2356_subset$cluster <- "2--3--5--6"
C2357_subset$cluster <- "2--3--5--7"
C2456_subset$cluster <- "2--4--5--6"
C23456_subset$cluster <- "2--3--4--5--6"
C234568_subset$cluster <- "2--3--4--5--6--8"

#Join subsetted dataframes
unique <- rbind(C2_subset, C3_subset, C4_subset, C5_subset, C6_subset, C7_subset, C8_subset,
                C13_subset, C23_subset, C24_subset, C25_subset, C27_subset, C28_subset, C34_subset, 
                C35_subset, C36_subset, C45_subset, C46_subset, C56_subset, C123_subset, C234_subset, 
                C235_subset, C236_subset, C245_subset, C246_subset, C248_subset, C256_subset, 
                C258_subset, C346_subset, C456_subset, C1236_subset, C2345_subset, C2346_subset,
                C2356_subset, C2357_subset, C2456_subset, C23456_subset, C234568_subset)
write.csv(unique, "R:/liddes01labspace/Jess/GRCh38_premRNA_CR4_objects/Astro_unique_disease-specific_downregulated_pathways_by_cluster.csv")

#-----------------------------------------------------------------------------------------------------
#GENERATE UpSetR PLOTS FOR OLIGODENDROCYTE DISEASE-SPECIFIC DEGs BY CLUSTER
#Read in csv files for oligodendrocyte disease-specific DEGs identified per cluster
Cluster0 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster0-0.25lfc.csv")
Cluster1 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster1-0.25lfc.csv")
Cluster2 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster2-0.25lfc.csv")
Cluster3 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster3-0.25lfc.csv")
Cluster4 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_sig_genes_cluster4-0.25lfc.csv")

##Subset oligodendrocyte disease-specific upregulated DEGs per cluster
C0_up <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC >= 0.25)
C1_up <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC >= 0.25)
C2_up <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC >= 0.25)
C3_up <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC >= 0.25)
C4_up <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC >= 0.25)
##Subset oligodendrocyte disease-specific downregulated DEGs per cluster
C0_down <- subset(Cluster0, select = -c(X,Pval), Cluster0$logFC <= -0.25)
C1_down <- subset(Cluster1, select = -c(X,Pval), Cluster1$logFC <= -0.25)
C2_down <- subset(Cluster2, select = -c(X,Pval), Cluster2$logFC <= -0.25)
C3_down <- subset(Cluster3, select = -c(X,Pval), Cluster3$logFC <= -0.25)
C4_down <- subset(Cluster4, select = -c(X,Pval), Cluster4$logFC <= -0.25)

#Subset row names as a vector
##For oligodendrocyte disease-specific upregulated DEGs
C0_U <- as.vector(C0_up$row.names.genes.)
C1_U <- as.vector(C1_up$row.names.genes.)
C2_U <- as.vector(C2_up$row.names.genes.)
C3_U <- as.vector(C3_up$row.names.genes.)
C4_U <- as.vector(C4_up$row.names.genes.)
##For oligodendrocyte disease-specific downregulated DEGs
C0_D <- as.vector(C0_down$row.names.genes.)
C1_D <- as.vector(C1_down$row.names.genes.)
C2_D <- as.vector(C2_down$row.names.genes.)
C3_D <- as.vector(C3_down$row.names.genes.)
C4_D <- as.vector(C4_down$row.names.genes.)

#Make a set list
##For oligodendrocyte disease-specific upregulated DEGs
upregulated_read_sets = list(Cluster_0 = C0_U, Cluster_1 = C1_U, Cluster_2 = C2_U, Cluster_3 = C3_U,
Cluster_4 = C4_U)
##For oligodendrocyte disease-specific downregulated DEGs
downregulated_read_sets = list(Cluster_0 = C0_D, Cluster_1 = C1_D, Cluster_2 = C2_D, Cluster_3 = C3_D,
Cluster_4 = C4_D)

#Visualize with UpSetR
##For oligodendrocyte disease-specific upregulated DEGs
upset(fromList(upregulated_read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"),order.by = "freq", keep.order = TRUE)
##For oligodendrocyte disease-specific downregulated DEGs
upset(fromList(downregulated_read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"),order.by = "freq", keep.order = TRUE)

#---------------------------------------------------------------------------------------------------
#GENERATE UpSetR PLOTS FOR OLIGODENDROCYTES DISEASE-SPECIFIC UPREGULATED GO TERMS BY CLUSTER
#Read in csv files for oligodendrocyte disease-specific upregulated GO terms identified per cluster
Cluster0 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster0_upregulated_pathways.csv")
Cluster1 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster1_upregulated_pathways.csv")
Cluster2 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster2_upregulated_pathways.csv")
Cluster3 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster3_upregulated_pathways.csv")
Cluster4 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster4_upregulated_pathways.csv")

#Pull GO column as a vector
C0 <- as.vector(Cluster0$ID)
C1 <- as.vector(Cluster1$ID)
C2 <- as.vector(Cluster2$ID)
C3 <- as.vector(Cluster3$ID)
C4 <- as.vector(Cluster4$ID)

#Make a set list
read_sets = list(Cluster_1 = C1, Cluster_0 = C0, Cluster_2 = C2, Cluster_3 = C3,
Cluster_4 = C4)

#Visualize with UpSetR
upset(fromList(read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"),order.by = "freq", keep.order = TRUE)

#In order to identify overlapping GO terms
#ref: https://github.com/hms-dbmi/UpSetR/issues/85
fromList <- function (input) {
    # Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
}

upreg_pathways <- fromList(read_sets)

#Identify unique upregulated pathways per oligodendrocyte cluster(s) by using logical operators
##Oligodendrocyte cluster1
unique_C1 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster2
unique_C2 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 1 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster3
unique_C3 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 1
                             & upreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster4
unique_C4 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 1), ]
##Oligodendrocyte cluster02
unique_C02 <- upreg_pathways[(upreg_pathways$Cluster_0 == 1 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 1 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster12
unique_C12 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                             & upreg_pathways$Cluster_2 == 1 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster14
unique_C14 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 1 
                             & upreg_pathways$Cluster_2 == 0 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 1), ]
##Oligodendrocyte cluster24
unique_C24 <- upreg_pathways[(upreg_pathways$Cluster_0 == 0 & upreg_pathways$Cluster_1 == 0 
                             & upreg_pathways$Cluster_2 == 1 & upreg_pathways$Cluster_3 == 0
                             & upreg_pathways$Cluster_4 == 1), ]

#Make a dataframe using rownames
C1 <- as.data.frame(row.names(unique_C1))
colnames(C1)<- c("GO")
C2 <- as.data.frame(row.names(unique_C2))
colnames(C2)<- c("GO")
C3 <- as.data.frame(row.names(unique_C3))
colnames(C3)<- c("GO")
C4 <- as.data.frame(row.names(unique_C4))
colnames(C4)<- c("GO")
C02 <- as.data.frame(row.names(unique_C02))
colnames(C02)<- c("GO")
C12 <- as.data.frame(row.names(unique_C12))
colnames(C12)<- c("GO")
C14 <- as.data.frame(row.names(unique_C14))
colnames(C14)<- c("GO")
C24 <- as.data.frame(row.names(unique_C24))
colnames(C24)<- c("GO")

#Read in unique list of oligodendrocyte disease-specific upregulated pathways
total_upreg_pathways <- read.csv(file = "./Oligo_disease-specific_upregulated_pathways.csv")
colnames(total_upreg_pathways)<- c("GO", "Description")

#Subset total pathways by oligodendrocyte cluster
C1_subset <- inner_join(total_upreg_pathways, C1 %>% select(GO))
C2_subset <- inner_join(total_upreg_pathways, C2 %>% select(GO))
C3_subset <- inner_join(total_upreg_pathways, C3 %>% select(GO))
C4_subset <- inner_join(total_upreg_pathways, C4 %>% select(GO))
C02_subset <- inner_join(total_upreg_pathways, C02 %>% select(GO))
C12_subset <- inner_join(total_upreg_pathways, C12 %>% select(GO))
C14_subset <- inner_join(total_upreg_pathways, C14 %>% select(GO))
C24_subset <- inner_join(total_upreg_pathways, C24 %>% select(GO))

#Add oligodendrocyte cluster(s) identifier to each subset
C1_subset$cluster <- "1"
C2_subset$cluster <- "2"
C3_subset$cluster <- "3"
C4_subset$cluster <- "4"
C02_subset$cluster <- "0--2"
C12_subset$cluster <- "1--2"
C14_subset$cluster <- "1--4"
C24_subset$cluster <- "2--4"

#Join subsetted dataframes
unique <- rbind(C1_subset, C2_subset, C3_subset, C4_subset, C02_subset, C12_subset, C14_subset,
                C24_subset)
write.csv(unique, "R:/liddes01labspace/Jess/GRCh38_premRNA_CR4_objects/Oligo_unique_disease-specific_upregulated_pathways_by_cluster.csv")

#---------------------------------------------------------------------------------------------------
#GENERATE UpSetR PLOTS FOR OLIGODENDROCYTES DISEASE-SPECIFIC DOWNREGULATED GO TERMS BY CLUSTER
#Read in csv files for oligodendrocyte disease-specific downregulated GO terms identified per cluster
Cluster0 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster0_downregulated_pathways.csv")
Cluster1 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster1_downregulated_pathways.csv")
Cluster2 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster2_downregulated_pathways.csv")
Cluster3 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster3_downregulated_pathways.csv")
Cluster4 <- read.csv(file = "./CR4_LEN_e23_noD5D9_so_oligo_15PC_zinb_cluster4_downregulated_pathways.csv")

#Pull GO column as a vector
C0 <- as.vector(Cluster0$ID)
C1 <- as.vector(Cluster1$ID)
C2 <- as.vector(Cluster2$ID)
C3 <- as.vector(Cluster3$ID)
C4 <- as.vector(Cluster4$ID)

#Make a set list
read_sets = list( Cluster_0 = C0, Cluster_1 = C1, Cluster_2 = C2, Cluster_3 = C3,
Cluster_4 = C4)

#Visualize with UpSetR
upset(fromList(read_sets),sets = c("Cluster_0", "Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4"),order.by = "freq", keep.order = TRUE)

#In order to identify overlapping GO terms
#ref: https://github.com/hms-dbmi/UpSetR/issues/85
fromList <- function (input) {
    # Same as original fromList()...
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    # ... Except now it conserves your original value names!
    row.names(data) <- elements
    return(data)
}

downreg_pathways <- fromList(read_sets)

#Identify unique downregulated pathways per oligodendrocyte cluster(s) by using logical operators
##Oligodendrocyte cluster0
unique_C0 <- downreg_pathways[(downreg_pathways$Cluster_0 == 1 & downreg_pathways$Cluster_1 == 0 
                               & downreg_pathways$Cluster_2 == 0 & downreg_pathways$Cluster_3 == 0 
                               & downreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster3
unique_C3 <- downreg_pathways[(downreg_pathways$Cluster_0 == 0 & downreg_pathways$Cluster_1 == 0 
                               & downreg_pathways$Cluster_2 == 0 & downreg_pathways$Cluster_3 == 1 
                               & downreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster4
unique_C4 <- downreg_pathways[(downreg_pathways$Cluster_0 == 0 & downreg_pathways$Cluster_1 == 0 
                               & downreg_pathways$Cluster_2 == 0 & downreg_pathways$Cluster_3 == 0 
                               & downreg_pathways$Cluster_4 == 1), ]
##Oligodendrocyte cluster02
unique_C02 <- downreg_pathways[(downreg_pathways$Cluster_0 == 1 & downreg_pathways$Cluster_1 == 0 
                               & downreg_pathways$Cluster_2 == 1 & downreg_pathways$Cluster_3 == 0 
                               & downreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster03
unique_C03 <- downreg_pathways[(downreg_pathways$Cluster_0 == 1 & downreg_pathways$Cluster_1 == 0 
                               & downreg_pathways$Cluster_2 == 0 & downreg_pathways$Cluster_3 == 1 
                               & downreg_pathways$Cluster_4 == 0), ]
##Oligodendrocyte cluster04
unique_C04 <- downreg_pathways[(downreg_pathways$Cluster_0 == 1 & downreg_pathways$Cluster_1 == 0 
                               & downreg_pathways$Cluster_2 == 0 & downreg_pathways$Cluster_3 == 0 
                               & downreg_pathways$Cluster_4 == 1), ]

#Make a dataframe using rownames
C0 <- as.data.frame(row.names(unique_C0))
colnames(C0)<- c("GO")
C3 <- as.data.frame(row.names(unique_C3))
colnames(C3)<- c("GO")
C4 <- as.data.frame(row.names(unique_C4))
colnames(C4)<- c("GO")
C02 <- as.data.frame(row.names(unique_C02))
colnames(C02)<- c("GO")
C03 <- as.data.frame(row.names(unique_C03))
colnames(C03)<- c("GO")
C04 <- as.data.frame(row.names(unique_C04))
colnames(C04)<- c("GO")

#Read in unique list of oligodendrocyte disease-specific downregulated pathways
total_downreg_pathways <- read.csv(file = "./Oligo_disease-specific_downregulated_pathways.csv")
colnames(total_downreg_pathways)<- c("GO", "Description")

#Subset total pathways by oligodendrocyte cluster
C0_subset <- inner_join(total_downreg_pathways, C0 %>% select(GO))
C3_subset <- inner_join(total_downreg_pathways, C3 %>% select(GO))
C4_subset <- inner_join(total_downreg_pathways, C4 %>% select(GO))
C02_subset <- inner_join(total_downreg_pathways, C02 %>% select(GO))
C03_subset <- inner_join(total_downreg_pathways, C03 %>% select(GO))
C04_subset <- inner_join(total_downreg_pathways, C04 %>% select(GO))

#Add oligodendrocyte cluster(s) identifier to each subset
C0_subset$cluster <- "0"
C3_subset$cluster <- "3"
C4_subset$cluster <- "4"
C02_subset$cluster <- "0--2"
C03_subset$cluster <- "0--3"
C04_subset$cluster <- "0--4"

#Join subsetted dataframes
unique <- rbind(C0_subset, C3_subset, C4_subset, C02_subset, C03_subset, C04_subset)
write.csv(unique, "R:/liddes01labspace/Jess/GRCh38_premRNA_CR4_objects/Oligo_unique_disease-specific_downregulated_pathways_by_cluster.csv")
