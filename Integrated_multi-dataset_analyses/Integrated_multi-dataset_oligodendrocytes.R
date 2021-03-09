#Integrated multi-dataset analysis of oligodendrocytes
#Goal: To directly compare oligodendrocytes captured in multiple AD snRNA-seq datasets, including Mathys et al. (2019), Grubman et al. (2019), Zhou et al. (2020), and our snRNA-seq (LHX2-positive/NeuN-negative (LEN) sorted nuclei) dataset.
#Pipeline prepared by Jessica S. Sadick
#Majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#---------------------------------------------------------------------------------------------------
#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load SeuratObjects
so_Mathys_oligo <- readRDS(file.path("file_path", "CR4_Mathys_so_oligo-ref_12PC.rds"))
so_Grubman_oligo <- readRDS(file.path("file_path", "CR4_Grubman_so_oligo_12PC.rds"))
so_Zhou_oligo <- readRDS(file.path("file_path", "CR4_Zhou_so_oligo-ref_15PC.rds"))
so_LEN_oligo <- readRDS(file.path("file_path", "CR4_LEN_e23_noD5D9_so_oligo_15PC.rds"))

#Add in metadata dataset identifiers
Meta_Mathys <- so_Mathys_oligo@meta.data
Meta_Mathys["dataset"] <- c("Mathys")
NewMeta_Mathys <- subset(Meta_Mathys, select = c("dataset"))
so_Mathys_oligo <- AddMetaData(so_Mathys_oligo, NewMeta_Mathys)
head(x = so_Mathys_oligo[[]])

Meta_Grubman <- so_Grubman_oligo@meta.data
Meta_Grubman["dataset"] <- c("Grubman")
NewMeta_Grubman <- subset(Meta_Grubman, select = c("dataset"))
so_Grubman_oligo <- AddMetaData(so_Grubman_oligo, NewMeta_Grubman)
head(x = so_Grubman_oligo[[]])

Meta_Zhou <- so_Zhou_oligo@meta.data
Meta_Zhou["dataset"] <- c("Zhou")
NewMeta_Zhou <- subset(Meta_Zhou, select = c("dataset"))
so_Zhou_oligo <- AddMetaData(so_Zhou_oligo, NewMeta_Zhou)
head(x = so_Zhou_oligo[[]])

Meta_Sadick <- so_LEN_oligo@meta.data
Meta_Sadick["dataset"] <- c("Sadick")
NewMeta_Sadick <- subset(Meta_Sadick, select = c("dataset"))
so_LEN_oligo <- AddMetaData(so_LEN_oligo, NewMeta_Sadick)
head(x = so_LEN_oligo[[]])

#Merge seurat objects
so_oligo_merge <- merge(x = so_LEN_oligo, y = c(so_Mathys_oligo, so_Grubman_oligo, so_Zhou_oligo))

#Remove integrated data associated with object and shift assay into RNA
DefaultAssay(so_oligo_merge) <- "RNA"
so_oligo_merge[['integrated']] <- NULL 
saveRDS(so_oligo_merge, file.path("file_path", "CR4_so_oligo_merge_all_unint.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTERING

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load SeuratObject
so_oligo_merge <- readRDS(file.path("file_path", "CR4_so_oligo_merge_all_unint.rds"))

#Prepare SCE object from merged object
sce <- as.SingleCellExperiment(so_oligo_merge, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_oligo_merge <- lapply(cells_by_sample, function(i)
  SubsetData(so_oligo_merge, cells = i))

#Normalize, find variable genes, and scale
so_oligo_merge <- lapply(so_oligo_merge, NormalizeData, verbose = FALSE)
so_oligo_merge <- lapply(so_oligo_merge, FindVariableFeatures, nfeatures = 2e3,
                         selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_oligo_merge <- lapply(so_oligo_merge, ScaleData, verbose = FALSE)

#Find anchors and integrate
##HAD TO DECREASE k.filter PARAMETER (minimum oligo nuclei/donor = 65)
#Used reference based integration because memory constraints
#Chose 1 donor from each dataset AND each condition with the most nuclei captured

#Mathys
D178757 <- which(names(so_oligo_merge) == "D17-8757")
D178759 <- which(names(so_oligo_merge) == "D17-8759")
D178780 <- which(names(so_oligo_merge) == "D17-8780")
D178768 <- which(names(so_oligo_merge) == "D17-8768")

#Grubman
AD1_AD2 <- which(names(so_oligo_merge) == "AD1_AD2")
AD3_AD4 <- which(names(so_oligo_merge) == "AD3_AD4")
Ct1_Ct2 <- which(names(so_oligo_merge) == "Ct1_Ct2")
Ct3_Ct4 <- which(names(so_oligo_merge) == "Ct3_Ct4")

#Zhou
AD5 <- which(names(so_oligo_merge) == "AD5")
AD12 <- which(names(so_oligo_merge) == "AD12")
C8 <- which(names(so_oligo_merge) == "C8")
C11 <- which(names(so_oligo_merge) == "C11")
P6 <- which(names(so_oligo_merge) == "P6")
P9 <- which(names(so_oligo_merge) == "P9")

#Sadick
D3 <- which(names(so_oligo_merge) == "D3")
D6 <- which(names(so_oligo_merge) == "D6")
D10 <- which(names(so_oligo_merge) == "D10")
D15 <- which(names(so_oligo_merge) == "D15")

as <- FindIntegrationAnchors(so_oligo_merge, reference = c(D178757, D178759, D178780, D178768, AD1_AD2, AD3_AD4, Ct1_Ct2, Ct3_Ct4, AD5, AD12, C8, C11, P6, P9, D3, D6, D10, D15), verbose=TRUE, k.filter = 65)
so_oligo_merge <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_oligo_merge) <- "integrated"
so_oligo_merge <- ScaleData(so_oligo_merge, display.progress = FALSE)

#DIMENSION REDUCTION
so_oligo_merge <- RunPCA(so_oligo_merge, npcs = 50, verbose = FALSE)
ElbowPlot(so_oligo_merge, ndims = 50)

so_oligo_merge <- RunTSNE(so_oligo_merge, reduction = "pca", dims = seq_len(12),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE, check_duplicates = FALSE)
so_oligo_merge <- RunUMAP(so_oligo_merge, reduction = "pca", dims = seq_len(12),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_oligo_merge <- FindNeighbors(so_oligo_merge, reduction = "pca", dims = seq_len(12), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_oligo_merge <- FindClusters(so_oligo_merge, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_oligo_merge, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_oligo_merge, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#Save SeuratObject
saveRDS(so_oligo_merge, file.path("file_path", "CR4_so_oligo_merge_all-ref_MGZS_12PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION

#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(inauguration)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so_oligo_merge <- readRDS(file.path("file_path", "CR4_so_oligo_merge_all-ref_MGZS_12PC.rds"))
sce <- as.SingleCellExperiment(so_oligo_merge, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_oligo_merge <- SetIdent(so_oligo_merge, value = "integrated_snn_res.0.3")
so_oligo_merge@meta.data$cluster_id <- Idents(so_oligo_merge)
sce$cluster_id <- Idents(so_oligo_merge)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/so_oligo_merge_all-ref_MGZS_12PC_res0.3_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green", "yellow", "purple"), gids)

#Generate relative cluster abundances
fqs <- prop.table(n_cells, margin = 2)
mat <- as.matrix(unclass(fqs))
Heatmap(mat,
        col = rev(brewer.pal(11, "RdGy")[-6]),
        name = "Frequency",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_title = "cluster_id",
        column_title = "sample_id",
        column_title_side = "bottom",
        rect_gp = gpar(col = "white"),
        cell_fun = function(i, j, x, y, width, height, fill)
          grid.text(round(mat[j, i] * 100, 2), x = x, y = y, 
                    gp = gpar(col = "white", fontsize = 8)))

#DR colored by cluster ID
cs <- sample(colnames(so_oligo_merge), 5e3)
.plot_dr <- function(so_oligo_merge, dr, id)
  DimPlot(so_oligo_merge, cells = cs, group.by = id, reduction = dr, pt.size = 0.1) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_oligo_merge, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_oligo_merge, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^MT-", x = rownames(so_oligo_merge@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_oligo_merge@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_oligo_merge@assays[["RNA"]])
so_oligo_merge$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^RP[SL]", x = rownames(so_oligo_merge@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so_oligo_merge@assays[["RNA"]][rb.genes, ])/Matrix::colSums(so_oligo_merge@assays[["RNA"]])
so_oligo_merge$percent.rb <- percent.rb

VlnPlot(object = so_oligo_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0, cols = (usecol("pal_unikn_pair", 7)))

VlnPlot(object = so_oligo_merge, features = c("nFeature_RNA"), pt.size = 0) + stat_summary(fun.y=median, geom="point", shape=23, size=2)

so_oligo_merge_copy <- so_oligo_merge
Idents(so_oligo_merge_copy) <- "dataset"
designated_levels <- c("Mathys", "Grubman", "Zhou", "Sadick")
Idents(so_oligo_merge_copy) <- factor(Idents(so_oligo_merge_copy), levels= designated_levels)
DefaultAssay(so_oligo_merge_copy) <- "RNA"
VlnPlot(object = so_oligo_merge_copy, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0, cols = inauguration("inauguration_2021", 4))

#Generate summary statistics for entire dataset
summary(so_oligo_merge$nFeature_RNA)
summary(so_oligo_merge$nCount_RNA)

library(pastecs)
stat.desc(so_oligo_merge$nFeature_RNA)
stat.desc(so_oligo_merge$nCount_RNA)

library(psych)
describe(so_oligo_merge$nFeature_RNA)
describe(so_oligo_merge$nCount_RNA)

#Generate summary statistics per sample
library(data.table)
library(psych)
feature_by_sample <- as.data.frame(so_oligo_merge$nFeature_RNA, row.names = so_oligo_merge$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so_oligo_merge$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "file_path/so_oligo_merge_all-ref_MGZS_12PC_res0.3_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so_oligo_merge$nCount_RNA, row.names = so_oligo_merge$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so_oligo_merge$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/so_oligo_merge_all-ref_MGZS_12PC_res0.3_QC_count_by_sample.csv")

#Generate ssummary statistics per cluster
library(data.table)
library(psych)
feature_by_cluster <- as.data.frame(so_oligo_merge$nFeature_RNA, row.names = so_oligo_merge$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so_oligo_merge$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/so_oligo_merge_all-ref_MGZS_12PC_res0.3_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so_oligo_merge$nCount_RNA, row.names = so_oligo_merge$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so_oligo_merge$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/so_oligo_merge_all-ref_MGZS_12PC_res0.3_count_by_cluster.csv")

#Assess oligo clusters
DefaultAssay(so_oligo_merge) <- "RNA"
DimPlot(so_oligo_merge, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair"))

#Find all markers
DefaultAssay(so_oligo_merge) <- "RNA"
so_oligo_merge.markers <- FindAllMarkers(so_oligo_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_oligo_merge.markers, "file_path/so_oligo_merge_all-ref_MGZS_12PC_res0.3_genes-RNA.csv")

#Visualize
VlnPlot(so_oligo_merge, features = c("ALDOC"), pt.size = 0)
DotPlot(so_oligo_merge, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(so_oligo_merge, features = c("PLP1", "MBP", "MOG", "PDGFRA", "GFAP", "AQP4"), reduction = "tsne")

#Evaluate meta-variables
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "dataset", pt.size = 0.001, cols = rev(inauguration("inauguration_2021", 4))) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "dataset", pt.size = 0.001, cols = rev(wes_palette(4, name = "Zissou1", type = "continuous"))) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "dataset", pt.size = 0.001, cols = rev(wes_palette(4, name = "Darjeeling1", type = "continuous"))) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "sample_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "sex_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "disease_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo_merge, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_oligo_merge) <- "RNA"
so_oligo_merge_cluster.averages <- AverageExpression(so_oligo_merge)
head(so_oligo_merge_cluster.averages[["RNA"]][, 1:5])
write.csv(so_oligo_merge_cluster.averages[["RNA"]], "file_path/so_oligo_merge_all-ref_MGZS_12PC_res0.3_avg_exp_by_cluster.csv")

orig.levels <- levels(so_oligo_merge)
Idents(so_oligo_merge) <- gsub(pattern = " ", replacement = "_", x = Idents(so_oligo_merge))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_oligo_merge) <- orig.levels
so_oligo_merge_cluster.averages <- AverageExpression(so_oligo_merge, return.seurat = TRUE)
so_oligo_merge_cluster.averages

CellScatter(so_oligo_merge_cluster.averages, cell1 = "2", cell2 = "3")

#Plot 'in silico' bulk datasets
DoHeatmap(so_oligo_merge_cluster.averages, features = unlist(TopFeatures(so_oligo_merge[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)

DoHeatmap(so_oligo_merge_cluster.averages, features = c("SVEP1","LINC01608","PLXDC2","DYSF","CTNNA2","CNDP1","ST3GAL6","QDPR","CRYAB","FP236383.3","MT-ND4","MT-ND3","MT-CO2","MT-ATP6","ACTN2","SLC5A11","RASGRF1","LINC00609","ANKRD18A","SGCZ","MDGA2","CNTN1","KCNIP4","FRY","RBFOX1","AFF3","ACSBG1","COL18A1","NRP2","LUCAT1","NAV2","CAMK2D","NEAT1"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))

##to visualize clusterXdataset features
DefaultAssay(so_oligo_merge) <- "RNA"
so_oligo_merge$cluster.dataset <- paste(Idents(so_oligo_merge), so_oligo_merge$dataset, sep = "_")
so_oligo_merge$cluster <- Idents(so_oligo_merge)
Idents(so_oligo_merge) <- "cluster.dataset"
designated_levels <- c("0_Mathys", "0_Grubman", "0_Zhou", "0_Sadick",
                       "1_Mathys", "1_Grubman", "1_Zhou", "1_Sadick",
                       "2_Mathys", "2_Grubman", "2_Zhou", "2_Sadick",
                       "3_Mathys", "3_Grubman", "3_Zhou", "3_Sadick",
                       "4_Mathys", "4_Grubman", "4_Zhou", "4_Sadick",
                       "5_Mathys", "5_Grubman", "5_Zhou", "5_Sadick",
                       "6_Mathys", "6_Grubman", "6_Zhou", "6_Sadick")
Idents(so_oligo_merge) <- factor(Idents(so_oligo_merge), levels= designated_levels)

DefaultAssay(so_oligo_merge) <- "RNA"
DimPlot(so_oligo_merge, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1)
