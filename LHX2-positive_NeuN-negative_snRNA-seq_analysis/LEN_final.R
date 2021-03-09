#LIM Homeobox 2 (LHX2)-positive/NeuN-negative snRNA-seq analysis for *FINAL* donors characterized
#Goal: To evaluate astrocyte and oligodendrocyte heterogeneity in APOE e2/3 Alzheimer's disease and age-matched non-symptomatic patient final cohort (i.e., 2 outlier donors removed).
#Pipeline prepared by Jessica S. Sadick
#Majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#---------------------------------------------------------------------------------------------------
#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load CR4_LEN_all_so_30PC SeuratObject
so <- readRDS(file.path("file_path", "CR4_LEN_all_so_30PC.rds"))
sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so <- SetIdent(so, value = "integrated_snn_res.0.1")
so@meta.data$cluster_id <- Idents(so)
sce$cluster_id <- Idents(so)
(n_cells <- table(sce$cluster_id, sce$sample_id))

DefaultAssay(so) <- "RNA"
DimPlot(so, reduction = "tsne") + theme(aspect.ratio = 1)

#Remove outlier donors (D5 and D9) from SeuratObject
so_noD5 <- subset(x = so, subset = sample_id == "D5", invert = TRUE)
so_noD5D9 <- subset(x = so_noD5, subset = sample_id == "D9", invert = TRUE)

saveRDS(so_noD5D9, file.path("file_path", "CR4_LEN_final_so_unint.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTERING

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load data
so <- readRDS(file.path("file_path", "CR4_LEN_final_so_unint.rds"))

#Convert SCE to SeuratObject
sce <- as.SingleCellExperiment(so, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so <- lapply(cells_by_sample, function(i)
  SubsetData(so, cells = i))

#Normalize, find variable genes, and scale
so <- lapply(so, NormalizeData, verbose = FALSE)
so <- lapply(so, FindVariableFeatures, nfeatures = 2e3,
             selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so <- lapply(so, ScaleData, verbose = FALSE)

#Find anchors and integrate
as <- FindIntegrationAnchors(so, verbose = TRUE)
so <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so) <- "integrated"
so <- ScaleData(so, display.progress = FALSE)

#DIMENSION REDUCTION 
so <- RunPCA(so, npcs = 100, verbose = FALSE)
ElbowPlot(so, ndims = 100)
#Update number of PCs used
so <- RunTSNE(so, reduction = "pca", dims = seq_len(30),
              seed.use = 1, do.fast = TRUE, verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = seq_len(30),
              seed.use = 1, verbose = FALSE)

#CLUSTERING 
so <- FindNeighbors(so, reduction = "pca", dims = seq_len(30), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so <- FindClusters(so, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#SAVE SeuratObject
saveRDS(so, file.path("file_path", "CR4_LEN_final_so_30PC.rds"))

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
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so <- readRDS(file.path("file_path", "CR4_LEN_final_so_30PC.rds"))
sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so <- SetIdent(so, value = "integrated_snn_res.0.1")
so@meta.data$cluster_id <- Idents(so)
sce$cluster_id <- Idents(so)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/LEN_final_so_PC30_res0.1_cluster_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

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
cs <- sample(colnames(so), 5e3)
.plot_dr <- function(so, dr, id)
  DimPlot(so, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^MT-", x = rownames(so@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so@assays[["RNA"]])
so$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^RP[SL]", x = rownames(so@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so@assays[["RNA"]][rb.genes, ])/Matrix::colSums(so@assays[["RNA"]])
so$percent.rb <- percent.rb

VlnPlot(object = so, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb"), ncol = 4, pt.size = 0)

VlnPlot(object = so, features = c("nFeature_RNA"), pt.size = 0) + stat_summary(fun.y=median, geom="point", shape=23, size=2)

#Generate summary statistics for entire dataset
summary(so$nFeature_RNA)
summary(so$nCount_RNA)

library(pastecs)
stat.desc(so$nFeature_RNA)
stat.desc(so$nCount_RNA)

library(psych)
describe(so$nFeature_RNA)
describe(so$nCount_RNA)

#Generate summary statistics per sample
library(data.table)
library(psych)
feature_by_sample <- as.data.frame(so$nFeature_RNA, row.names = so$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "file_path/LEN_final_so_PC30_res0.1_cluster_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so$nCount_RNA, row.names = so$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/LEN_final_so_PC30_res0.1_cluster_QC_count_by_sample.csv")

#Generate summary statistics per cluster
feature_by_cluster <- as.data.frame(so$nFeature_RNA, row.names = so$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/LEN_final_so_PC30_res0.1_cluster_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so$nCount_RNA, row.names = so$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/LEN_final_so_PC30_res0.1_cluster_QC_count_by_cluster.csv")

#DETERMINE WHAT CELL TYPES ARE PRESENT IN DATASET BEFORE MAKING NEW SEURAT OBJECTS
DefaultAssay(so) <- "RNA"
DimPlot(so, reduction = "tsne") + theme(aspect.ratio = 1)

#Find all markers
DefaultAssay(so) <- "RNA"
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers, "file_path/LEN_final_so_PC30_res0.1_cluster_genes-RNA.csv")

#Visualize
VlnPlot(so, features = c("ALDOC"), pt.size = 0)
DotPlot(so, features = c("AQP4", "GFAP", "FGFR3", "CLDN10", "GJA1", "ALDH1L1", "SLC1A3", "SLC1A2", "HIST1H3E", "TPX2", "NUSAP1", "PPDPF"))
FeaturePlot(so, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#MAKE CELL TYPE-SPECIFIC SEURAT OBJECTS
#Assign cell type identity to clusters
so.renamed <- RenameIdents(so, `0` = "Astro_A", `1` = "Oligo", `2` = "Astro_B", `3` = "Neuro_A", `4` = "Neuro_B", `5`= "Micro", `6` = "OPC", `7`= "Endo-like")
DefaultAssay(so.renamed) <- "RNA"
DimPlot(so.renamed, reduction = "tsne") + theme(aspect.ratio = 1)

#Reorder levels of ident
so_copy <- so.renamed
designated_levels <- c("Astro_A", "Astro_B", "OPC", "Oligo", "Neuro_A", "Neuro_B", "Micro", "Endo-like")
Idents(so_copy) <- factor(Idents(so_copy), levels= designated_levels)
DefaultAssay(so_copy) <- "RNA"
DimPlot(so_copy, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair"))

#Evaluate meta-variables
DimPlot(so_copy, reduction = "tsne", group.by = "sample_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "sex_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "disease_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_copy, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_copy) <- "RNA"
so_copy_cluster.averages <- AverageExpression(so_copy)
head(so_copy_cluster.averages[["RNA"]][, 1:5])
write.csv(so_copy_cluster.averages[["RNA"]], "file_path/LEN_final_so_PC30_res0.1_avg_exp_by_cluster.csv")

orig.levels <- levels(so_copy)
Idents(so_copy) <- gsub(pattern = " ", replacement = "_", x = Idents(so_copy))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_copy) <- orig.levels
so_copy_cluster.averages <- AverageExpression(so_copy, return.seurat = TRUE)
so_copy_cluster.averages

CellScatter(so_copy_cluster.averages, cell1 = "Astro_A", cell2 = "Astro_B")

#Plot 'in silico' bulk datasets
DoHeatmap(so_copy_cluster.averages, features = unlist(TopFeatures(so_copy[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)

DoHeatmap(so_copy_cluster.averages, features = c("AQP4", "GJA1", "ALDOC", "GFAP", "CLU", "ALDH1L1", "SLC1A3", "GLUL", "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "PLP1", "MBP", "MOBP", "MOG", "SNAP25", "STMN2", "SYT1", "SLC17A7", "GRIN2A", "GAD1", "GAD2", "C1QA", "C1QB", "C1QC", "TYROBP", "P2RY12", "HEXB", "TREM2", "CTSS", "CLDN5", "NOSTRIN", "FLT1", "ITM2A", "KLF2", "BSG"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))

#Subset clusters to make cell type-specific objects
so_astro <- subset(x = so.renamed, idents = c("Astro_A", "Astro_B"), invert = FALSE)
so_oligo <- subset(x = so.renamed, idents = c("Oligo"), invert = FALSE)

saveRDS(so_astro, file.path("file_path", "CR4_LEN_final_so_astro_unint.rds"))
saveRDS(so_oligo, file.path("file_path", "CR4_LEN_final_so_oligo_unint.rds"))

#---------------------------------------------------------------------------------------------------
#RECLUSTER ASTROCYTE-SPECIFIC SeuratObject

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
#memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load data
so_astro <- readRDS(file.path("file_path", "CR4_LEN_final_so_astro_unint.rds"))

#Convert SCE to SeuratObject
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_astro <- lapply(cells_by_sample, function(i)
  SubsetData(so_astro, cells = i))

#Normalize, find variable genes, and scale
so_astro <- lapply(so_astro, NormalizeData, verbose = FALSE)
so_astro <- lapply(so_astro, FindVariableFeatures, nfeatures = 2e3,
                   selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_astro <- lapply(so_astro, ScaleData, verbose = FALSE)

#Find anchors and integrate
as <- FindIntegrationAnchors(so_astro, verbose = TRUE)
so_astro <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_astro) <- "integrated"
so_astro <- ScaleData(so_astro, display.progress = FALSE)

#DIMENSION REDUCTION
so_astro <- RunPCA(so_astro, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro, ndims = 50)
#Update number of PCs used
so_astro <- RunTSNE(so_astro, reduction = "pca", dims = seq_len(20),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_astro <- RunUMAP(so_astro, reduction = "pca", dims = seq_len(20),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_astro <- FindNeighbors(so_astro, reduction = "pca", dims = seq_len(20), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_astro <- FindClusters(so_astro, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_astro, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_astro, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#Save SeuratObject
saveRDS(so_astro, file.path("file_path", "CR4_LEN_final_so_astro_20PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION: ASTROCYTES

#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so_astro <- readRDS(file.path("file_path", "CR4_LEN_final_so_astro_20PC.rds"))
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_astro <- SetIdent(so_astro, value = "integrated_snn_res.0.3")
so_astro@meta.data$cluster_id <- Idents(so_astro)
sce$cluster_id <- Idents(so_astro)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/LEN_final_so_astro_20PC_res0.3_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

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
cs <- sample(colnames(so_astro), 5e3)
.plot_dr <- function(so_astro, dr, id)
  DimPlot(so_astro, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_astro, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_astro, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^MT-", x = rownames(so_astro@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_astro@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_astro@assays[["RNA"]])
so_astro$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^RP[SL]", x = rownames(so_astro@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so_astro@assays[["RNA"]][rb.genes, ])/Matrix::colSums(so_astro@assays[["RNA"]])
so_astro$percent.rb <- percent.rb

VlnPlot(object = so_astro, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb"), ncol = 4, pt.size = 0)

#Assess astrocyte clusters
DefaultAssay(so_astro) <- "RNA"
DimPlot(so_astro, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair"))

#Find all markers
DefaultAssay(so_astro) <- "RNA"
so_astro.markers <- FindAllMarkers(so_astro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro.markers, "file_path/LEN_final_so_astro_20PC_res0.3_genes-RNA.csv")

#Visualize
VlnPlot(so_astro, features = c("ALDOC"), pt.size = 0)
DotPlot(so_astro, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(so_astro, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#Evaluate meta-variables
DimPlot(so_astro, reduction = "tsne", group.by = "sample_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "sex_id", cols =  c("seagreen3", "gold1"), pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "disease_id", cols = c("firebrick3", "steelblue3"), pt.size = 0.0001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_astro) <- "RNA"
so_astro_cluster.averages <- AverageExpression(so_astro)
head(so_astro_cluster.averages[["RNA"]][, 1:5])
write.csv(so_astro_cluster.averages[["RNA"]], "file_path/LEN_final_so_astro_20PC_res0.3_avg_exp_by_cluster.csv")

orig.levels <- levels(so_astro)
Idents(so_astro) <- gsub(pattern = " ", replacement = "_", x = Idents(so_astro))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_astro) <- orig.levels
so_astro_cluster.averages <- AverageExpression(so_astro, return.seurat = TRUE)
so_astro_cluster.averages

CellScatter(so_astro_cluster.averages, cell1 = "1", cell2 = "2")

#Plot 'in silico' bulk datasets
DoHeatmap(so_astro_cluster.averages, features = unlist(TopFeatures(so_astro[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE)

DoHeatmap(so_astro_cluster.averages, features = c("HPSE2","CABLES1","ZNF98","ARHGAP24","TRPM3","PLCB1","F11-AS1","CACNB2","SYNE1","SGCD","CST3","APOE","FTH1","HEPN1","ITM2C","FTL","CKB","DLC1","TUBB2B","CLU","AC012405.1","ADAMTSL3","DPP10","LINC00609","L3MBTL4","CACNA2D3","AC073941.1","AQP1","VCAN","MAP1B","CHI3L1","SERPINA3","TPST1","ARHGEF3","CD44","RASGEF1B","SAMD4A","ELL2","GNA14","EMP1","DCLK1","SLC38A1","DPP6","LINC01411","FOS","HSP90AA1","CRYAB","GRIA1","UBC","HSPA1A","ID3","DDIT4","DPP10-AS3","JUN","AQP4","JUNB","LINC01088","KAZN","TTN","TNC","PTCHD4","KCNMB2-AS1","MIR34AHG","AL353138.1","ASCC3","ZMAT3","CDKN1A","DDB2","TNFRSF10B","KCNMB2","MOXD1","DACH1","CREB5","LTBP1","AL589740.1","WIF1","TMEM132C","HGF","ADAMTS17","GRID2","KCNIP4","IL1RAPL1","RBFOX1","CTNNA3","CSMD1","SLC24A2","ST18","MBP","ANK3"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = rev(brewer.pal(8,name = "RdBu")))

#Subset out Astro Cluster 9 due to high MT and ribosomal contamination
so_astro_r2 <- subset(x = so_astro, idents = c("9"), invert = TRUE)
saveRDS(so_astro_r2, file.path("file_path", "CR4_LEN_final_so_astro_r2_unint.rds"))

#---------------------------------------------------------------------------------------------------
#RECLUSTER ASTROCYTE-SPECIFIC SeuratObject ROUND 2

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
#memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load data
so_astro_r2 <- readRDS(file.path("file_path", "CR4_LEN_final_so_astro_r2_unint.rds"))

#Convert SCE to SeuratObject
sce <- as.SingleCellExperiment(so_astro_r2, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_astro_r2 <- lapply(cells_by_sample, function(i)
  SubsetData(so_astro_r2, cells = i))

#Normalize, find variable genes, and scale
so_astro_r2 <- lapply(so_astro_r2, NormalizeData, verbose = FALSE)
so_astro_r2 <- lapply(so_astro_r2, FindVariableFeatures, nfeatures = 2e3,
                      selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_astro_r2 <- lapply(so_astro_r2, ScaleData, verbose = FALSE)

#Find anchors and integrate
as <- FindIntegrationAnchors(so_astro_r2, verbose = TRUE)
so_astro_r2 <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_astro_r2) <- "integrated"
so_astro_r2 <- ScaleData(so_astro_r2, display.progress = FALSE)

#DIMENSION REDUCTION
so_astro_r2 <- RunPCA(so_astro_r2, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro_r2, ndims = 50)
#Update number of PCs used
so_astro_r2 <- RunTSNE(so_astro_r2, reduction = "pca", dims = seq_len(20),
                       seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_astro_r2 <- RunUMAP(so_astro_r2, reduction = "pca", dims = seq_len(20),
                       seed.use = 1, verbose = FALSE)

#CLUSTERING
so_astro_r2 <- FindNeighbors(so_astro_r2, reduction = "pca", dims = seq_len(20), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_astro_r2 <- FindClusters(so_astro_r2, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_astro_r2, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_astro_r2, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#Save SeuratObject
saveRDS(so_astro_r2, file.path("file_path", "CR4_LEN_final_so_astro_r2_20PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION: ASTROCYTES ROUND 2

#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so_astro_r2 <- readRDS(file.path("file_path", "CR4_LEN_final_so_astro_r2_20PC.rds"))
sce <- as.SingleCellExperiment(so_astro_r2, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_astro_r2 <- SetIdent(so_astro_r2, value = "integrated_snn_res.0.3")
so_astro_r2@meta.data$cluster_id <- Idents(so_astro_r2)
sce$cluster_id <- Idents(so_astro_r2)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/LEN_final_so_astro_r2_20PC_res0.3_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

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
cs <- sample(colnames(so_astro_r2), 5e3)
.plot_dr <- function(so_astro_r2, dr, id)
  DimPlot(so_astro_r2, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_astro_r2, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_astro_r2, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^MT-", x = rownames(so_astro_r2@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_astro_r2@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_astro_r2@assays[["RNA"]])
so_astro_r2$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^RP[SL]", x = rownames(so_astro_r2@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so_astro_r2@assays[["RNA"]][rb.genes, ])/Matrix::colSums(so_astro_r2@assays[["RNA"]])
so_astro_r2$percent.rb <- percent.rb

VlnPlot(object = so_astro_r2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb"), ncol = 4, pt.size = 0)

#Assess astrocyte clusters
DefaultAssay(so_astro_r2) <- "RNA"
DimPlot(so_astro_r2, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair"))

#Find all markers
DefaultAssay(so_astro_r2) <- "RNA"
so_astro_r2.markers <- FindAllMarkers(so_astro_r2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro_r2.markers, "file_path/LEN_final_so_astro_r2_20PC_res0.3_genes-RNA.csv")

#Visualize
VlnPlot(so_astro_r2, features = c("ALDOC"), pt.size = 0)
DotPlot(so_astro_r2, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(so_astro_r2, features = c("C3"), reduction = "tsne", pt.size = 0.001, cols = c("lightgrey", "#611ac6"), split.by = "disease_id") + theme(aspect.ratio = 1)

#Evaluate meta-variables
DimPlot(so_astro_r2, reduction = "tsne", group.by = "sample_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro_r2, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro_r2, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro_r2, reduction = "tsne", group.by = "sex_id", cols = c("#035F72", "#EFDC60"), pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro_r2, reduction = "tsne", group.by = "disease_id", cols = c("#d42027", "#094eb2"), pt.size = 0.0001) + theme(aspect.ratio = 1)
DimPlot(so_astro_r2, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 11))
DimPlot(so_astro_r2, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro_r2, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro_r2, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 11))
DimPlot(so_astro_r2, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_astro_r2) <- "RNA"
so_astro_r2_cluster.averages <- AverageExpression(so_astro_r2)
head(so_astro_r2_cluster.averages[["RNA"]][, 1:5])
write.csv(so_astro_r2_cluster.averages[["RNA"]], "file_path/LEN_final_so_astro_r2_20PC_res0.3_avg_exp_by_cluster.csv")

orig.levels <- levels(so_astro_r2)
Idents(so_astro_r2) <- gsub(pattern = " ", replacement = "_", x = Idents(so_astro_r2))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_astro_r2) <- orig.levels
so_astro_r2_cluster.averages <- AverageExpression(so_astro_r2, return.seurat = TRUE)
so_astro_r2_cluster.averages

CellScatter(so_astro_r2_cluster.averages, cell1 = "1", cell2 = "2")

#Plot 'in silico' bulk datasets
DoHeatmap(so_astro_r2_cluster.averages, features = unlist(TopFeatures(so_astro_r2[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE)

DoHeatmap(so_astro_r2_cluster.averages, features = c("HPSE2","CABLES1","ZNF98","ARHGAP24","TRPM3","CST3","APOE","FTH1","HEPN1","ITM2C","DPP10","AC012405.1","L3MBTL4","LINC00609","ADAMTSL3","CHI3L1","SERPINA3","TPST1","SAMD4A","ARHGEF3","DCLK1","SLC38A1","DPP6","LINC01411","FOS","LINC01088","KAZN","HSPA1A","ID3","JUN","DDIT4","PTCHD4","KCNMB2-AS1","MIR34AHG","ZMAT3","AL353138.1","MOXD1","DACH1","SGCD","LTBP1","CREB5"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))

#---------------------------------------------------------------------------------------------------
#RECLUSTER OLIGODENDROCYTE-SPECIFIC SeuratObject

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
#memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load data
so_oligo <- readRDS(file.path("file_path", "CR4_LEN_final_so_oligo_unint.rds"))

#Convert SCE to SeuratObject
sce <- as.SingleCellExperiment(so_oligo, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_oligo <- lapply(cells_by_sample, function(i)
  SubsetData(so_oligo, cells = i))

#Normalize, find variable genes, and scale
so_oligo <- lapply(so_oligo, NormalizeData, verbose = FALSE)
so_oligo <- lapply(so_oligo, FindVariableFeatures, nfeatures = 2e3,
                   selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_oligo <- lapply(so_oligo, ScaleData, verbose = FALSE)

#Find anchors and integrate
##Decrease k.filter to minimize number of oligodendrocytes identified per sample
as <- FindIntegrationAnchors(so_oligo, verbose = TRUE, k.filter = 70)
so_oligo <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_oligo) <- "integrated"
so_oligo <- ScaleData(so_oligo, display.progress = FALSE)

#DIMENSION REDUCTION
so_oligo <- RunPCA(so_oligo, npcs = 50, verbose = FALSE)
ElbowPlot(so_oligo, ndims = 50)
#Update number of PCs used
so_oligo <- RunTSNE(so_oligo, reduction = "pca", dims = seq_len(15),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_oligo <- RunUMAP(so_oligo, reduction = "pca", dims = seq_len(15),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_oligo <- FindNeighbors(so_oligo, reduction = "pca", dims = seq_len(15), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_oligo <- FindClusters(so_oligo, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_oligo, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_oligo, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#Save SeuratObject
saveRDS(so_oligo, file.path("file_path", "CR4_LEN_final_so_oligo_15PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION: OLIGODENDROCYTES

#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so_oligo <- readRDS(file.path("file_path", "CR4_LEN_final_so_oligo_15PC.rds"))
sce <- as.SingleCellExperiment(so_oligo, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_oligo <- SetIdent(so_oligo, value = "integrated_snn_res.0.1")
so_oligo@meta.data$cluster_id <- Idents(so_oligo)
sce$cluster_id <- Idents(so_oligo)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/LEN_final_so_oligo_20PC_res0.1_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

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
cs <- sample(colnames(so_oligo), 5e3)
.plot_dr <- function(so_oligo, dr, id)
  DimPlot(so_oligo, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_oligo, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_oligo, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^MT-", x = rownames(so_oligo@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_oligo@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_oligo@assays[["RNA"]])
so_oligo$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^RP[SL]", x = rownames(so_oligo@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so_oligo@assays[["RNA"]][rb.genes, ])/Matrix::colSums(so_oligo@assays[["RNA"]])
so_oligo$percent.rb <- percent.rb

VlnPlot(object = so_oligo, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb"), ncol = 4, pt.size = 0)

VlnPlot(object = so_oligo, features = c("nFeature_RNA"), pt.size = 0) + stat_summary(fun.y=median, geom="point", shape=23, size=2)

#Assess oligo clusters
DefaultAssay(so_oligo) <- "RNA"
DimPlot(so_oligo, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair"))

#Find all markers
DefaultAssay(so_oligo) <- "RNA"
so_oligo.markers <- FindAllMarkers(so_oligo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_oligo.markers,"file_path/LEN_final_so_oligo_20PC_res0.2_genes-RNA.csv")

#Visualize
VlnPlot(so_oligo, features = c("PLP1"), pt.size = 0)
DotPlot(so_oligo, features = c("PLP1", "MBP", "MOG", "OLIG2"))
FeaturePlot(so_oligo, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#Evaluate meta-variables
DimPlot(so_oligo, reduction = "tsne", group.by = "sample_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "sex_id", cols =  c("seagreen3", "gold1"), pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "disease_id", cols = c("firebrick3", "steelblue3"), pt.size = 0.0001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 11))
DimPlot(so_oligo, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 11))
DimPlot(so_oligo, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1) 

#Calculate average gene expression within a cluster
DefaultAssay(so_oligo) <- "RNA"
so_oligo_cluster.averages <- AverageExpression(so_oligo)
head(so_oligo_cluster.averages[["RNA"]][, 1:5])
write.csv(so_oligo_cluster.averages[["RNA"]], "file_path/LEN_final_so_oligo_20PC_res0.1_avg_exp_by_cluster.csv")

orig.levels <- levels(so_oligo)
Idents(so_oligo) <- gsub(pattern = " ", replacement = "_", x = Idents(so_oligo))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_oligo) <- orig.levels
so_oligo_cluster.averages <- AverageExpression(so_oligo, return.seurat = TRUE)
so_oligo_cluster.averages

CellScatter(so_oligo_cluster.averages, cell1 = "1", cell2 = "2")

#Plot 'in silico' bulk datasets
DoHeatmap(so_oligo_cluster.averages, features = unlist(TopFeatures(so_oligo[["pca"]], balanced = TRUE)), size = 3, draw.lines = FALSE)

DoHeatmap(so_oligo_cluster.averages, features = c("LURAP1L-AS1","SLC5A11","ANKRD18A","HIP1","FP236383.3","FTH1","PLP1","CRYAB","DBNDD2","SELENOP","RBFOX1","RASGRF1","ACSBG1","AFF3","COL18A1","NRG3","ADGRV1","GPM6A","SLC1A2","DPP10","ZC3HAV1","IFIT2","BIRC3","CAMK2D","NAV2"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))
