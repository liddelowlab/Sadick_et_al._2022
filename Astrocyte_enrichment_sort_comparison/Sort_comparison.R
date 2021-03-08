#Astrocyte enrichment sort strategy comparison (SOX9+ and SOX9- versus LHX2+/NeuN-)
#Goal: To evaluate two potential astrocytic nuclei enrichment methods using a single donor
#Pipeline prepared by Jessica S. Sadick
#Majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#---------------------------------------------------------------------------------------------------
#Load libraries
library(cowplot)
library(ggplot2)
library(scater)
library(scds)
library(SingleCellExperiment)

#LOAD AND REFORMAT DATA
#Load raw counts
fastq_dirs <- list.dirs("file_path/GRCh38_premRNA_CR4_sort_comparison", recursive = FALSE, full.names = TRUE)
names(fastq_dirs) <- basename(fastq_dirs)
sce <- DropletUtils::read10xCounts(fastq_dirs)

#Rename row/colData colnames & SCE dimnames
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))
dimnames(sce) <- list(
  with(rowData(sce), paste(SYMBOL, sep = ".")),
  with(colData(sce), paste(barcode, sample_id, sep = ".")))

#Load metadata 
md_dir <- file.path("file_path/GRCh38_premRNA_CR4_sort_comparison", "metadata_sort_comparison.xlsx")
md <- readxl::read_excel(md_dir)
m <- match(sce$sample_id, md$`Sample ID`)
sce$group_id <- md$Characteristics[m]
sce$donor_id <- md$Donor[m]
sce$sex_id <- md$Sex[m]
sce$disease_id <- md$Disease[m]
sce$age_id <- md$Age[m]
sce$enrichment_id <- md$Enrichment[m]
sce$RIN_id <- md$RIN[m]
sce$APOE_id <- md$APOE[m]

#Remove undetected genes
sce <- sce[Matrix::rowSums(counts(sce) > 0) > 0, ]
dim(sce)

#DOUBLET REMOVAL
#Split by sample
cs_by_s <- split(colnames(sce), sce$sample_id)
sce_by_s <- lapply(cs_by_s, function(cs) sce[, cs])

#Run 'scds'
sce_by_s <- lapply(sce_by_s, function(u) 
  cxds_bcds_hybrid(bcds(cxds(u))))

#Remove doublets
sce_by_s <- lapply(sce_by_s, function(u) {
  #Compute expected nb. of doublets (10x)
  n_dbl <- ceiling(0.01 * ncol(u)^2 / 1e3)
  #Remove 'n_dbl' cells w/ highest doublet score
  o <- order(u$hybrid_score, decreasing = TRUE)
  u[, -o[seq_len(n_dbl)]]
})

#Merge back into single SCE
sce <- do.call(cbind, sce_by_s)

#CALCULATE QC METRICS
(mito <- grep("MT-", rownames(sce), value = TRUE))
sce <- calculateQCMetrics(sce, feature_controls = list(Mt = mito))
#plotHighestExprs(sce, n = 20)

#FILTERING
#Get sample-specific outliers
cols <- c("total_counts", "total_features_by_counts", "pct_counts_Mt")
log <- c(TRUE, TRUE, FALSE)
type <- c("both", "both", "higher")

drop_cols <- paste0(cols, "_drop")
for (i in seq_along(cols))
  colData(sce)[[drop_cols[i]]] <- isOutlier(sce[[cols[i]]], 
                                            nmads = 2.5, type = type[i], log = log[i], batch = sce$sample_id)

sapply(drop_cols, function(i) 
  sapply(drop_cols, function(j)
    sum(sce[[i]] & sce[[j]])))

cd <- data.frame(colData(sce))
ps <- lapply(seq_along(cols), function (i) {
  p <- ggplot(cd, aes_string(x = cols[i], alpha = drop_cols[i])) +
    geom_histogram(bins = 100, show.legend = FALSE) +
    scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.4)) +
    facet_wrap(~sample_id, ncol = 1, scales = "free") + 
    theme_classic() + theme(strip.background = element_blank())
  if (log[i]) 
    p <- p + scale_x_log10()
  return(p)
})
plot_grid(plotlist = ps, ncol = 3)

layout(matrix(1:2, nrow = 1))
ol <- Matrix::rowSums(as.matrix(colData(sce)[drop_cols])) != 0
x <- sce$total_counts
y <- sce$total_features_by_counts
LSD::heatscatter(x, y, log="xy", main = "unfiltered", 
                 xlab = "Total counts", ylab = "Non-zero features")
LSD::heatscatter(x[!ol], y[!ol], log="xy", main = "filtered", 
                 xlab = "Total counts", ylab = "Non-zero features")

#Generate summary of cells kept
ns <- table(sce$sample_id)
ns_fil <- table(sce$sample_id[!ol])
print(rbind(
  unfiltered = ns, filtered = ns_fil, 
  "%" = ns_fil / ns * 100), digits = 0)

#Drop outlier cells
sce <- sce[, !ol]
dim(sce)

#Require count > 1 in at least 20 cells
sce <- sce[Matrix::rowSums(counts(sce) > 1) >= 20, ]
dim(sce)

#SAVE SCE
saveRDS(sce, file.path("file_path", "CR4_sort_comparison_SCE.rds"))

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
sce <- readRDS(file.path("file_path", "CR4_sort_comparison_SCE.rds"))

#INTEGRATE
#Create SeuratObject
so <- CreateSeuratObject(
  counts = counts(sce),
  meta.data = data.frame(colData(sce)),
  project = "sort_comparison_data")

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
saveRDS(so, file.path("file_path", "CR4_sort_comparison_so_30PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION

#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(viridis)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so <- readRDS(file.path("file_path", "CR4_sort_comparison_so_30PC.rds"))
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
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/sort_comparison_so_PC30_res0.1_cluster_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue"), gids)

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
write.csv(feature_by_sample_table, "file_path/sort_comparison_PC30_res0.1_cluster_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so$nCount_RNA, row.names = so$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/sort_comparison_PC30_res0.1_cluster_QC_count_by_sample.csv")

#Generate summary statistics per cluster
feature_by_cluster <- as.data.frame(so$nFeature_RNA, row.names = so$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/sort_comparison_PC30_res0.1_cluster_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so$nCount_RNA, row.names = so$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/sort_comparison_PC30_res0.1_cluster_QC_count_by_cluster.csv")

#DETERMINE WHAT CELL TYPES ARE PRESENT IN DATASET BEFORE MAKING NEW SEURAT OBJECTS
DefaultAssay(so) <- "RNA"
DimPlot(so, reduction = "tsne") + theme(aspect.ratio = 1)

#Find all markers
DefaultAssay(so) <- "RNA"
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers, "file_path/sort_comparison_so_PC30_res0.1_genes-RNA.csv")

#Visualize
VlnPlot(so, features = c("ALDOC"), pt.size = 0)
DotPlot(so, features = c("AQP4", "GFAP", "FGFR3", "CLDN10", "GJA1", "ALDH1L1", "SLC1A3", "SLC1A2", "HIST1H3E", "TPX2", "NUSAP1", "PPDPF"))
FeaturePlot(so, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#MAKE CELL TYPE-SPECIFIC SEURAT OBJECTS
#Assign cell type identity to clusters
so.renamed <- RenameIdents(so, `0` = "Astro", `1` = "Oligo", `2` = "OPC", `3` = "Neuro_A", `4` = "Micro", `5`= "Neuro_B", `6` = "Neuro_C", `7`= "Neuro_D", `8`= "Endo-like_A", `9` = "Neuro_E", `10` = "Neuro_F", `11` = "Neuro_G", `12` = "Endo-like_B", `13` = "Neuro_H", `14` = "Endo")
DefaultAssay(so.renamed) <- "RNA"
DimPlot(so.renamed, reduction = "tsne") + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair", 15))

#Reorder levels of ident
so_copy <- so.renamed
designated_levels <- c("Astro", "OPC", "Oligo", "Neuro_A", "Neuro_B", "Neuro_C", "Neuro_D", "Neuro_E", "Neuro_F", "Neuro_G", "Neuro_H", "Micro", "Endo", "Endo-like_A", "Endo-like_B")
Idents(so_copy) <- factor(Idents(so_copy), levels= designated_levels)
DefaultAssay(so_copy) <- "RNA"
DimPlot(so_copy, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair",22))

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
DimPlot(so_copy, reduction = "tsne", group.by = "enrichment_id", pt.size = 0.001, cols = c("#df3c22", "#5698af", "#e6cc4f")) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_copy) <- "RNA"
so_copy_cluster.averages <- AverageExpression(so_copy)
head(so_copy_cluster.averages[["RNA"]][, 1:5])
write.csv(so_copy_cluster.averages[["RNA"]], "file_path/sort_comparison_so_PC30_res0.1_avg_exp_by_cluster.csv")

orig.levels <- levels(so_copy)
Idents(so_copy) <- gsub(pattern = " ", replacement = "_", x = Idents(so_copy))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_copy) <- orig.levels
so_copy_cluster.averages <- AverageExpression(so_copy, return.seurat = TRUE)
so_copy_cluster.averages

CellScatter(so_copy_cluster.averages, cell1 = "Astro", cell2 = "Oligo")

#Plot 'in silico' bulk datasets
DoHeatmap(so_copy_cluster.averages, features = unlist(TopFeatures(so_copy[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)

DoHeatmap(so_copy_cluster.averages, features = c("AQP4", "GFAP", "ALDOC", "GFAP", "CLU", "ALDH1L1", "SLC1A3", "GLUL", "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "PLP1", "MBP", "MOBP", "MOG", "SNAP25", "STMN2", "SYT1", "SLC17A7", "GRIN2A", "GAD1", "GAD2", "C1QA", "C1QB", "C1QC", "TYROBP", "P2RY12", "HEXB", "TREM2", "CTSS", "CLDN5", "NOSTRIN", "FLT1", "ITM2A", "KLF2", "BSG"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))

#Subset clusters to make cell type-specific objects 
so_astro <- subset(x = so.renamed, idents = c("Astro"), invert = FALSE)

#---------------------------------------------------------------------------------------------------
#RECLUSTER ASTROCYTE-SPECIFIC SeuratObject

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

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
##Decrease k.filter to minimize number of astrocytes identified per sample
as <- FindIntegrationAnchors(so_astro, verbose = TRUE, k.filter = 125)
so_astro <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_astro) <- "integrated"
so_astro <- ScaleData(so_astro, display.progress = FALSE)

#DIMENSION REDUCTION
so_astro <- RunPCA(so_astro, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro, ndims = 50)
#Update number of PCs used
so_astro <- RunTSNE(so_astro, reduction = "pca", dims = seq_len(15),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_astro <- RunUMAP(so_astro, reduction = "pca", dims = seq_len(15),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_astro <- FindNeighbors(so_astro, reduction = "pca", dims = seq_len(15), verbose = FALSE)
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
saveRDS(so_astro, file.path("file_path/", "CR4_sort_comparison_so_astro_15PC.rds"))

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
so_astro <- readRDS(file.path("file_path", "CR4_sort_comparison_so_astro_15PC.rds"))
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_astro <- SetIdent(so_astro, value = "integrated_snn_res.0.4")
so_astro@meta.data$cluster_id <- Idents(so_astro)
sce$cluster_id <- Idents(so_astro)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/sort_comparison_so_astro_15PC_res0.4_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue"), gids)

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
cs <- sample(colnames(so_astro), 1e3)
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
DimPlot(so_astro, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1)  +scale_color_manual(values = usecol("pal_unikn_pair"))

#Find all markers
DefaultAssay(so_astro) <- "RNA"
so_astro.markers <- FindAllMarkers(so_astro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro.markers, "file_path/sort_comparison_so_astro_15PC_res0.4_genes-RNA.csv")

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
DimPlot(so_astro, reduction = "tsne", group.by = "enrichment_id", pt.size = 0.001, cols = c("#df3c22", "#5698af", "#e6cc4f")) + theme(aspect.ratio = 1)

#Visualize
VlnPlot(so_astro, features = c("ALDOC"), pt.size = 0)
DotPlot(so_astro, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(so_astro, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#Calculate average gene expression within a cluster
DefaultAssay(so_astro) <- "RNA"
so_astro_cluster.averages <- AverageExpression(so_astro)
head(so_astro_cluster.averages[["RNA"]][, 1:5])
write.csv(so_astro_cluster.averages[["RNA"]], "file_path/sort_comparison_so_astro_15PC_res0.4_avg_exp_by_cluster.csv")

orig.levels <- levels(so_astro)
Idents(so_astro) <- gsub(pattern = " ", replacement = "_", x = Idents(so_astro))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_astro) <- orig.levels
so_astro_cluster.averages <- AverageExpression(so_astro, return.seurat = TRUE)
so_astro_cluster.averages

CellScatter(so_astro_cluster.averages, cell1 = "1", cell2 = "2")

#Plot 'in silico' bulk datasets
DoHeatmap(so_astro_cluster.averages, features = unlist(TopFeatures(so_astro[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)

DoHeatmap(so_astro_cluster.averages, features = c("RERG","AL137139.2","LINC01727","CABLES1","SGCD","APOE","FTH1","CST3","HEPN1","ITM2C","AC012405.1","DPP10","AC073941.1","L3MBTL4","AQP1","CHI3L1","GNA14","SERPINA3","ZFP36","SOCS3","DCLK1","KAZN","CD44","SLC38A1","HSP90AA1","LINC01411","CNTN5","IL1RAPL1","ANK3","NRXN3","CTNNA3"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))
