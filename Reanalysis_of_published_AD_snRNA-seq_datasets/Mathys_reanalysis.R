#Pipeline prepared and Mathys et al. (2019) data reanalyzed by Jessica S. Sadick
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
fastq_dirs <- list.dirs("file_path/GRCh38_premRNA_CR4_Mathys", recursive = FALSE, full.names = TRUE)
names(fastq_dirs) <- basename(fastq_dirs)
sce <- DropletUtils::read10xCounts(fastq_dirs)

#Rename row/colData colnames and SCE dimnames
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))
dimnames(sce) <- list(
  with(rowData(sce), paste(SYMBOL, sep = ".")),
  with(colData(sce), paste(barcode, sample_id, sep = ".")))

#Load metadata 
md_dir <- file.path("file_path/GRCh38_premRNA_CR4_Mathys", "metadata_Mathys.xlsx")
md <- readxl::read_excel(md_dir)
m <- match(sce$sample_id, md$`Sample ID`)
sce$group_id <- md$Characteristics[m]
sce$donor_id <- md$Donor[m]
sce$sex_id <- md$Sex[m]
sce$disease_id <- md$Disease[m]
sce$age_id <- md$Age[m]
sce$APOE_id <- md$APOE[m]
sce$braak_id <- md$braaksc[m]
sce$cerad_id <- md$ceradsc[m]
sce$MMSE_id <- md$MMSE[m]
sce$PMI_id <- md$PMI[m]

#Remove undetected genes
sce <- sce[Matrix::rowSums(counts(sce) > 0) > 0, ]
dim(sce)

#DOUBLET REMOVAL
#Split SCE by sample
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

#SAVE
saveRDS(sce, file.path("file_path", "CR4_Mathys_SCE.rds"))

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
sce <- readRDS(file.path("file_path", "CR4_Mathys_SCE.rds"))

#INTEGRATE
#Create SeuratObject
so <- CreateSeuratObject(
  counts = counts(sce),
  meta.data = data.frame(colData(sce)),
  project = "Mathys_10x_data")

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
##Had to change to reference-based integration because of memory constraints
##Chose 1 donor from each condition with the most nuclei captured
#NS/F refs
D178757 <- which(names(so) == "D17-8757")
#NS/M refs
D178767 <- which(names(so) == "D17-8767")
#AD/F refs
D178762 <- which(names(so) == "D17-8762")
#AD/M refs
D178760 <- which(names(so) == "D17-8760")

as <- FindIntegrationAnchors(so, reference = c(D178757, D178767, D178762, D178760), verbose=TRUE)
so <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so) <- "integrated"
so <- ScaleData(so, display.progress = TRUE)

#DIMENSION REDUCTION
so <- RunPCA(so, npcs = 100, verbose = TRUE)
ElbowPlot(so, ndims = 100)
#Update number of PCs used
so <- RunTSNE(so, reduction = "pca", dims = seq_len(30),
              seed.use = 1, do.fast = TRUE, verbose = TRUE)
so <- RunUMAP(so, reduction = "pca", dims = seq_len(30),
              seed.use = 1, verbose = TRUE)

#CLUSTERING
so <- FindNeighbors(so, reduction = "pca", dims = seq_len(30), verbose = TRUE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1))
  so <- FindClusters(so, resolution = res, random.seed = 1, verbose = TRUE)

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

#Save SeuratObject
saveRDS(so, file.path("file_path", "CR4_Mathys_so-ref_30PC.rds"))

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
so <- readRDS(file.path("file_path", "CR4_Mathys_so-ref_30PC.rds"))
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
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Mathys_so-ref_PC30_res0.1_cluster_numbers.csv")

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
write.csv(feature_by_sample_table, "file_path/Mathys_so-ref_PC30_res0.1_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so$nCount_RNA, row.names = so$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/Mathys_so-ref_PC30_res0.1_QC_count_by_sample.csv")

#Generate summary statistics per cluster
library(data.table)
library(psych)
feature_by_cluster <- as.data.frame(so$nFeature_RNA, row.names = so$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/Mathys_so-ref_PC30_res0.1_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so$nCount_RNA, row.names = so$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/Mathys_so-ref_PC30_res0.1_QC_count_by_cluster.csv")

#DETERMINE WHAT CELL TYPES ARE PRESENT IN DATASET BEFORE MAKING NEW SEURAT OBJECTS
DefaultAssay(so) <- "RNA"
DimPlot(so, reduction = "tsne") + theme(aspect.ratio = 1)

#Find all markers
DefaultAssay(so) <- "RNA"
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers, "file_path/Mathys_so-ref_PC30_res0.1_genes-RNA.csv")

#Visualize
VlnPlot(so, features = c("ALDOC"), pt.size = 0)
DotPlot(so, features = c("AQP4", "GFAP", "FGFR3", "CLDN10", "GJA1", "ALDH1L1", "SLC1A3", "SLC1A2", "HIST1H3E", "TPX2", "NUSAP1", "PPDPF"))
FeaturePlot(so, features = c("ALDOC", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#MAKE CELL TYPE-SPECIFIC SEURAT OBJECTS
#Assign cell type identity to clusters
so.renamed <- RenameIdents(so, `0` = "Oligo", `1` = "Neuro_A", `2` = "Neuro_B", 
                           `3` = "Neuro_C", `4` = "Astro", `5`= "Neuro_D", 
                           `6` = "Neuro_E", `7`= "Neuro_F", `8`= "OPC", `9`= "Micro", 
                           `10`= "Neuro_G", `11`= "Endo-like", `12`= "Neuro_H", 
                           `13`= "Neuro_I", `14`= "Neuro_J", `15`= "Neuro_K", 
                           `16`= "Neuro_L", `17`= "Neuro_M")
DefaultAssay(so.renamed) <- "RNA"
DimPlot(so.renamed, reduction = "tsne") + theme(aspect.ratio = 1)

#Reorder levels of ident
so_copy <- so.renamed
designated_levels <- c("Astro", "OPC", "Oligo", "Neuro_A", "Neuro_B", "Neuro_C", 
                       "Neuro_D", "Neuro_E", "Neuro_F", "Neuro_G", "Neuro_H", 
                       "Neuro_I", "Neuro_J", "Neuro_K", "Neuro_L", "Neuro_M", "Micro", "Endo-like")
Idents(so_copy) <- factor(Idents(so_copy), levels= designated_levels)
DefaultAssay(so_copy) <- "RNA"
DimPlot(so_copy, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair", n = 18))

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
write.csv(so_copy_cluster.averages[["RNA"]], "file_path/Mathys_so-ref_PC30_res0.1_avg_exp_by_cluster.csv")

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

DoHeatmap(so_copy_cluster.averages, features = c("AQP4", "GJA1", "ALDOC", "GFAP", "CLU", "ALDH1L1", "SLC1A3", "GLUL", "PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "PLP1", "MBP", "MOBP", "MOG", "SNAP25", "STMN2", "SYT1", "SLC17A7", "GRIN2A", "GAD1", "GAD2", "C1QA", "C1QB", "C1QC", "TYROBP", "P2RY12", "HEXB", "TREM2", "CTSS", "CLDN5", "NOSTRIN", "FLT1", "ITM2A", "KLF2", "BSG"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair", n = 18)) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))

#Subset clusters to make cell type-specific objects 
#Removed donors with < 50 astrocytes captured (17 donors)
so_1 <- subset(x = so.renamed, subset = sample_id == "D17-8753", invert = TRUE)
so_2 <- subset(x = so_1, subset = sample_id == "D17-8754", invert = TRUE)
so_3 <- subset(x = so_2, subset = sample_id == "D17-8766", invert = TRUE)
so_4 <- subset(x = so_3, subset = sample_id == "D17-8768", invert = TRUE)
so_5 <- subset(x = so_4, subset = sample_id == "D17-8775", invert = TRUE)
so_6 <- subset(x = so_5, subset = sample_id == "D17-8776", invert = TRUE)
so_7 <- subset(x = so_6, subset = sample_id == "D17-8781", invert = TRUE)
so_8 <- subset(x = so_7, subset = sample_id == "D17-8782", invert = TRUE)
so_9 <- subset(x = so_8, subset = sample_id == "D17-8783", invert = TRUE)
so_10 <- subset(x = so_9, subset = sample_id == "D17-8785", invert = TRUE)
so_11 <- subset(x = so_10, subset = sample_id == "D17-8786", invert = TRUE)
so_12 <- subset(x = so_11, subset = sample_id == "D17-8788", invert = TRUE)
so_13 <- subset(x = so_12, subset = sample_id == "D17-8789", invert = TRUE)
so_14 <- subset(x = so_13, subset = sample_id == "D17-8791", invert = TRUE)
so_15 <- subset(x = so_14, subset = sample_id == "D17-8795", invert = TRUE)
so_16 <- subset(x = so_15, subset = sample_id == "D17-8797", invert = TRUE)
so_17 <- subset(x = so_16, subset = sample_id == "D17-8800", invert = TRUE)

so_astro <- subset(x = so_17, idents = c("Astro"), invert = FALSE)
saveRDS(so_astro, file.path("file_path", "CR4_Mathys_so_astro_unint.rds"))

#Removed donors with < 50 oligodendrocytes captured (2 donors)
so_A <- subset(x = so.renamed, subset = sample_id == "D17-8783", invert = TRUE)
so_B <- subset(x = so_A, subset = sample_id == "D17-8800", invert = TRUE)

so_oligo <- subset(x = so_B, idents = c("Oligo"), invert = FALSE)
saveRDS(so_oligo, file.path("file_path", "CR4_Mathys_so_oligo_unint.rds"))

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

#Load SeuratObject
so_astro <- readRDS(file.path("file_path", "CR4_Mathys_so_astro_unint.rds"))

#INTEGRATE
#Generate new SCE because of donor removal
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
##Had to change to reference-based integration because of memory constraints
##Chose 1 donor from each condition with the most nuclei captured
#NS/F refs
D178757 <- which(names(so_astro) == "D17-8757")
#NS/M refs
D178767 <- which(names(so_astro) == "D17-8767")
#AD/F refs
D178770 <- which(names(so_astro) == "D17-8770")
#AD/M refs
D178756 <- which(names(so_astro) == "D17-8756")

##HAD TO DECREASE k.filter PARAMETER (because of donors with low astro capture rates)
as <- FindIntegrationAnchors(so_astro, reference = c(D178757, D178767, D178770, D178756), verbose = TRUE, k.filter = 50)
so_astro <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_astro) <- "integrated"
so_astro <- ScaleData(so_astro, display.progress = TRUE)

#DIMENSION REDUCTION
so_astro <- RunPCA(so_astro, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro, ndims = 50)

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
saveRDS(so_astro, file.path("file_path", "CR4_Mathys_so_astro-ref_15PC.rds"))

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
library(inauguration)
library(wesanderson)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so_astro <- readRDS(file.path("file_path", "CR4_Mathys_so_astro-ref_15PC.rds"))
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
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Mathys_so_astro-ref_PC15_res0.4_numbers.csv")

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
##Changed from 5e3 to 3e3
cs <- sample(colnames(so_astro), 3e3)
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
DimPlot(so_astro, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pair"))

#Find all markers
DefaultAssay(so_astro) <- "RNA"
so_astro.markers <- FindAllMarkers(so_astro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro.markers, "file_path/Mathys_so_astro-ref_PC15_res0.4_genes-RNA.csv")

#Visualize
VlnPlot(so_astro, features = c("C3"), pt.size = 0)
DotPlot(so_astro, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(so_astro, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#Evaluate meta-variables
DimPlot(so_astro, reduction = "tsne", group.by = "sample_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 31))
DimPlot(so_astro, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "sex_id", pt.size = 0.001, cols = c("#035F72", "#EFDC60")) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "disease_id", cols = c("#d42027", "#094eb2"), pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 22))
DimPlot(so_astro, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_astro, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 28))
DimPlot(so_astro, reduction = "tsne", group.by = "APOE_id", pt.size = 0.001, cols = c("#008ece", "#e0607e", "#eca0b2", "#077187", "#8e2043")) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_astro) <- "RNA"
so_astro_cluster.averages <- AverageExpression(so_astro)
head(so_astro_cluster.averages[["RNA"]][, 1:5])
write.csv(so_astro_cluster.averages[["RNA"]], "file_path/Mathys_so_astro-ref_PC15_res0.4_avg_exp_by_cluster.csv")

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

DoHeatmap(so_astro_cluster.averages, features = c("AC012404.1","LSAMP-AS1","FTH1","CST3","APOE","MT-ND3","MT-ND2","MT-CYB","MT-ND1","MT-ATP6","DPP10","AC012405.1","DCLK1","LINC00609","LINC01088","KCNIP4","CSMD1","CNTNAP2","RBFOX1","OPCML","IL1RAPL1","PLP1","ST18","CTNNA3","MBP"), size = 3,
          draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))

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
so_oligo <- readRDS(file.path("file_path", "CR4_Mathys_so_oligo_unint.rds"))

#Create SCE from SeuratObject
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
##Had to change to reference-based integration because of memory constraints
##Chose 1 donor from each condition with the most nuclei captured
#NS/F refs
D178757 <- which(names(so_oligo) == "D17-8757")
#NS/M refs
D178759 <- which(names(so_oligo) == "D17-8759")
#AD/F refs
D178780 <- which(names(so_oligo) == "D17-8780")
#AD/M refs
D178768 <- which(names(so_oligo) == "D17-8768")

##HAD TO DECREASE k.filter PARAMETER (because of donors with low oligo capture rates)
as <- FindIntegrationAnchors(so_oligo, reference = c(D178757, D178759, D178780, D178768), verbose = TRUE, k.filter = 80)
so_oligo <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE)

#Scale integrated data
DefaultAssay(so_oligo) <- "integrated"
so_oligo <- ScaleData(so_oligo, display.progress = TRUE)

#DIMENSION REDUCTION
so_oligo <- RunPCA(so_oligo, npcs = 50, verbose = FALSE)
ElbowPlot(so_oligo, ndims = 50)

so_oligo <- RunTSNE(so_oligo, reduction = "pca", dims = seq_len(12),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_oligo <- RunUMAP(so_oligo, reduction = "pca", dims = seq_len(12),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_oligo <- FindNeighbors(so_oligo, reduction = "pca", dims = seq_len(12), verbose = FALSE)
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
saveRDS(so_oligo, file.path("file_path", "CR4_Mathys_so_oligo-ref_12PC.rds"))

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
so_oligo <- readRDS(file.path("file_path", "CR4_Mathys_so_oligo-ref_12PC.rds"))
sce <- as.SingleCellExperiment(so_oligo, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_oligo <- SetIdent(so_oligo, value = "integrated_snn_res.0.3")
so_oligo@meta.data$cluster_id <- Idents(so_oligo)
sce$cluster_id <- Idents(so_oligo)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Mathys_so_oligo-ref_12PC_res0.3_numbers.csv")

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
write.csv(so_oligo.markers, "file_path/Mathys_so_oligo-ref_12PC_res0.3_genes-RNA.csv")

#Visualize
VlnPlot(so_oligo, features = c("PLP1"), pt.size = 0)
DotPlot(so_oligo, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(so_oligo, features = c("PLP1", "MBP", "MOG", "OLIG2", "GFAP", "SNAP25"), reduction = "tsne")

#Evaluate meta-variables
DimPlot(so_oligo, reduction = "tsne", group.by = "sample_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "sex_id", cols =  c("#035F72", "#EFDC60"), pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "disease_id", cols = c("#d42027", "#094eb2"), pt.size = 0.0001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_oligo, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 35))
DimPlot(so_oligo, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1) +scale_color_manual(values = usecol("pal_unikn_pref", n = 40))
DimPlot(so_oligo, reduction = "tsne", group.by = "APOE_id", pt.size = 0.001, cols = c("#008ece", "#e0607e", "#eca0b2", "#077187", "#8e2043")) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_oligo) <- "RNA"
so_oligo_cluster.averages <- AverageExpression(so_oligo)
head(so_oligo_cluster.averages[["RNA"]][, 1:5])
write.csv(so_oligo_cluster.averages[["RNA"]], "file_path/Mathys_so_oligo-ref_12PC_res0.3_avg_exp_by_cluster.csv")

orig.levels <- levels(so_oligo)
Idents(so_oligo) <- gsub(pattern = " ", replacement = "_", x = Idents(so_oligo))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_oligo) <- orig.levels
so_oligo_cluster.averages <- AverageExpression(so_oligo, return.seurat = TRUE)
so_oligo_cluster.averages

CellScatter(so_oligo_cluster.averages, cell1 = "1", cell2 = "2")

#Plot 'in silico' bulk datasets
DoHeatmap(so_oligo_cluster.averages, features = unlist(TopFeatures(so_oligo[["pca"]], balanced = TRUE)), size = 3, 
          draw.lines = FALSE)

DoHeatmap(so_oligo_cluster.averages, features = c("FTH1","CRYAB","DBNDD2","RNASE1","PPP1R14A","XIST","RBFOX1","RASGRF1","SLC5A11","LINC00609","ANKRD18A","SGCZ","MDGA2","CNTN1","FRY","KCNIP4","AL033523.1","PLA2G4C","AATK","DPYSL5","GRID2"), size = 3,
          draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))
