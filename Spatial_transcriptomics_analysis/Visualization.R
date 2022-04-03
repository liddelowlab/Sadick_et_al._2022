# This script provides the code for the visualization of astrocyte cluster gene signature modules in the Hasel et al. (2021) and Maynard et al. (2021) spatial transcriptomics datasets  

# load required packages
library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(readxl)
library(Seurat)
library(cowplot)
library(ggpubr)

### get package requirements file
pkgs <- loadedNamespaces()
for (i in 1:length(pkgs)) {
  cat(print(paste(pkgs[i], getNamespaceVersion(pkgs[i]), sep = " ")), file = "Visualization_packages.txt", append = TRUE)
  cat("\n", file = "Visualization_packages.txt", append = TRUE)
}

# load representative Maynard Visium section from spaceranger output files
hu.visium <- Load10X_Spatial("file_path/Maynard_2021/outs/")

# normalize count data
hu.visium <- SCTransform(hu.visium, assay = "Spatial", verbose = TRUE, return.only.var.genes = FALSE)

# create SingleCellExperiment object from visium data in bayesspace format
colData <- read.csv("/file_path/Maynard_2021/outs/tissue_positions_list.csv", header=FALSE)
colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
rownames(colData) <- colData$spot
colData <- colData[colData$in_tissue > 0, ] 

rowData <- read.table("/file_path/Maynard_2021/filtered_feature_bc_matrix/features.tsv.gz", header=FALSE)
colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
rowData <- rowData[, c("gene_id", "gene_name")]
rownames(rowData) <- make.unique(rowData$gene_name) 
genes.list = rownames(hu.visium)

rowData <- subset(rowData, rownames(rowData) %in% genes.list)

logcounts <- as.matrix(hu.visium@assays$SCT@data)
barcodes <- colnames(hu.visium)
rownames(logcounts) <- rownames(rowData)

colnames(logcounts) <- barcodes
logcounts <- logcounts[, rownames(colData)]

hu.sce <- SingleCellExperiment(assays=list(logcounts=as(logcounts, "dgCMatrix")),
                               rowData=rowData,
                               colData=colData)

metadata(hu.sce)$BayesSpace.data <- list()
metadata(hu.sce)$BayesSpace.data$platform <- "Visium"
metadata(hu.sce)$BayesSpace.data$is.enhanced <- FALSE

# BayeSpace preprocessing
hu.sce <- spatialPreprocess(hu.sce, platform="Visium", n.PCs=15, n.HVGs=2000, log.normalize=FALSE) 
hu.sce <- qTune(hu.sce, qs=seq(1, 20), platform="Visium", d=15)

# plot q plot
tiff("human_spatial_cluster_qplot.tiff", units = "in", width = 5, height = 5, res = 300)
qPlot(hu.sce) # picking 7 clusters
dev.off()

# spatial clustering
set.seed(167)
hu.sce <- spatialCluster(hu.sce, q=7, platform="Visium", d=15, save.chain=TRUE)

# plot spatial clusters
tiff("human_spatial_clusters.tiff", units = "in", width = 8, height = 8, res = 300)
clusterPlot(hu.sce) + theme(aspect.ratio = 1) + scale_y_reverse() + coord_flip()
dev.off()

# spatial cluster resolution enhancement
hu.sce.enhanced <- spatialEnhance(hu.sce, q=7, platform="Visium", d=15, save.chain=TRUE, verbose = TRUE) 

# plot resolution enhanced spatial cluster
tiff("human_spatial_clustering_enhanced.tiff", units = "in", height = 8, width = 8, res= 300)
clusterPlot(hu.sce.enhanced) + theme(aspect.ratio = 1) + scale_y_reverse() + coord_flip()
dev.off()

# save BayesSpace SCE objects, standard and enhanced resolution 
saveRDS(hu.sce, "hu_visium_bayes_sct.rds")
saveRDS(hu.sce.enhanced, "hu_visium_bayes_enhanced_sct.rds")

# save Visium seurat object
saveRDS(hu.visium, "hu_visium_seurat_object.rds")

# get vector of gene names
all_genes = rownames(hu.sce)

# spatial gene expression resolution enhancement
hu.sce.enhanced.allfeatures <- enhanceFeatures(hu.sce.enhanced, hu.sce, nrounds=0, feature_names = all_genes) 

# save SCE object with resolution enhanced gene expression
saveRDS(hu.sce.enhanced.allfeatures, "hu_visium_bayes_enhanced_sct_allfeatures.rds")

# load representative mouse Visium section (LPS)
mus.visium <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-113/count-B1/outs/")

# normalize
mus.visium <- SCTransform(mus.visium, assay = "Spatial", verbose = TRUE, return.only.var.genes = FALSE)

# save Visium seurat object
saveRDS(mus.visium, "mus_visium_sct.rds")

# creating SingleCellExperiment object from visium data in the bayesspace format
colData <- read.csv("file_path/Hasel_2021/spaceranger-113/count-B1/outs/spatial/tissue_positions_list.csv", header=FALSE)
colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
rownames(colData) <- colData$spot
colData <- colData[colData$in_tissue > 0, ] 

rowData <- read.table("file_path/Hasel_2021/spaceranger-113/count-B1/outs/filtered_feature_bc_matrix/features.tsv.gz", header=FALSE)
colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
rowData <- rowData[, c("gene_id", "gene_name")]
rownames(rowData) <- make.unique(rowData$gene_name) 
genes.list = rownames(mus.visium)

rowData <- subset(rowData, rownames(rowData) %in% genes.list)

logcounts <- as.matrix(mus.visium@assays$SCT@data)
barcodes <- colnames(mus.visium)
rownames(logcounts) <- rownames(rowData)

colnames(logcounts) <- barcodes
logcounts <- logcounts[, rownames(colData)]

mus.sce <- SingleCellExperiment(assays=list(logcounts=as(logcounts, "dgCMatrix")),
                                rowData=rowData,
                                colData=colData)

metadata(mus.sce)$BayesSpace.data <- list()
metadata(mus.sce)$BayesSpace.data$platform <- "Visium"
metadata(mus.sce)$BayesSpace.data$is.enhanced <- FALSE

# save SCE object
saveRDS(mus.sce, "mus_visium_sct_sce.rds")

# load Visium seurat object from spaceranger output file
mus.visium.cnt <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-114/count-C1/outs")

# normalize count data
mus.visium.cnt <- SCTransform(mus.visium.cnt, assay = "Spatial", verbose = TRUE, return.only.var.genes = FALSE)

# save seurat object
saveRDS(mus.visium.cnt, "mus_visium_sct_cnt.rds")

# creating SingleCellExperiment object from visium data in the bayesspace format
colData <- read.csv("file_path/Hasel_2021/spaceranger-114/count-C1/outs/spatial/tissue_positions_list.csv", header=FALSE)
colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
rownames(colData) <- colData$spot
colData <- colData[colData$in_tissue > 0, ] 

rowData <- read.table("file_path/Hasel_2021/spaceranger-114/count-C1/outs/filtered_feature_bc_matrix/features.tsv.gz", header=FALSE)
colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
rowData <- rowData[, c("gene_id", "gene_name")]
rownames(rowData) <- make.unique(rowData$gene_name) 
genes.list = rownames(mus.visium.cnt)
rowData <- subset(rowData, rownames(rowData) %in% genes.list)
logcounts <- as.matrix(mus.visium.cnt@assays$SCT@data)
barcodes <- colnames(mus.visium.cnt)
rownames(logcounts) <- rownames(rowData)
colnames(logcounts) <- barcodes
logcounts <- logcounts[, rownames(colData)]

mus.sce.cnt <- SingleCellExperiment(assays=list(logcounts=as(logcounts, "dgCMatrix")),
                                    rowData=rowData,
                                    colData=colData)

metadata(mus.sce.cnt)$BayesSpace.data <- list()
metadata(mus.sce.cnt)$BayesSpace.data$platform <- "Visium"
metadata(mus.sce.cnt)$BayesSpace.data$is.enhanced <- FALSE

# save SCE object
saveRDS(mus.sce.cnt, "mus_visium_sct_sce_cnt.rds")

# BayesSpace preprocessing
mus.sce <- spatialPreprocess(mus.sce, platform="Visium", n.PCs=20, n.HVGs=2000, log.normalize=FALSE)

# choosing cluster number
mus.sce <- qTune(mus.sce, qs=seq(1, 30), platform="Visium", d=20)

# plot qplot
tiff("mouse_spatial_cluster_elbowplot_lps.tiff", units = "in", height = 5, width = 5, res = 300)
qPlot(mus.sce)
dev.off()

# spatial clustering 
set.seed(149)
mus.sce <- spatialCluster(mus.sce, q=17, platform="Visium", d=20, save.chain=TRUE)

# plot spatial clusters
tiff("mouse_spatial_clusters_lps.tiff", units = "in", width = 5, height = 5, res = 300)
clusterPlot(mus.sce) + theme(aspect.ratio = 1) + scale_y_reverse() + coord_flip()
dev.off()

# save SCE object
saveRDS(mus.sce, "mus_visium_bayes_sct.rds")

# cluster resolution enhancement
mus.sce.enhanced <- spatialEnhance(mus.sce, q=17, platform="Visium", d=20, save.chain=TRUE, verbose = TRUE)

# plot resolution enhanced spatial clusters
tiff("mouse_spatial_clustering_enhanced_lps.tiff", units = "in", height = 5, width = 5, res= 300)
clusterPlot(mus.sce.enhanced) + theme(aspect.ratio = 1) + scale_y_reverse() + coord_flip()
dev.off()

# save cluster resolution enhanced SCE object
saveRDS(mus.sce.enhanced, "mus_visium_bayes_enhanced_sct.rds")

# read in saline SCE object
mus.sce = readRDS("mus_visium_sct_sce_cnt.rds")

# BayesSpace preprocessing
mus.sce <- spatialPreprocess(mus.sce, platform="Visium", n.PCs=20, n.HVGs=2000, log.normalize=FALSE)

# choosing spatial clusters
mus.sce <- qTune(mus.sce, qs=seq(1, 30), platform="Visium", d=20)

# plot qplot
tiff("mouse_spatial_cluster_elbowplot_saline.tiff", units = "in", height = 5, width = 5, res = 300)
qPlot(mus.sce)
dev.off()

# spatial clustering
set.seed(667)
mus.sce <- spatialCluster(mus.sce, q=18, platform="Visium", d=20, save.chain=TRUE)

# plot spatial clusters
tiff("mouse_spatial_clusters_saline.tiff", units = "in", width = 5, height = 5, res = 300)
clusterPlot(mus.sce) + theme(aspect.ratio = 1) + scale_y_reverse() + coord_flip()
dev.off()

# save saline SCE object
saveRDS(mus.sce, "mus_visium_bayes_sct_cnt.rds")

# spatial cluster resolution enhancement
mus.sce.enhanced <- spatialEnhance(mus.sce, q=18, platform="Visium", d=20, save.chain=TRUE, verbose = TRUE) 

# plot resolution enhanced spatial clusters
tiff("mouse_spatial_clustering_enhanced_saline.tiff", units = "in", height = 5, width = 5, res= 300)
clusterPlot(mus.sce.enhanced) + theme(aspect.ratio = 1) + scale_y_reverse() + coord_flip()
dev.off()

# save SCE object
saveRDS(mus.sce.enhanced, "mus_visium_bayes_enhanced_sct_cnt.rds")

# extract gene names
all_genes = rownames(mus.sce)

# spatial gene expression resolution enhancement 
mus.sce.enhanced.allfeatures <- enhanceFeatures(mus.sce.enhanced, mus.sce, nrounds=0, feature_names = all_genes) 

# save resolution enhanced SCE object
saveRDS(mus.sce.enhanced.allfeatures, "mouse_visium_bayes_enhanced_sct_cnt_allfeatures.rds")

# load cluster gene modules
mouse.modules <- readRDS("file_path/final_mouse_gene_modules_forhaselsections.rds")
human.modules <- readRDS("file_path/final_human_gene_modules_formaynardsections.rds")


# create enhanced seurat object for module score calculation
hu.enhanced.seurat <- as.Seurat(hu.sce.enhanced.allfeatures, data = "logcounts", counts = "logcounts") 

# calculate gene module scores
hu.enhanced.seurat <- AddModuleScore(hu.enhanced.seurat, features = human.modules, name = "Sadick_human_gene_module", assay = "originalexp")

# create spatial enrichment plot legend
leg <- as_ggplot(get_legend(featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module1) + 
                              scale_fill_viridis_c(option = "inferno", guide = guide_colorbar(title.position = "top", title.hjust = 0.5, ticks = FALSE, label = FALSE)) + 
                              theme(aspect.ratio = 1, legend.text=element_text(size=8), legend.position = "top") + labs(fill = "Relative Gene\nModule Enrichment")))

# plot legend
postscript("relative_gene_module_enrichment_legend.ps", width = 6, height = 4)
leg
dev.off()

# plot enrichment scores for each cluster module
postscript("final_bayesspace_figure_panel_cluster0.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module1) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 0") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster1.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module2) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 1") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster2.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module3) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 2") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster3.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module4) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 3") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster4.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module5) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 4") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster5.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module6) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 5") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster6.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module7) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 6") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster7.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module8) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 7") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster8.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = hu.enhanced.seurat$Sadick_human_gene_module9) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 8") + NoLegend()
dev.off()

# save enhanced resolution seurat object with cluster enrichment scores
saveRDS(hu.enhanced.seurat, "maynard_representative_slice_bayesspace_enhanced_seurat_object_with_cluster_module_scores.rds")

# load resolution enhanced SCE objects
cnt.sce.enhanced.allfeatures <- readRDS("file_path/mouse_visium_bayes_enhanced_sct_cnt_allfeatures.rds")
lps.sce.enhanced.allfeatures <- readRDS("file_path/mouse_visium_bayes_enhanced_sct_allfeatures.rds")

# create seurat objects from SCE objects
cnt.sce.enhanced.allfeatures.seurat  <- as.Seurat(cnt.sce.enhanced.allfeatures , data = "logcounts", counts = "logcounts") 
cnt.sce.enhanced.allfeatures.seurat <- AddModuleScore(cnt.sce.enhanced.allfeatures.seurat, features = mouse.modules, name = "Sadick_mouse_gene_module", assay = "originalexp")

lps.sce.enhanced.allfeatures.seurat  <- as.Seurat(lps.sce.enhanced.allfeatures , data = "logcounts", counts = "logcounts") 
lps.sce.enhanced.allfeatures.seurat <- AddModuleScore(lps.sce.enhanced.allfeatures.seurat, features = mouse.modules, name = "Sadick_mouse_gene_module", assay = "originalexp")

# create enrichment plot legend
leg <- as_ggplot(get_legend(featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module1) + 
                              scale_fill_viridis_c(option = "inferno", guide = guide_colorbar(title.position = "top", title.hjust = 0.5, ticks = FALSE, label = FALSE)) + 
                              theme(aspect.ratio = 1, legend.text=element_text(size=8), legend.position = "top") + labs(fill = "Relative Enrichment")))

# plot legend 
postscript("relative_gene_module_enrichment_legend_mouse.ps", width = 6, height = 4)
leg
dev.off()

# cluster gene module enrichment score plots
postscript("final_bayesspace_figure_panel_cluster0_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module1) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 0") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster1_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module2) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 1") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster2_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module3) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 2") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster3_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module4) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 3") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster4_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module5) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 4") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster5_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module6) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 5") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster6_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module7) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 6") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster7_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module8) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 7") + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster8_mouse.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cnt.sce.enhanced.allfeatures.seurat$Sadick_mouse_gene_module9) + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "Cluster 8") + NoLegend()
dev.off()

# creating a merged seurat object for comparative module scoring between LPS and saline sections
cnt.sce.enhanced.allfeatures.seurat  <- as.Seurat(cnt.sce.enhanced.allfeatures , data = "logcounts", counts = "logcounts") 
lps.sce.enhanced.allfeatures.seurat  <- as.Seurat(lps.sce.enhanced.allfeatures , data = "logcounts", counts = "logcounts") 

# add sample metadata
cnt.sce.enhanced.allfeatures.seurat$sample <- "CNT_A"
lps.sce.enhanced.allfeatures.seurat$sample <- "LPS_A"

# merge seurat objects
cntA.plus.lpsA <- merge(cnt.sce.enhanced.allfeatures.seurat, lps.sce.enhanced.allfeatures.seurat)

# calculate cluster enrichment scores
cntA.plus.lpsA <- AddModuleScore(cntA.plus.lpsA, features = mouse.modules, name = "Sadick_mouse_gene_module", assay = "originalexp")

# create data frame of cluster enrichment scores
cntA.plus.lpsA.scores.df <- cntA.plus.lpsA[[c("Sadick_mouse_gene_module1", "Sadick_mouse_gene_module2", "Sadick_mouse_gene_module3",
                                              "Sadick_mouse_gene_module4", "Sadick_mouse_gene_module5", "Sadick_mouse_gene_module6",
                                              "Sadick_mouse_gene_module7", "Sadick_mouse_gene_module8", "Sadick_mouse_gene_module9",
                                              "Sadick_mouse_gene_module10", "Sadick_mouse_gene_module11", "Sadick_mouse_gene_module12",
                                              "Sadick_mouse_gene_module13", "Sadick_mouse_gene_module14", "Sadick_mouse_gene_module15",
                                              "Sadick_mouse_gene_module16", "Sadick_mouse_gene_module17", "Sadick_mouse_gene_module18", "sample")]]
cntA.plus.lpsA.scores.df <- cntA.plus.lpsA.scores.df %>% mutate(barcode = gsub("_[1-2]$", "", rownames(.)))

# plot cluster 3 enrichment score: saline and lps
postscript("final_bayesspace_figure_panel_cluster3_mouse_split_cnt_vs_lps__CNT.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cntA.plus.lpsA.scores.df %>% filter(sample == "CNT_A") %>% .$Sadick_mouse_gene_module4) + 
  scale_fill_viridis_c(option = "inferno", lim = range(c(cntA.plus.lpsA.scores.df %>% filter(sample == "CNT_A") %>% .$Sadick_mouse_gene_module4, cntA.plus.lpsA.scores.df %>% filter(sample == "LPS_A") %>% .$Sadick_mouse_gene_module4))) + 
  theme(aspect.ratio = 1, legend.text=element_text(size=8)) + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster3_mouse_split_cnt_vs_lps__LPS.ps", width = 4, height = 4)
featurePlot(lps.sce.enhanced.allfeatures, feature = cntA.plus.lpsA.scores.df %>% filter(sample == "LPS_A") %>% .$Sadick_mouse_gene_module4) + 
  scale_fill_viridis_c(option = "inferno", lim = range(c(cntA.plus.lpsA.scores.df %>% filter(sample == "CNT_A") %>% .$Sadick_mouse_gene_module4, cntA.plus.lpsA.scores.df %>% filter(sample == "LPS_A") %>% .$Sadick_mouse_gene_module4))) + 
  theme(aspect.ratio = 1, legend.text=element_text(size=8)) + NoLegend()
dev.off()

# plot cluster 5 AD gene module enrichment: saline vs lps
postscript("final_bayesspace_figure_panel_cluster5_ad_genes_mouse_split_cnt_vs_lps__CNT.ps", width = 4, height = 4)
featurePlot(cnt.sce.enhanced.allfeatures, feature = cntA.plus.lpsA.scores.df %>% filter(sample == "CNT_A") %>% .$Sadick_mouse_gene_module15) + 
  scale_fill_viridis_c(option = "inferno", lim = range(c(cntA.plus.lpsA.scores.df %>% filter(sample == "CNT_A") %>% .$Sadick_mouse_gene_module15, cntA.plus.lpsA.scores.df %>% filter(sample == "LPS_A") %>% .$Sadick_mouse_gene_module15))) + 
  theme(aspect.ratio = 1, legend.text=element_text(size=8)) + NoLegend()
dev.off()

postscript("final_bayesspace_figure_panel_cluster5_ad_genes_mouse_split_cnt_vs_lps__LPS.ps", width = 4, height = 4)
featurePlot(lps.sce.enhanced.allfeatures, feature = cntA.plus.lpsA.scores.df %>% filter(sample == "LPS_A") %>% .$Sadick_mouse_gene_module15) + 
  scale_fill_viridis_c(option = "inferno", lim = range(c(cntA.plus.lpsA.scores.df %>% filter(sample == "CNT_A") %>% .$Sadick_mouse_gene_module15, cntA.plus.lpsA.scores.df %>% filter(sample == "LPS_A") %>% .$Sadick_mouse_gene_module15))) + 
  theme(aspect.ratio = 1, legend.text=element_text(size=8)) + NoLegend()
dev.off()

# save final enrichment scores
cnt.A.module.scores.df <- cntA.plus.lpsA[[c("Sadick_mouse_gene_module1", "Sadick_mouse_gene_module2", "Sadick_mouse_gene_module3",
                                            "Sadick_mouse_gene_module4", "Sadick_mouse_gene_module5", "Sadick_mouse_gene_module6",
                                            "Sadick_mouse_gene_module7", "Sadick_mouse_gene_module8", "Sadick_mouse_gene_module9",
                                            "Sadick_mouse_gene_module10", "Sadick_mouse_gene_module11", "Sadick_mouse_gene_module12",
                                            "Sadick_mouse_gene_module13", "Sadick_mouse_gene_module14", "Sadick_mouse_gene_module15",
                                            "Sadick_mouse_gene_module16", "Sadick_mouse_gene_module17", "Sadick_mouse_gene_module18", "sample")]]

saveRDS(cnt.A.module.scores.df, "cntA_bayesspace_enhanced_modulescores_dataframe_for_final_figures_enrichment_plots.rds")
write.csv(cnt.A.module.scores.df, "cntA_bayesspace_enhanced_modulescores_dataframe_for_final_figures_enrichment_plots.csv")

saveRDS(cntA.plus.lpsA.scores.df, "cntA_plus_lpsA_merged_bayesspace_enhanced_module_scores_dataframe_for_final_figure.rds")
write.csv(cntA.plus.lpsA.scores.df, "cntA_plus_lpsA_merged_bayesspace_enhanced_module_scores_dataframe_for_final_figure.csv")

# save seurat objects with module scores
saveRDS(cnt.sce.enhanced.allfeatures.seurat, "cnt_A_final_figures_bayesspace_enhanced_seurat_object_with_cluster_module_scores.rds")
saveRDS(lps.sce.enhanced.allfeatures.seurat, "lps_A_final_figures_bayesspace_enhanced_seurat_object_with_cluster_module_scores.rds")
saveRDS(cntA.plus.lpsA, "cntA_plus_lpsA_merged_final_figures_bayesspace_enhanced_seurat_object_with_cluster_module_scores_for_direct_comparison.rds")

# cluster 6 unenhanced enrichment score plot
leg = as_ggplot(get_legend(SpatialFeaturePlot(maynard.combo, features = "Sadick_human_gene_module7", images = "slice1.8") + 
                             scale_fill_viridis_c(option = "inferno", name = "Gene Module Score", guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
                             theme(aspect.ratio = 1, legend.position = "top")))

postscript("cluster6_unenhanced.ps", width = 4, height = 4)
SpatialFeaturePlot(maynard.combo, features = "Sadick_human_gene_module7", images = "slice1.8") + 
  scale_fill_viridis_c(option = "inferno", name = "Gene Module Score", guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
  theme(aspect.ratio = 1, legend.position = "top") + NoLegend()
dev.off()

postscript("cluster6_unenhanced_legend.ps", width = 4, height = 4)
leg
dev.off()

# plot resolution enhanced Aqp4 expression on human Visium section
postscript("eg_enhanced_human_plot_aqp4.ps", width = 4, height = 4)
featurePlot(hu.sce.enhanced.allfeatures, feature = "AQP4") + scale_fill_viridis_c(option = "inferno") + theme(aspect.ratio = 1, legend.text=element_text(size=8)) + labs(fill = "AQP4") + NoLegend()
dev.off()
