# This script provides the code necessary for the preprocessing of the Hasel et al. (2021) mouse spatial transcriptomics data

# load required packages
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

### get package requirements file
pkgs <- loadedNamespaces()
for (i in 1:length(pkgs)) {
  cat(print(paste(pkgs[i], getNamespaceVersion(pkgs[i]), sep = " ")), file = "Hasel_visium_preprocessing_packages.txt", append = TRUE)
  cat("\n", file = "Hasel_visium_preprocessing_packages.txt", append = TRUE)
}

## loading Visium data from Cell Ranger outputs as Seurat objects; sequencing data is publicly available at NCBI GEO Series GSE165098
lps.1 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-113/count-B1/outs/")
lps.2 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-113/count-A1/outs/")
lps.3 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-114/count-A1/outs/")
lps.4 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-114/count-B1/outs/")

cnt.1 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-114/count-C1/outs/")
cnt.2 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-113/count-C1/outs/")
cnt.3 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-113/count-D1/outs/")
cnt.4 <- Load10X_Spatial("file_path/Hasel_2021/spaceranger-114/count-D1/outs/")

# add sample metadata
lps.1$sample_id <- "LPS_1"
lps.2$sample_id <- "LPS_2"
lps.3$sample_id <- "LPS_3"
lps.4$sample_id <- "LPS_4"

cnt.1$sample_id <- "CNT_1"
cnt.2$sample_id <- "CNT_2"
cnt.3$sample_id <- "CNT_3"
cnt.4$sample_id <- "CNT_4"

# merge Seurat objects
merged.spatial <- merge(lps.1, c(lps.2, lps.3, lps.4, cnt.1, cnt.2, cnt.3, cnt.4))

# plot QC metrics
VlnPlot(merged.spatial, features = "nCount_Spatial", group.by = "sample_id", cols = viridis::viridis_pal()(8), pt.size = 0) + theme(aspect.ratio = 1) + labs(x = NULL, title = NULL, y = "Number of UMIs per spot") + NoLegend()

lim <- range(lps.1@meta.data$nCount_Spatial, lps.2@meta.data$nCount_Spatial, lps.3@meta.data$nCount_Spatial, lps.4@meta.data$nCount_Spatial)
leg <- get_legend(SpatialFeaturePlot(lps.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim, name = "UMIs") + theme(legend.position = "right", aspect.ratio = 1))
plot1 <- SpatialFeaturePlot(lps.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot2 <- SpatialFeaturePlot(lps.2, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot3 <- SpatialFeaturePlot(lps.3, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot4 <- SpatialFeaturePlot(lps.4, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot_grid(plot1, plot2, plot3, plot4, leg, nrow = 1, rel_widths = c(1, 1, 1, 1, 0.2))

lim <- range(cnt.1@meta.data$nCount_Spatial, cnt.2@meta.data$nCount_Spatial, cnt.3@meta.data$nCount_Spatial, cnt.4@meta.data$nCount_Spatial)
leg <- get_legend(SpatialFeaturePlot(cnt.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim, name = "UMIs") + theme(legend.position = "right", aspect.ratio = 1))
plot1 <- SpatialFeaturePlot(cnt.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot2 <- SpatialFeaturePlot(cnt.2, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot3 <- SpatialFeaturePlot(cnt.3, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot4 <- SpatialFeaturePlot(cnt.4, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot_grid(plot1, plot2, plot3, plot4, leg, nrow = 1, rel_widths = c(1, 1, 1, 1, 0.2))

# excluding one CNT and one LPS sample which have low gene/UMI counts
# renaming remaining samples
lps.1$sample_id <- "LPS_A"
lps.2$sample_id <- "LPS_B"
lps.4$sample_id <- "LPS_C"

cnt.1$sample_id <- "CNT_A"
cnt.3$sample_id <- "CNT_B"
cnt.4$sample_id <- "CNT_C"

# merge seurat objects
merged.spatial <- merge(lps.1, c(lps.2, lps.4, cnt.1, cnt.3, cnt.4))

# plot QC metrics
VlnPlot(merged.spatial, features = "nCount_Spatial", group.by = "sample_id", 
        cols = c(viridis::viridis_pal()(8)[1], viridis::viridis_pal()(8)[3], viridis::viridis_pal()(8)[4], viridis::viridis_pal()(8)[5], viridis::viridis_pal()(8)[6], viridis::viridis_pal()(8)[8]), 
        pt.size = 0) + theme(aspect.ratio = 1) + labs(x = NULL, title = NULL, y = "Number of UMIs per spot") + NoLegend()

lim <- range(lps.1@meta.data$nCount_Spatial, lps.2@meta.data$nCount_Spatial, lps.4@meta.data$nCount_Spatial)
leg <- get_legend(SpatialFeaturePlot(lps.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim, name = "UMIs") + theme(legend.position = "right", aspect.ratio = 1))
plot1 <- SpatialFeaturePlot(lps.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot2 <- SpatialFeaturePlot(lps.2, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot4 <- SpatialFeaturePlot(lps.4, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot_grid(plot1, plot2, plot4, leg, nrow = 1, rel_widths = c(1, 1, 1, 0.2))

lim <- range(cnt.1@meta.data$nCount_Spatial, cnt.3@meta.data$nCount_Spatial, cnt.4@meta.data$nCount_Spatial)
leg <- get_legend(SpatialFeaturePlot(cnt.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim, name = "UMIs") + theme(legend.position = "right", aspect.ratio = 1))
plot1 <- SpatialFeaturePlot(cnt.1, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot3 <- SpatialFeaturePlot(cnt.3, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot4 <- SpatialFeaturePlot(cnt.4, features = "nCount_Spatial") + scale_fill_viridis_c(option = "magma", lim = lim) + theme(legend.position = "right", aspect.ratio = 1) + NoLegend()
plot_grid(plot1, plot3, plot4, leg, nrow = 1, rel_widths = c(1, 1, 1, 0.2))

# log-normalize count data
lps.1 <- NormalizeData(lps.1, verbose = TRUE, assay = "Spatial")
dim(lps.1)
lps.2 <- NormalizeData(lps.2, verbose = TRUE, assay = "Spatial")
dim(lps.2)
lps.4 <- NormalizeData(lps.4, verbose = TRUE, assay = "Spatial")
dim(lps.4)

cnt.1 <- NormalizeData(cnt.1, verbose = TRUE, assay = "Spatial")
dim(cnt.1)
cnt.3 <- NormalizeData(cnt.3, verbose = TRUE, assay = "Spatial")
dim(cnt.3)
cnt.4 <- NormalizeData(cnt.4, verbose = TRUE, assay = "Spatial")
dim(cnt.4)

# save Seurat objects for each sample
saveRDS(lps.1, "file_path/hasel_lpsA_spatial_seurat_object.rds")
saveRDS(lps.2, "file_path/hasel_lpsB_spatial_seurat_object.rds")
saveRDS(lps.4, "file_path/hasel_lpsC_spatial_seurat_object.rds")

saveRDS(cnt.1, "file_path/hasel_cntA_spatial_seurat_object.rds")
saveRDS(cnt.3, "file_path/hasel_cntB_spatial_seurat_object.rds")
saveRDS(cnt.4, "file_path/hasel_cntC_spatial_seurat_object.rds")
