# This script provides code used for the spatial annotation of cortical regions in the mouse spatial transcriptomics data from Hasel et al. (2021)

# load required packages
library(dplyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(qs)
library(ComplexHeatmap)
library(purrr)
library(stringr)
library(HGC)
library(dendextend)
library(tibble)
library(ggridges)
library(biomaRt)
library(corrplot)

### get package requirements file
pkgs <- loadedNamespaces()
for (i in 1:length(pkgs)) {
  cat(print(paste(pkgs[i], getNamespaceVersion(pkgs[i]), sep = " ")), file = "Regional_annotation_packages.txt", append = TRUE)
  cat("\n", file = "Regional_annotation_packages.txt", append = TRUE)
}

# load Seurat objects
cnt.A <- readRDS("file_path/hasel_cntA_spatial_seurat_object.rds")
cnt.B <- readRDS("file_path/hasel_cntB_spatial_seurat_object.rds")
cnt.C <- readRDS("file_path/hasel_cntC_spatial_seurat_object.rds")

lps.A <- readRDS("file_path/hasel_lpsA_spatial_seurat_object.rds")
lps.B <- readRDS("file_path/hasel_lpsB_spatial_seurat_object.rds")
lps.C <- readRDS("file_path/hasel_lpsC_spatial_seurat_object.rds")

# add sample metadata
cnt.A$group_id <- "CNT"
cnt.B$group_id <- "CNT"
cnt.C$group_id <- "CNT"

lps.A$group_id <- "LPS"
lps.B$group_id <- "LPS"
lps.C$group_id <- "LPS"

cnt.A$sample <- "CNT_A"
cnt.B$sample <- "CNT_B"
cnt.C$sample <- "CNT_C"

lps.A$sample <- "LPS_A"
lps.B$sample <- "LPS_B"
lps.C$sample <- "LPS_C"

# merge Seurat objects into one combined Seurat object
combo <- merge(cnt.A, c(cnt.B, cnt.C, lps.A, lps.B, lps.C))

object.list <- c(cnt.A, cnt.B, cnt.C, lps.A, lps.B, lps.C)

# find common variable features 
features <- SelectIntegrationFeatures(object.list, nfeatures = 3000)

# log normalize data
combo <- NormalizeData(combo) 

# calculate QC metrics
combo[["percent.mt"]] <- PercentageFeatureSet(combo, pattern = "^mt-")

# Scale variable features and regress technical variables
combo <- ScaleData(combo, features = features, vars.to.regress = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) # regressing out UMI counts, genes detected counts, and percent.mt 

# run PCA
combo <- RunPCA(combo, features = features, npcs = 50, verbose = TRUE)

# get vector of all genes
all_genes <- Reduce(intersect, lapply(object.list, rownames))

# Integrate between samples and conditions using Harmony
set.seed(120) # make sure to set a random seed so integration results are 100% reproducible
combo <- combo %>% RunHarmony(c("sample", "group_id"), plot_convergence = TRUE, assay.use = "Spatial")

# perform UMAP using top 30 principal components 
combo <- RunUMAP(combo, reduction = "harmony", dims = 1:30)

# save Seurat object
saveRDS(combo, "combo.rds")

# identify white matter (WM) and gray matter spots
# create gene module of myelin associated genes
wm.genes <- c("Mbp", "Mobp", "Plp1", "Mag", "Mog", "Mal")
# calculate spot module scores
combo <- AddModuleScore(combo, features = list(wm.genes), name = "module_score")
combo$WM <- combo$module_score1

# label each spot exceeding threshold as WM
tmp.df <- data.frame(barcode = colnames(combo), wm = combo$WM)
tmp.df <- tmp.df %>% mutate(label = ifelse(wm > 1.480081, "WM", "other"))
combo$orig_annotation <- tmp.df$label
Idents(combo) <- combo$orig_annotation

# Isocortex gray matter spots were next manually labeled using the CellSelector function; white matter spots were adjusted using CellSelector as well
# extract seurat object of only Isocortex spots
isocortex.combo <- subset(combo, idents = "Isocortex")

# split object
iso.combo.list <- SplitObject(isocortex.combo, split.by = "sample")

# select common variable features across samples
layer.hvgs <- SelectIntegrationFeatures(iso.combo.list, nfeatures = 3000)

# integrate and cluster the isocortex seurat object and regress technical variables
isocortex.combo <- ScaleData(isocortex.combo, features = layer.hvgs, vars.to.regress = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) # regressing out UMI counts, genes detected counts, and percent.mt 

# run PCA
isocortex.combo <- RunPCA(isocortex.combo, features = layer.hvgs, npcs = 50, verbose = TRUE)

# Integrate across samples and condition using Harmony
set.seed(666) 
isocortex.combo <- isocortex.combo %>% RunHarmony(c("sample", "group_id"), plot_convergence = TRUE, assay.use = "Spatial")

# plot PCA elbow plot
tiff("layer_harmony_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(isocortex.combo, ndims = 50, reduction = "harmony")
dev.off()

tiff("layer_uncorrected_pca_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(isocortex.combo, ndims = 50, reduction = "pca") 
dev.off()

# run UMAP
isocortex.combo <- RunUMAP(isocortex.combo, reduction = "harmony", dims = 1:15)

# plot UMAP
DimPlot(isocortex.combo, reduction = "umap", group.by = "group_id")

# construct nearest neighbors graph and cluster using several resolutions
isocortex.combo <- FindNeighbors(isocortex.combo, reduction = "harmony", dims = 1:20)
isocortex.combo <- FindClusters(isocortex.combo, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

# generate color palette
library(rcartocolor)
pal <- carto_pal(nColor, "Safe")

# plot UMAP of original cluster identities using resolution of 0.6
tiff("isocortex_clustering_umap.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(isocortex.combo, reduction = "umap", group.by = "Spatial_snn_res.0.6") + theme(aspect.ratio = 1)
dev.off()

# visualize across all sections
SpatialDimPlot(isocortex.combo, group.by = "Spatial_snn_res.0.6", images = "slice1.1") 
SpatialDimPlot(isocortex.combo, group.by = "Spatial_snn_res.0.6", images = "slice1.2")
SpatialDimPlot(isocortex.combo, group.by = "Spatial_snn_res.0.6", images = "slice1.3")
SpatialDimPlot(isocortex.combo, group.by = "Spatial_snn_res.0.6", images = "slice1.4")
SpatialDimPlot(isocortex.combo, group.by = "Spatial_snn_res.0.6", images = "slice1.5")

# make cluster assignments new identifiers
Idents(isocortex.combo) <- isocortex.combo$Spatial_snn_res.0.6

# run differential expression testing; marker genes for each cluster will be used along with H&E staining to identify the cortical layer of each spot; 
iso.res0.6.markers <- FindAllMarkers(isocortex.combo)
iso.res0.6.markers <- iso.res0.6.markers %>% arrange(cluster, desc(avg_log2FC))

# assign temporary layer labels based on marker gene expression and general anatomical location of cluster; note: spots in the somotomotor and retrosplenial areas were less well defined by the original clustering; we will subcluster these later; 
current.cluster.ids <- levels(Idents(isocortex.combo))
new.cluster.ids <- c("L6",
                     "L4",
                     "L2/3",
                     "L5",
                     "L5",
                     "L6",
                     "RA L5/6",
                     "L1",
                     "RA L1",
                     "RA L2/3",
                     "RA L2/3")
Idents(isocortex.combo) <- plyr::mapvalues(x = Idents(isocortex.combo), from = current.cluster.ids, to = new.cluster.ids)
Idents(isocortex.combo) <- factor(Idents(isocortex.combo), levels = c("L1", "L2/3", "L4", "L5", "L6", "RA L1", "RA L2/3", "RA L5/6"))

# plot new cluster labels on umap
tiff("isocortex_assignments_umap.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(isocortex.combo, reduction = "umap") + scale_color_manual(values = pal) + theme(aspect.ratio = 1)
dev.off()

# evaluate across each section
SpatialDimPlot(isocortex.combo, images = "slice1") + scale_fill_manual(values = pal)
SpatialDimPlot(isocortex.combo, images = "slice1.1")+ scale_fill_manual(values = pal)
SpatialDimPlot(isocortex.combo, images = "slice1.2")+ scale_fill_manual(values = pal)
SpatialDimPlot(isocortex.combo, images = "slice1.3")+ scale_fill_manual(values = pal)
SpatialDimPlot(isocortex.combo, images = "slice1.4")+ scale_fill_manual(values = pal)
SpatialDimPlot(isocortex.combo, images = "slice1.5")+ scale_fill_manual(values = pal)

# save temporary labels as cortex_region in metadata
isocortex.combo$cortex_region <- Idents(isocortex.combo)
# also save original cortex_region annotation and cluster labels in full combo seurat object (containing spots from all regions, not just isocortex)
combo$cortex_region <- Idents(isocortex.combo)
combo$cortex_region_clustering <- isocortex.combo$Spatial_snn_res.0.6

# plot original clustering results on section in full combo seurat object
tiff("isocortex_clustering_results.tiff", units = "in", width = 6, height = 6, res = 300)
SpatialDimPlot(combo, group.by = "cortex_region_clustering", images = "slice1") + theme(aspect.ratio = 1)
dev.off()

# subset combo seurat object to contain both Isocortex and white matter spots, in new seurat object called combo.cortex
Idents(combo) <- combo$orig_annotation
combo.cortex <- subset(combo, idents = c("WM", "Isocortex"))

# make new cortex_region labels the cell identifiers
Idents(combo.cortex) <- combo.cortex$cortex_region

# plot on one section
SpatialDimPlot(combo.cortex, images = "slice1", group.by = "cortex_region") + scale_fill_manual(values = c(pal, "#FF7F00"))

# create data frame of spot barcodes, and cortex_region label 
temp.df <- data.frame(barcode = colnames(combo.cortex), cortex_region = as.character(combo.cortex$cortex_region))

# if spot is unlabeled, add "WM" white matter label
temp.df <- temp.df %>% mutate(new_region = ifelse(is.na(cortex_region), "WM", cortex_region))

# add this new identifier as cortex_region2 to the combo.cortex seurat object
combo.cortex$cortex_region2 <- temp.df$new_region

# plot cortex_region2 labels on section (showing both White matter (cortical and subcortical) and isocortex preliminary identifiers)
tiff("slice1_isocortex_region_spatialdimplot.tiff", units = "in", width = 6, height = 6, res = 300)
SpatialDimPlot(combo.cortex, images = "slice1", group.by = "cortex_region2") + scale_fill_manual(values = c(pal[1:8], "#FF7F00"), name = "Isocortex Region")
dev.off()

# save seurat objects; 
saveRDS(isocortex.combo, "isocortex_combo.rds") # this seurat object includes isocortex spots and their original regional labels based on first round of clustering 
saveRDS(combo.cortex, "combo_cortex.rds") # this seurat object includes isocortex labeled spots + (cortical + subcortical) WM spots and their regional labels based on first round of clustering

# now, we'll manually adjust the cortical layer assignments based on H&E image to fix isolated spots which were labeled with incorrect region identifiers: 
# this was done using the following code, repeatedly for each identifier label, until remaining spots appeared to be in correct location
# in this example, L1-labeled spots will be edited to assign misidentified spots to new labels

# set identifiers as spot idents
Idents(combo.cortex) <- combo.cortex$cortex_region2

# view regional identifiers for each section; note L1 labeled spots which appear out of place
SpatialDimPlot(combo.cortex, images = "slice1", group.by = "cortex_region2") 
SpatialDimPlot(combo.cortex, images = "slice1.1", group.by = "cortex_region2")
SpatialDimPlot(combo.cortex, images = "slice1.2", group.by = "cortex_region2")
SpatialDimPlot(combo.cortex, images = "slice1.3", group.by = "cortex_region2")
SpatialDimPlot(combo.cortex, images = "slice1.4", group.by = "cortex_region2")
SpatialDimPlot(combo.cortex, images = "slice1.5", group.by = "cortex_region2")

# highlight L1 spots in each section to help with visualization
cells <- WhichCells(combo.cortex, idents = "L1")
SpatialDimPlot(combo.cortex, images = "slice1", cells.highlight = cells) + NoLegend()
SpatialDimPlot(combo.cortex, images = "slice1.1", cells.highlight = cells) + NoLegend()
SpatialDimPlot(combo.cortex, images = "slice1.2", cells.highlight = cells) + NoLegend()
SpatialDimPlot(combo.cortex, images = "slice1.3", cells.highlight = cells) + NoLegend()
SpatialDimPlot(combo.cortex, images = "slice1.4", cells.highlight = cells) + NoLegend()
SpatialDimPlot(combo.cortex, images = "slice1.5", cells.highlight = cells) + NoLegend()

## select all the L1 labeled spots which should be reassigned to L2/3; these appear in slices 1.1, 1.2, & 1.5

# subset out L1 spots 
tmp.subset <- subset(combo.cortex, idents = "L1")
# view slice 1.1 in interactive cell selector and choose spots to be relabeled
p <- ggplot(tmp.subset@images$slice1.1@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
tmp.subset <- CellSelector(p, tmp.subset)
#  view slice 1.2 in interactive cell selector and choose spots to be relabeled
p <- ggplot(tmp.subset@images$slice1.2@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
# several rounds of interactive cell selection may be necessary for each slice if mislabeled spots are not contiguous; in this case, we select groups of mislabeled cells three times
tmp.subset <- CellSelector(p, tmp.subset)
tmp.subset <- CellSelector(p, tmp.subset) 
tmp.subset <- CellSelector(p, tmp.subset)
#  view slice 1.5 in interactive cell selector and choose spots to be relabeled
p <- ggplot(tmp.subset@images$slice1.5@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
tmp.subset <- CellSelector(p, tmp.subset)
# now extract vector of spot barcodes corresponding to all the selected spots; 
L1.to.L2.3 <- WhichCells(tmp.subset, idents = "SelectedCells") 

# visualize the selected spots in each section to confirm correct selection for reassignment; 
SpatialDimPlot(combo.cortex, images = "slice1.1", cells.highlight = L1.to.L2.3) + NoLegend()
SpatialDimPlot(combo.cortex, images = "slice1.2", cells.highlight = L1.to.L2.3) + NoLegend()
SpatialDimPlot(combo.cortex, images = "slice1.5", cells.highlight = L1.to.L2.3) + NoLegend()

# repeated as above for all other reassignments; omitted here

# next create a list containing all the selected spot vectors for reassignment:
cortical.corrections1 <- list(L1.to.L2.3, L1.to.L5, L1.to.L6,
                              L2.3.to.L4, L2.3.to.L5, L2.3.to.L6,
                              L4.to.L2.3, L4.to.L5, L4.to.L6,
                              L5.to.L1, L5.to.L2.3, L5.to.L4, L5.to.L6, L5.to.RA.L1, L5.to.RA.L2.3, L5.to.RA.L5.6,
                              L6.to.L5,
                              RA.L1.to.L1, RA.L1.to.L2.3, RA.L1.to.L4, RA.L1.to.L5, RA.L1.to.L6, RA.L1.to.RA.L2.3, RA.L1.to.RA.L5.6,
                              RA.L2.3.to.L4, RA.L2.3.to.L5,
                              RA.L5.6.to.L5, RA.L5.6.to.L6)

# save corrections list
saveRDS(cortical.corrections1, "cortical_corrections1_list.rds")

# create a data frame containing spot identifier labels
tmp.df <- data.frame(spot = colnames(combo.cortex), layer = as.character(combo.cortex@active.ident))

# mutate spot identifier labels to reassign spots based on vectors 
tmp.df <- tmp.df %>% mutate(layer_correction1 = ifelse(spot %in% c(L1.to.L2.3, L4.to.L2.3, L5.to.L2.3, RA.L1.to.L2.3), "L2/3",
                                                       ifelse(spot %in% c(L5.to.L1, RA.L1.to.L1), "L1",
                                                              ifelse(spot %in% c(L2.3.to.L4, L5.to.L4, RA.L1.to.L4, RA.L2.3.to.L4), "L4",
                                                                     ifelse(spot %in% c(L1.to.L5, L2.3.to.L5, L6.to.L5, RA.L1.to.L5, RA.L2.3.to.L5, RA.L5.6.to.L5), "L5",
                                                                            ifelse(spot %in% c(L1.to.L6, L2.3.to.L6, L4.to.L6, L5.to.L6, RA.L1.to.L6, RA.L5.6.to.L6), "L6",
                                                                                   ifelse(spot %in% L5.to.RA.L1, "RA L1",
                                                                                          ifelse(spot %in% c(L5.to.RA.L2.3, RA.L1.to.RA.L2.3), "RA L2/3",
                                                                                                 ifelse(spot %in% c(L5.to.RA.L5.6, RA.L1.to.RA.L5.6), "RA L5/6", layer)))))))))

# add new identifiers in metadata slot in seurat object
combo.cortex$layer_correction1 <- tmp.df$layer_correction1

# visualize corrected spot labels on section
SpatialDimPlot(combo.cortex, images = "slice1", group.by = "layer_correction1") 

# assign new labels as spot identifiers
Idents(combo.cortex) <- combo.cortex$layer_correction1

######################
# we're now going to combine the RA and other isocortex layers into one annotation; to do this, we'll subcluster the RA/somatomotor clusters to match isocortical layers 

# subset RA/somatomotor clusters from combo seurat object
ra.subset <- subset(combo.cortex, idents = c("RA L1", "RA L2/3", "RA L5/6"))
ra.subset # there are 424 spots in this area

# split seurat object by sample
ra.subset.list <- SplitObject(ra.subset, split.by = "sample")

# find top 2000 variable features
ra.subset.features <- SelectIntegrationFeatures(ra.subset.list, nfeatures = 2000)

# scale data and regress technical variables
ra.subset <- ScaleData(ra.subset, features = ra.subset.features, vars.to.regress = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"))

# run PCA
ra.subset <- RunPCA(ra.subset, features = ra.subset.features, npcs = 50, verbose = TRUE)

# integrate across sample and condition using Harmony
set.seed(576) 
ra.subset <- ra.subset %>% RunHarmony(c("sample", "group_id"), plot_convergence = TRUE, assay.use = "Spatial")

# plot PC elbow plots (harmony corrected and un-corrected)
tiff("harmony_elbowplot_ra_subset.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(combo, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot_ra_subset.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(combo, ndims = 50, reduction = "pca") 
dev.off()

# run UMAP 
ra.subset <- RunUMAP(ra.subset, reduction = "harmony", dims = 1:20)

# create nearest neighbors graph
ra.subset <- FindNeighbors(ra.subset, reduction = "harmony", dims = 1:20)

# cluster
ra.subset <- FindClusters(ra.subset, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))

# plot clusters on UMAP 
DimPlot(ra.subset, group.by = "Spatial_snn_res.0.5") 

# visualize clusters on each Visium slice; clusters appear to map well to layers
SpatialDimPlot(ra.subset, group.by = "Spatial_snn_res.0.5", images = "slice1")
SpatialDimPlot(ra.subset, group.by = "Spatial_snn_res.0.5", images = "slice1.1")
SpatialDimPlot(ra.subset, group.by = "Spatial_snn_res.0.5", images = "slice1.2")
SpatialDimPlot(ra.subset, group.by = "Spatial_snn_res.0.5", images = "slice1.3")
SpatialDimPlot(ra.subset, group.by = "Spatial_snn_res.0.5", images = "slice1.4")
SpatialDimPlot(ra.subset, group.by = "Spatial_snn_res.0.5", images = "slice1.5")

# assign cluster labels as spot idents
Idents(ra.subset) <- ra.subset$Spatial_snn_res.0.5

# save RA/somatomotor seurat object with cluster identities
saveRDS(ra.subset, "ra_subset_seurat_object.rds")

# create data frame will all cortex spots and their corrected layer annotations
tmp.df <- data.frame(spot = colnames(combo.cortex), layer = as.character(combo.cortex$layer_correction1))

# make new column assigning new subcluster layer identities to RA/somatomotor spots
tmp.df <- tmp.df %>% mutate(layer_correction2 = ifelse(spot %in% WhichCells(ra.subset, idents = 0), "L1",
                                                       ifelse(spot %in% WhichCells(ra.subset, idents = 1), "L5",
                                                              ifelse(spot %in% WhichCells(ra.subset, idents = 2), "L2/3",
                                                                     ifelse(spot %in% WhichCells(ra.subset, idents = 3), "L6", layer)))))

# add new identities to cortex seurat object
combo.cortex$layer_correction2 <- tmp.df$layer_correction2

# plot new identities on an example slice
SpatialDimPlot(combo.cortex, group.by = "layer_correction2", images = "slice1.5")

# make edits to individual spots from the RA/MO assignments which were incorrectly assigned based on clustering
# layer adjustment was performed manually as above; specific adjustments omitted here

# create list of vectors of each spot being reassigned
cortical.corrections2 <- list(final.L1.to.L2.3, final.L1.to.L5,
                              final.L5.to.L1, final.L5.to.L2.3, final.L5.to.L6,
                              final.L6.to.L1, final.L6.to.L5)

# save correction list
saveRDS(cortical.corrections2, "cortical_corrections2_list.rds")

# create data frame with spot identities
tmp.df <- data.frame(spot = colnames(combo.cortex), layer = as.character(combo.cortex@active.ident))

# create new column renaming mislabeled spots to their new layer labels
tmp.df <- tmp.df %>% mutate(layer_correction3 = ifelse(spot %in% c(final.L1.to.L2.3, final.L5.to.L2.3), "L2/3",
                                                       ifelse(spot %in% c(final.L5.to.L1, final.L6.to.L1), "L1",
                                                              ifelse(spot %in% c(final.L1.to.L5, final.L6.to.L5), "L5",
                                                                     ifelse(spot %in% final.L5.to.L6, "L6", layer)))))

# assign corrected spot labels to seurat object
combo.cortex$layer_correction3 <- tmp.df$layer_correction3

# set identities to new corrected layer annotations
Idents(combo.cortex) <- combo.cortex$layer_correction3

# now identify subcortical white matter spots from each slice so they can be removed from final cortex object
subcortex.wm.subset <- combo.cortex
p <- ggplot(subcortex.wm.subset@images$slice1@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
p <- ggplot(subcortex.wm.subset@images$slice1.1@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
p <- ggplot(subcortex.wm.subset@images$slice1.2@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
p <- ggplot(subcortex.wm.subset@images$slice1.3@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
p <- ggplot(subcortex.wm.subset@images$slice1.4@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
p <- ggplot(subcortex.wm.subset@images$slice1.5@coordinates, aes(x = imagerow, y = imagecol)) +
  geom_point()
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)
subcortex.wm.subset <- CellSelector(p, subcortex.wm.subset)

# create vector of spots
subcortical.wm <- WhichCells(subcortex.wm.subset, idents = "SelectedCells")

# create data frame of layer annotations and barcodes
tmp.df <- data.frame(spot = colnames(combo.cortex), layer_annotation = as.character(combo.cortex$layer_correction3))

# if barcode is in subcortical WM vector, change identity to NA
tmp.df <- tmp.df %>% mutate(final_layer_annotation = ifelse(spot %in% subcortical.wm, NA, layer_annotation))

# apply new layer annotations to cortex seurat object
combo.cortex$final_layer_annotation <- tmp.df$final_layer_annotation

# create vector of spots labeled NA as subcortical WM 
other.na <- setdiff(colnames(combo), colnames(combo.cortex))

# create new data frame of spot barcodes and final annotations
tmp.df <- data.frame(spot = c(colnames(combo.cortex), other.na), layer_annotation = c(as.character(combo.cortex$layer_correction3), rep(NA, length(other.na))))
tmp.df <- tmp.df %>% mutate(final_layer_annotation = ifelse(spot %in% subcortical.wm, NA, layer_annotation))

# reorder to match barcode order in combo seurat object
tmp.df <- tmp.df[match(colnames(combo), tmp.df$spot),]

# add final annotations to original combo seurat object
combo$final_layer_annotation <- tmp.df$final_layer_annotation

# get color palette
new.pal <- ggsci::pal_npg()(5)

# plot final layer annotations on representative sections
postscript("cnt_A_final_layer_annotation.ps", width = 6, height = 6)
SpatialDimPlot(combo, images = "slice1", group.by = "final_layer_annotation") + scale_fill_manual(values = c(new.pal, "#FF7F00"), name = "Layer") + theme(aspect.ratio = 1) + NoLegend()
dev.off()

postscript("lps_A_final_layer_annotation.ps", width = 6, height = 6)
SpatialDimPlot(combo, images = "slice1.3", group.by = "final_layer_annotation") + scale_fill_manual(values = c(new.pal, "#FF7F00"), name = "Layer") + theme(aspect.ratio = 1) + NoLegend()
dev.off()

# visualize on all sections
SpatialDimPlot(combo, images = "slice1.1", group.by = "final_layer_annotation") + scale_fill_manual(values = c(new.pal, "#FF7F00"), name = "Layer")
SpatialDimPlot(combo, images = "slice1.2", group.by = "final_layer_annotation") + scale_fill_manual(values = c(new.pal, "#FF7F00"), name = "Layer")
SpatialDimPlot(combo, images = "slice1.3", group.by = "final_layer_annotation") + scale_fill_manual(values = c(new.pal, "#FF7F00"), name = "Layer")
SpatialDimPlot(combo, images = "slice1.4", group.by = "final_layer_annotation") + scale_fill_manual(values = c(new.pal, "#FF7F00"), name = "Layer")
SpatialDimPlot(combo, images = "slice1.5", group.by = "final_layer_annotation") + scale_fill_manual(values = c(new.pal, "#FF7F00"), name = "Layer")

# plot representative slice with cortical spots highlighted 
postscript("combo_cortex_highlight.ps", width = 4, height = 4)
SpatialDimPlot(combo, group.by = "final_layer_annotation", images = "slice1", 
               cols = c("red", "red", "red", "red", "red", "red")) + 
  theme(aspect.ratio = 1) + labs(fill = "") + NoLegend()
dev.off()

# save cortex and original seurat objects with final annotations
saveRDS(combo.cortex, "combo_cortex_seurat_object_final.rds")
saveRDS(combo, "combo_seurat_object_final.rds")

# now, let's plot original WM thresholding values for visualization
# for all cortical spots, plot white matter gene module scores; 
tmp.df <- data.frame(barcode = colnames(combo), wm = combo$WM, region = combo$final_layer_annotation)
tmp.df <- tmp.df %>% mutate("wm2" = ifelse(!is.na(region), wm, NA))
combo$WM2 <- tmp.df$wm2

postscript("combo_white_matter_module_onlycortex.ps", width = 4, height = 4)
SpatialFeaturePlot(combo, features = "WM2", images = "slice1") + scale_fill_distiller(palette = "RdBu") + theme(aspect.ratio = 1) + NoLegend()
dev.off()

# plot legend object
leg <- as_ggplot(get_legend(SpatialFeaturePlot(combo, features = "WM2", images = "slice1") + scale_fill_distiller(palette = "RdBu", lim = range(combo$WM2, na.rm = TRUE)) + theme(aspect.ratio = 1)))
postscript("combo_white_matter_module_onlycortex_legend.ps", width = 4, height = 4)
leg
dev.off()

# remove non-cortical spots; 
tmp.df.trimmed <- tmp.df %>% filter(!is.na(region))

# plot white matter module scores for cortical spots as a density plot
tmp.tmp.df <- data.frame(WM2 = combo$WM2, sample = rep("X", ncol(combo)))
postscript("combo_whitematter_cortexonly_densityplot.ps", width = 4, height = 4)
ggplot(tmp.tmp.df, aes(x = `WM2`, y = `sample`, fill = stat(x))) +
  scale_x_continuous(limits = range(combo$WM2, na.rm = TRUE), expand = c(0, 0)) + 
  geom_density_ridges_gradient(scale = 50, rel_min_height = 0.01) +
  scale_fill_distiller(palette = "RdBu", name = "White Matter Score", lim = range(combo$WM2, na.rm = TRUE),
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
  labs(x = "White Matter Score") + 
  theme(legend.position = "top", plot.background = element_blank(), panel.background = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + NoLegend() + 
  geom_vline(xintercept = 1.480081, linetype = "longdash", size = 1)
dev.off()

# add final layer annotation metadata to isocortex.combo seurat object
isocortex.combo$final_layer_annotation <- combo$final_layer_annotation

# make color palette
new.pal <- ggsci::pal_npg()(5)

# plot final layer annotations on UMAP (excluding white matter spots)
postscript("final_layer_annotation_umap_all_isocortex_no_wm.ps", width = 4, height = 4)
DimPlot(isocortex.combo, group.by = "final_layer_annotation", cols = new.pal) + theme(aspect.ratio = 1) + ggtitle(NULL) + NoLegend()
dev.off()

# plot final layer annotation legend
postscript("final_layer_annotation_legend.ps", width = 8, height = 4)
as_ggplot(get_legend(DimPlot(combo, group.by = "final_layer_annotation", cols = c(new.pal, "#FF7F00"), pt.size = 8) + theme(aspect.ratio = 1, legend.position = "top") +
                       guides(color = guide_legend(nrow = 1))))
dev.off()


# to evaluate the accuracy of our clustering and our ability to compare between the human and mouse regions annotated, we will evaluate the correlation between expression of highly variable genes in the human and mouse Visium data across each annotated region

# extract mouse visium highly variable genes
mouse.hvgs <- features
# extract highly variable genes from Maynard et al Visium data
human.hvgs <- SelectIntegrationFeatures(SplitObject(maynard.combo, split.by = "sample"), nfeatures = 3000)

# get ensembl mart objects
human <- useEnsembl(biomart = "ensembl", 
                    dataset = "hsapiens_gene_ensembl", 
                    mirror = "useast")
mouse <- useEnsembl(biomart = "ensembl", 
                    dataset = "mmusculus_gene_ensembl", 
                    mirror = "useast")

# this function converts human to mouse gene names; retrieves one-to-one or one-to-many orthologs
convertHumanGeneList <- function(x){
  require(biomaRt)
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T, )
  return(genesV2)
}

# get orthologs of human HVGs
human.hvgs.df <- convertHumanGeneList(human.hvgs)
colnames(human.hvgs.df) <- c("gene", "mus_gene")

# evaluate how many genes had successfully identified mouse orthologs
colSums(!is.na(human.hvgs.df)) # found mouse orthologs for 2,463/3,000 genes

human.hvgs.df <- human.hvgs.df[!(duplicated(human.hvgs.df$gene) | duplicated(human.hvgs.df$gene, fromLast = TRUE)), ] # remove all many to one orthologs
human.hvgs.df <- human.hvgs.df[!(duplicated(human.hvgs.df$mus_gene) | duplicated(human.hvgs.df$mus_gene, fromLast = TRUE)), ] # remove all one to many orthlogs
human.hvgs.orthologs <- human.hvgs.df$mus_gene
human.hvgs.orthologs <- human.hvgs.orthologs[!is.na(human.hvgs.orthologs)]

# select orthologous genes for correlation that were highly variable across both the mouse and human data
common.hvgs <- intersect(human.hvgs.orthologs, mouse.hvgs) 

length(common.hvgs) # 707 genes were common between datasets as highly variable, with one-to-one orthologs

hvg.ortholog.df <- human.hvgs.df[human.hvgs.df$mus_gene %in% common.hvgs,]

# get average expression matrices for spatial spots from both datasets 
human.avg.exp <- as.matrix(as.data.frame(AverageExpression(maynard.combo, group.by = "final_layer")))
mouse.avg.exp <- as.matrix(as.data.frame(AverageExpression(combo, group.by = "final_layer_annotation")))

# extract only the genes of interest
human.avg.exp <- human.avg.exp[hvg.ortholog.df$gene,]
mouse.avg.exp <- mouse.avg.exp[hvg.ortholog.df$mus_gene,]

# rename columns to reflect regions
colnames(human.avg.exp) <- c("L1", "L2/3", "L4", "L5", "L6", "WM")
colnames(mouse.avg.exp) <- c("L1", "L2/3", "L4", "L5", "L6", "WM")

# order the expression matrices to match between datasets
human.avg.exp <- human.avg.exp[match(hvg.ortholog.df$gene, rownames(human.avg.exp)),]
mouse.avg.exp <- mouse.avg.exp[match(hvg.ortholog.df$mus_gene, rownames(mouse.avg.exp)),]

# z-score
human.avg.exp <- t(scale(t(human.avg.exp)))
mouse.avg.exp <- t(scale(t(mouse.avg.exp)))

# run Spearman correlation
layer.hvg.cor.results <- cor(human.avg.exp, mouse.avg.exp, method = "spearman") # output: human in rows, mouse in columns

# plot correlations between human and mouse regions 
postscript("layer_hvg_corrplot.ps", width = 4, height = 4)
corrplot(layer.hvg.cor.results, col = rev(COL2("RdBu", 200)), tl.col = 'black', addgrid.col = NA)
dev.off()

# plot legend for correlation plot
pdf("layer_hvg_corrplot.pdf", width = 4, height = 4)
corrplot(layer.hvg.cor.results, col = rev(COL2("RdBu", 200)), tl.col = 'black', addgrid.col = NA)
dev.off()
