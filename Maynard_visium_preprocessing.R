# This script provides the code necessary for the preprocessing of the Maynard et al. (2021) human spatial transcriptomics data

# load required packages
library(dplyr)
library(ggplot2)
library(Seurat)
library(spatialLIBD)

### get package requirements file
pkgs <- loadedNamespaces()
for (i in 1:length(pkgs)) {
  cat(print(paste(pkgs[i], getNamespaceVersion(pkgs[i]), sep = " ")), file = "Maynard_visium_preprocessing_packages.txt", append = TRUE)
  cat("\n", file = "Maynard_visium_preprocessing_packages.txt", append = TRUE)
}

# retrieve Maynard SingleCellExperiment object from ExperimentHub
ehub <- ExperimentHub::ExperimentHub()
sce <- fetch_data(type = "sce", eh = ehub)

# save SCE object
saveRDS(sce, "file_path/spatial_LIBD_sce.rds")

# create Seurat object from Cell Ranger outputs files for all samples from Maynard et al. dataset
sample_151507 = Load10X_Spatial("file_path/Maynard_2021/151507",
                                filename = "151507_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

# assign layer identities to each spot from annotations in the retrieved SCE object
sample_151507$layer <- sce[,sce@colData$sample_name == 151507]$layer_guess_reordered

sample_151508 = Load10X_Spatial("file_path/Maynard_2021/151508",
                                filename = "151508_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151508$layer <- sce[,sce@colData$sample_name == 151508]$layer_guess_reordered

sample_151509 = Load10X_Spatial("file_path/Maynard_2021/151509",
                                filename = "151509_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151509$layer <- sce[,sce@colData$sample_name == 151509]$layer_guess_reordered

sample_151510 = Load10X_Spatial("file_path/Maynard_2021/151510",
                                filename = "151510_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151510$layer <- sce[,sce@colData$sample_name == 151510]$layer_guess_reordered

sample_151669 = Load10X_Spatial("file_path/Maynard_2021/151669",
                                filename = "151669_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151669$layer <- sce[,sce@colData$sample_name == 151669]$layer_guess_reordered

sample_151670 = Load10X_Spatial("file_path/Maynard_2021/151670",
                                filename = "151670_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151670$layer <- sce[,sce@colData$sample_name == 151670]$layer_guess_reordered

sample_151671 = Load10X_Spatial("file_path/Maynard_2021/151671",
                                filename = "151671_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151671$layer <- sce[,sce@colData$sample_name == 151671]$layer_guess_reordered

sample_151672 = Load10X_Spatial("file_path/Maynard_2021/151672",
                                filename = "151672_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151672$layer <- sce[,sce@colData$sample_name == 151672]$layer_guess_reordered

sample_151673 = Load10X_Spatial("file_path/Maynard_2021/151673",
                                filename = "151673_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151673$layer <- sce[,sce@colData$sample_name == 151673]$layer_guess_reordered

sample_151674 = Load10X_Spatial("file_path/Maynard_2021/151674",
                                filename = "151674_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151674$layer <- sce[,sce@colData$sample_name == 151674]$layer_guess_reordered

sample_151675 = Load10X_Spatial("file_path/Maynard_2021/151675",
                                filename = "151675_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151675$layer <- sce[,sce@colData$sample_name == 151675]$layer_guess_reordered

sample_151676 = Load10X_Spatial("file_path/Maynard_2021/151676",
                                filename = "151676_filtered_feature_bc_matrix.h5",
                                assay = "Spatial",
                                slice = "slice1"
)

sample_151676$layer <- sce[,sce@colData$sample_name == 151676]$layer_guess_reordered

# log-normalize count data for each sample
sample_151507 <- NormalizeData(sample_151507)
sample_151508 <- NormalizeData(sample_151508)
sample_151509 <- NormalizeData(sample_151509)
sample_151510 <- NormalizeData(sample_151510)
sample_151669 <- NormalizeData(sample_151669)
sample_151670 <- NormalizeData(sample_151670)
sample_151671 <- NormalizeData(sample_151671)
sample_151672 <- NormalizeData(sample_151672)
sample_151673 <- NormalizeData(sample_151673)
sample_151674 <- NormalizeData(sample_151674)
sample_151675 <- NormalizeData(sample_151675)
sample_151676 <- NormalizeData(sample_151676)

# save Seurat objects for each sample
saveRDS(sample_151507, "file_path/sample_151507_seuratobject.rds")
saveRDS(sample_151508, "file_path/sample_151508_seuratobject.rds")
saveRDS(sample_151509, "file_path/sample_151509_seuratobject.rds")
saveRDS(sample_151510, "file_path/sample_151510_seuratobject.rds")
saveRDS(sample_151669, "file_path/sample_151669_seuratobject.rds")
saveRDS(sample_151670, "file_path/sample_151670_seuratobject.rds")
saveRDS(sample_151671, "file_path/sample_151671_seuratobject.rds")
saveRDS(sample_151672, "file_path/sample_151672_seuratobject.rds")
saveRDS(sample_151673, "file_path/sample_151673_seuratobject.rds")
saveRDS(sample_151674, "file_path/sample_151674_seuratobject.rds")
saveRDS(sample_151675, "file_path/sample_151675_seuratobject.rds")
saveRDS(sample_151676, "file_path/sample_151676_seuratobject.rds")
