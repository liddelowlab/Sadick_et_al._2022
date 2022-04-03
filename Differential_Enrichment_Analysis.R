# This script provides the code necessary for the differential enrichment testing of astrocyte cluster gene signature modules based on our single-nucleus RNA-seq results between layers and conditions from the Hasel et al. (2021) and Maynard et al. (2021) spatial transcriptomics datasets  

# load required packages
library(dplyr)
library(ggplot2)
library(Seurat)
library(dendextend)
library(tibble)
library(ggdendro)
library(tidyr)
library(cowplot)
library(ComplexHeatmap)
library(ggpubr)

### get package requirements file
pkgs <- loadedNamespaces()
for (i in 1:length(pkgs)) {
  cat(print(paste(pkgs[i], getNamespaceVersion(pkgs[i]), sep = " ")), file = "Differential_Enrichment_Analysis_packages.txt", append = TRUE)
  cat("\n", file = "Differential_Enrichment_Analysis_packages.txt", append = TRUE)
}

# read in Hasel seurat object
combo <- readRDS("file_path/combo_seurat_object_final.rds")

# read in mouse gene modules
mouse.modules <- readRDS("file_path/final_mouse_gene_modules_forhaselsections.rds")

# gene modules
combo <- AddModuleScore(combo, features = mouse.modules, name = "Sadick_mouse_gene_modules")

# create data frame of cluster gene module scores, layer identities, sample identities, and conditions
test = combo[[c("Sadick_mouse_gene_modules1", "Sadick_mouse_gene_modules2", "Sadick_mouse_gene_modules3", "Sadick_mouse_gene_modules4", 
                "Sadick_mouse_gene_modules5", "Sadick_mouse_gene_modules6", "Sadick_mouse_gene_modules7", "Sadick_mouse_gene_modules8", 
                "Sadick_mouse_gene_modules9", 
                "Sadick_mouse_gene_modules10", "Sadick_mouse_gene_modules11", "Sadick_mouse_gene_modules12", "Sadick_mouse_gene_modules13",
                "Sadick_mouse_gene_modules14", "Sadick_mouse_gene_modules15", "Sadick_mouse_gene_modules16", "Sadick_mouse_gene_modules17",
                "Sadick_mouse_gene_modules18", "final_layer_annotation", "sample", "group_id")]]

# remove spots that are not cortical
test <- test %>% filter(!is.na(final_layer_annotation))

# spot cluster module score plots vs cortical layer (tiffs)
tiff("cluster0_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules1, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster1_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules2, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster2_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules3, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster3_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules4, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster4_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules5, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster5_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules6, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster6_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules7, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster7_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules8, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("cluster8_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules9, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

# same plots as postscript files
postscript("cluster0_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules1, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster1_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules2, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster2_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules3, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster3_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules4, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster4_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules5, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster5_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules6, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster6_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules7, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster7_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules8, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("cluster8_gene_set_score.ps", width = 5, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules9, fill = final_layer_annotation)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=final_layer_annotation), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

# spot AD-associated gene module scores by cortical layer, split by condition (tiffs)
tiff("cluster0_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules10, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster1_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules11, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster2_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules12, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster3_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules13, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster4_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules14, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster5_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules15, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster6_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules16, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster7_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules17, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster8_ad_gene_set_score.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules18, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# postscript files
postscript("cluster0_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules10, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster1_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules11, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster2_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules12, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster3_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules13, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster4_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules14, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster5_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules15, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster6_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules16, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster7_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules17, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster8_ad_gene_set_score.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules18, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) + 
  ggnewscale::new_scale("fill") + 
  ggnewscale::new_scale("color") + 
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) + 
  scale_fill_manual(values = c("navyblue", "red3")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8 AD Genes") + theme_pubr() + NoLegend() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# cluster gene module scores split by condition by cortical layer (tiffs)
tiff("cluster0_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules1, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster1_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules2, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster2_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules3, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster3_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules4, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster4_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules5, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


tiff("cluster5_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules6, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster6_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules7, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster7_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules8, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

tiff("cluster8_gene_set_score_split_by_condition.tiff", units = "in", width = 8, height = 4, res = 300)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules9, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

# postscript files
postscript("cluster0_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules1, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster1_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules2, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster2_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules3, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster3_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules4, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster4_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules5, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


postscript("cluster5_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules6, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster6_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules7, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster7_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules8, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()

postscript("cluster8_gene_set_score_split_by_condition.ps",  width = 8, height = 4)
ggplot(test, aes(x = final_layer_annotation, y = Sadick_mouse_gene_modules9, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


# cluster 3 genes split by condition; coord flipped 
postscript("cluster3_gene_set_score_split_by_condition_coordflip.ps",  width = 4, height = 8)
ggplot(test, aes(x = forcats::fct_rev(final_layer_annotation), y = Sadick_mouse_gene_modules4, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3 Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + coord_flip()
dev.off()

# cluster 5 AD genes split by condition; coord flipped 
postscript("cluster5_ad_gene_set_score_split_by_condition_coordflip.ps",  width = 4, height = 8)
ggplot(test, aes(x = forcats::fct_rev(final_layer_annotation), y = Sadick_mouse_gene_modules15, color = group_id, fill = final_layer_annotation)) +
  geom_boxplot(aes(fill = final_layer_annotation, color = group_id), notch  = TRUE, width = .35, outlier.shape = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + scale_color_manual(values = c("navyblue", "red3")) +
  ggnewscale::new_scale("fill") +
  ggnewscale::new_scale("color") +
  ggdist::stat_halfeye(aes(color = "black", fill=group_id), adjust = .5, width = .35, .width = 0, justification = -0.75, point_colour = NA, position=position_dodge(0.75)) +
  scale_fill_manual(values = c("navyblue", "red3")) +
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5 AD Genes") + theme_pubr() + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + coord_flip()
dev.off()

# Statistical testing: 

# Kruskal-Wallis tests for each Cluster Module across layers: 
kruskal.test.results.mouse <- data.frame(statistic = c(), df = c(), p_value = c(), method = c(), data = c())
count <- 0
for(x in c(paste0("Sadick_mouse_gene_modules", 1:9))){
  count <- count + 1
  print(count)
  tmp <- kruskal.test(as.formula(paste0(x, " ~ final_layer_annotation")), data = test)
  tmp.df <- data.frame(statistic = tmp$statistic, df = tmp$parameter, p_value = tmp$p.value, method = tmp$method, data = tmp$data.name)
  kruskal.test.results.mouse <- rbind(kruskal.test.results.mouse, tmp.df)
}
kruskal.test.results.mouse
rownames(kruskal.test.results.mouse) <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8")
kruskal.test.results.mouse <- kruskal.test.results.mouse %>% mutate(adj_p_val = p_value*9) # apply Bonferroni correction

saveRDS(kruskal.test.results.mouse, "kruskal_wallis_results_mouse.rds")
write.csv(kruskal.test.results.mouse, "kruskal_wallis_results_mouse.csv")

## wilcoxon tests for each layer vs all other layers pooled; then correct for multiple comparisons using Bonferonni correction
# bonferroni correction: number of cluster modules tested: 9; number of comparisons 6; 9*6 = 54
wilcoxon.test.results.mouse <- data.frame(region = c(), cluster_module = c(), statistic = c(), p_value = c(), alternative = c(), method = c(), conf_int = c(), estimate = c(), cluster_average = c(), other_average = c())
count <- 0
for(x in colnames(test)[1:9]){ # looping through first 9 columns of test (cluster modules scores 0-9)
  for(y in c("L1", "L2/3", "L4", "L5", "L6", "WM")){
    count <- count + 1
    print(count)
    tmp <- wilcox.test(test %>% filter(final_layer_annotation == y) %>% .[,x], test %>% filter(final_layer_annotation != y) %>% .[,x], alternative = "two.sided", conf.int = TRUE)
    tmp.df <- data.frame(region = y, cluster_module = x, statistic = tmp$statistic, p_value = tmp$p.value, alternative = tmp$alternative, method = tmp$method,
                         conf_int_lower = tmp$conf.int[1],
                         conf_int_upper = tmp$conf.int[2],
                         estimate = tmp$estimate,
                         cluster_average = test %>% filter(final_layer_annotation == y) %>% .[,x] %>% mean(.), 
                         other_average = test %>% filter(final_layer_annotation != y) %>% .[,x] %>% mean(.))
    wilcoxon.test.results.mouse <- rbind(wilcoxon.test.results.mouse, tmp.df)
  }
}
wilcoxon.test.results.mouse
wilcoxon.test.results.mouse <- wilcoxon.test.results.mouse %>% mutate(adj_p_val = p_value*54) # apply Bonferroni correction
wilcoxon.test.results.mouse <- wilcoxon.test.results.mouse %>% mutate(Cluster = as.numeric(gsub("Sadick_mouse_gene_modules", "", cluster_module))-1)
wilcoxon.test.results.mouse <- wilcoxon.test.results.mouse %>% mutate(enrich = ifelse(adj_p_val >= 0.05, "",
                                                                                      ifelse(estimate > 0, "+", "-")))
saveRDS(wilcoxon.test.results.mouse, "wilcoxon_test_results_mouse.rds")
write.csv(wilcoxon.test.results.mouse, "wilcoxon_test_results_mouse.csv")

wilcoxon.mouse.summary <- wilcoxon.test.results.mouse[,c("Cluster", "region", "enrich")]
saveRDS(wilcoxon.mouse.summary, "wilcoxon_mouse_summary.rds")
write.csv(wilcoxon.mouse.summary, "wilcoxon_mouse_summary.csv")

## wilcoxon tests for each cluster gene module; control versus LPS; (note: these are cluster, not AD, gene modules)
wilcoxon.test.results.mouse.lps.v.control <- data.frame(region = c(), cluster_module = c(), statistic = c(), p_value = c(), alternative = c(), method = c(), conf_int = c(), estimate = c(), lps_average = c(), cnt_average = c())
count <- 0
for(x in colnames(test)[1:9]){ 
  for(y in c("L1", "L2/3", "L4", "L5", "L6", "WM")){
    count <- count + 1
    print(count)
    tmp <- wilcox.test(test %>% filter(group_id == "LPS") %>% filter(final_layer_annotation == y) %>% .[,x], test %>% filter(group_id == "CNT") %>% filter(final_layer_annotation == y) %>% .[,x], alternative = "two.sided", conf.int = TRUE)
    tmp.df <- data.frame(region = y, cluster_module = x, statistic = tmp$statistic, p_value = tmp$p.value, alternative = tmp$alternative, method = tmp$method,
                         conf_int_lower = tmp$conf.int[1],
                         conf_int_upper = tmp$conf.int[2],
                         estimate = tmp$estimate,
                         lps_average = test %>% filter(group_id == "LPS") %>% filter(final_layer_annotation == y) %>% .[,x] %>% mean(.), 
                         cnt_average = test %>% filter(group_id == "CNT") %>% filter(final_layer_annotation == y) %>% .[,x] %>% mean(.))
    wilcoxon.test.results.mouse.lps.v.control <- rbind(wilcoxon.test.results.mouse.lps.v.control, tmp.df)
  }
}
wilcoxon.test.results.mouse.lps.v.control <- wilcoxon.test.results.mouse.lps.v.control %>% mutate(adj_p_val = p_value*54) # apply Bonferroni correction
wilcoxon.test.results.mouse.lps.v.control <- wilcoxon.test.results.mouse.lps.v.control %>% mutate(Cluster = as.numeric(gsub("Sadick_mouse_gene_modules", "", cluster_module))-1)
wilcoxon.test.results.mouse.lps.v.control <- wilcoxon.test.results.mouse.lps.v.control %>% mutate(enrich = ifelse(adj_p_val >= 0.05, "",
                                                                                                                  ifelse(estimate > 0, "+", "-")))
saveRDS(wilcoxon.test.results.mouse.lps.v.control, "wilcoxon_test_results_mouse_lps_v_control.rds")
write.csv(wilcoxon.test.results.mouse.lps.v.control, "wilcoxon_test_results_mouse_lps_v_control.csv")

wilcoxon.mouse.summary.lps.v.control <- wilcoxon.test.results.mouse.lps.v.control[,c("Cluster", "region", "enrich")]
saveRDS(wilcoxon.mouse.summary.lps.v.control, "wilcoxon_mouse_summary_lps_v_control.rds")
write.csv(wilcoxon.mouse.summary.lps.v.control, "wilcoxon_mouse_summary_lps_v_control.csv")
wilcoxon.mouse.summary.lps.v.control

## wilcoxon tests for each cluster gene module; control versus LPS; AD GENE MODULES
wilcoxon.test.results.mouse.lps.v.control.ad.genes <- data.frame(region = c(), cluster_module = c(), statistic = c(), p_value = c(), alternative = c(), method = c(), conf_int = c(), estimate = c(), lps_average = c(), cnt_average = c())
count <- 0
for(x in colnames(test)[10:18]){ 
  for(y in c("L1", "L2/3", "L4", "L5", "L6", "WM")){
    count <- count + 1
    print(count)
    tmp <- wilcox.test(test %>% filter(group_id == "LPS") %>% filter(final_layer_annotation == y) %>% .[,x], test %>% filter(group_id == "CNT") %>% filter(final_layer_annotation == y) %>% .[,x], alternative = "two.sided", conf.int = TRUE)
    tmp.df <- data.frame(region = y, cluster_module = x, statistic = tmp$statistic, p_value = tmp$p.value, alternative = tmp$alternative, method = tmp$method,
                         conf_int_lower = tmp$conf.int[1],
                         conf_int_upper = tmp$conf.int[2],
                         estimate = tmp$estimate,
                         lps_average = test %>% filter(group_id == "LPS") %>% filter(final_layer_annotation == y) %>% .[,x] %>% mean(.), 
                         cnt_average = test %>% filter(group_id == "CNT") %>% filter(final_layer_annotation == y) %>% .[,x] %>% mean(.))
    wilcoxon.test.results.mouse.lps.v.control.ad.genes <- rbind(wilcoxon.test.results.mouse.lps.v.control.ad.genes, tmp.df)
  }
}
wilcoxon.test.results.mouse.lps.v.control.ad.genes
wilcoxon.test.results.mouse.lps.v.control.ad.genes <- wilcoxon.test.results.mouse.lps.v.control.ad.genes %>% mutate(adj_p_val = p_value*54) # apply Bonferroni correction
wilcoxon.test.results.mouse.lps.v.control.ad.genes <- wilcoxon.test.results.mouse.lps.v.control.ad.genes %>% mutate(Cluster = as.numeric(gsub("Sadick_mouse_gene_modules", "", cluster_module))-15)
wilcoxon.test.results.mouse.lps.v.control.ad.genes <- wilcoxon.test.results.mouse.lps.v.control.ad.genes %>% mutate(enrich = ifelse(adj_p_val >= 0.05, "",
                                                                                                                                    ifelse(estimate > 0, "+", "-")))
saveRDS(wilcoxon.test.results.mouse.lps.v.control.ad.genes, "wilcoxon_test_results_mouse_lps_v_control_ad_genes.rds")
write.csv(wilcoxon.test.results.mouse.lps.v.control.ad.genes, "wilcoxon_test_results_mouse_lps_v_control_ad_genes.csv")

wilcoxon.mouse.summary.lps.v.control.ad.genes <- wilcoxon.test.results.mouse.lps.v.control.ad.genes[,c("Cluster", "region", "enrich")]
saveRDS(wilcoxon.mouse.summary.lps.v.control.ad.genes, "wilcoxon_mouse_summary_lps_v_control_ad_genes.rds")
write.csv(wilcoxon.mouse.summary.lps.v.control.ad.genes, "wilcoxon_mouse_summary_lps_v_control_ad_genes.csv")
wilcoxon.mouse.summary.lps.v.control.ad.genes

#### Maynard, human:
# load seurat objects
sample_151507 <- readRDS("sample_151507_seuratobject.rds")
sample_151507$sample <- "151507"
sample_151508 <- readRDS("sample_151508_seuratobject.rds")
sample_151508$sample <- "151508"
sample_151509 <- readRDS("sample_151509_seuratobject.rds")
sample_151509$sample <- "151509"
sample_151510 <- readRDS("sample_151510_seuratobject.rds")
sample_151510$sample <- "151510"
sample_151669 <- readRDS("sample_151669_seuratobject.rds")
sample_151669$sample <- "151669"
sample_151670 <- readRDS("sample_151670_seuratobject.rds")
sample_151670$sample <- "151670"
sample_151671 <- readRDS("sample_151671_seuratobject.rds")
sample_151671$sample <- "151671"
sample_151672 <- readRDS("sample_151672_seuratobject.rds")
sample_151672$sample <- "151672"
sample_151673 <- readRDS("sample_151673_seuratobject.rds")
sample_151673$sample <- "151673"
sample_151674 <- readRDS("sample_151674_seuratobject.rds")
sample_151674$sample <- "151674"
sample_151675 <- readRDS("sample_151675_seuratobject.rds")
sample_151675$sample <- "151675"
sample_151676 <- readRDS("sample_151676_seuratobject.rds")
sample_151676$sample <- "151676"

# combine into one object
maynard.combo <- merge(sample_151507, c(sample_151508, 
                                        sample_151509,
                                        sample_151510,
                                        sample_151669,
                                        sample_151670,
                                        sample_151671,
                                        sample_151672,
                                        sample_151673,
                                        sample_151674,
                                        sample_151675,
                                        sample_151676))

rm(list = c("sample_151507", "sample_151508", 
            "sample_151509",
            "sample_151510",
            "sample_151669",
            "sample_151670",
            "sample_151671",
            "sample_151672",
            "sample_151673",
            "sample_151674",
            "sample_151675",
            "sample_151676"))
gc()

# load human gene modules
human.modules <- readRDS("file_path/final_human_gene_modules_formaynardsections.rds")

# calculate cluster module scores
maynard.combo <- AddModuleScore(maynard.combo, features = human.modules, name = "Sadick_human_gene_module")

# save Maynard combo seurat object
saveRDS(maynard.combo, "maynard_combo_seurat_object.rds")

# create data frame of human layer module scores, layer identities, and sample ids
human.layer.modules.df = maynard.combo[[c("Sadick_human_gene_module1", "Sadick_human_gene_module2", "Sadick_human_gene_module3", "Sadick_human_gene_module4", 
                                          "Sadick_human_gene_module5", "Sadick_human_gene_module6", "Sadick_human_gene_module7", "Sadick_human_gene_module8", 
                                          "Sadick_human_gene_module9", "layer", "sample")]]

human.layer.modules.df <- human.layer.modules.df %>% filter(!is.na(layer)) # getting rid of spots labeled NA for layer
human.layer.modules.df <- human.layer.modules.df %>% mutate(layer = ifelse(layer == "Layer2" | layer == "Layer3", "L2/3", # grouping 2 and 3 together (to match mouse) and renaming to match mouse
                                                                           ifelse(layer == "Layer1", "L1", 
                                                                                  ifelse(layer == "Layer4", "L4",
                                                                                         ifelse(layer == "Layer5", "L5",
                                                                                                ifelse(layer == "Layer6", "L6", layer))))))
# create color palette
new.pal <- ggsci::pal_npg()(5)

# gene module score plots; tiffs
tiff("maynard_cluster0_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module1, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster1_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module2, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster2_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module3, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster3_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module4, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster4_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module5, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster5_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module6, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster6_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module7, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster7_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module8, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

tiff("maynard_cluster8_gene_set_score.tiff", units = "in", width = 5, height = 4, res = 300)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module9, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()


# postscript files
postscript("maynard_cluster0_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module1, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 0") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster1_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module2, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 1") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster2_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module3, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 2") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster3_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module4, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 3") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster4_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module5, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 4") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster5_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module6, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 5") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster6_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module7, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 6") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster7_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module8, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 7") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

postscript("maynard_cluster8_gene_set_score.ps", width = 5, height = 4)
ggplot(human.layer.modules.df, aes(x = layer, y = Sadick_human_gene_module9, fill = layer)) +
  geom_boxplot(notch=TRUE, width=0.4, outlier.shape = NA) + 
  ggdist::stat_halfeye(aes(color = "black", fill=layer), adjust = .5, width = .4, .width = 0, justification = -0.7, point_colour = NA) + 
  scale_fill_manual(values = c(new.pal, "#FF7F00")) + 
  labs(x = NULL, y = "Gene Module Score", title = "Cluster 8") + theme_pubr() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + NoLegend()
dev.off()

# statistical testing

# Kruskal-Wallis tests for each Cluster Module across layers:
kruskal.test.results.human <- data.frame(statistic = c(), df = c(), p_value = c(), method = c(), data = c())
count <- 0
for(x in c(paste0("Sadick_human_gene_module", 1:9))){
  count <- count + 1
  print(count)
  tmp <- kruskal.test(as.formula(paste0(x, " ~ layer")), data = human.layer.modules.df)
  tmp.df <- data.frame(statistic = tmp$statistic, df = tmp$parameter, p_value = tmp$p.value, method = tmp$method, data = tmp$data.name)
  kruskal.test.results.human <- rbind(kruskal.test.results.human, tmp.df)
}
kruskal.test.results.human
rownames(kruskal.test.results.human) <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8")
kruskal.test.results.human <- kruskal.test.results.human %>% mutate(adj_p_val = p_value*9) # apply Bonferroni correction

saveRDS(kruskal.test.results.human, "kruskal_wallis_results_human.rds")
write.csv(kruskal.test.results.human, "kruskal_wallis_results_human.csv")

## wilcoxon tests for each layer vs all other layers pooled; then correct for multiple comparisons using Bonferonni correction
# bonferroni correction: number of cluster modules tested: 9; number of comparisons 6; 9*6 = 54
wilcoxon.test.results.human <- data.frame(region = c(), cluster_module = c(), statistic = c(), p_value = c(), alternative = c(), method = c(), conf_int = c(), estimate = c(), cluster_average = c(), other_average = c())
count <- 0
for(x in colnames(human.layer.modules.df)[1:9]){ 
  for(y in c("L1", "L2/3", "L4", "L5", "L6", "WM")){
    count <- count + 1
    print(count)
    tmp <- wilcox.test(human.layer.modules.df %>% filter(layer == y) %>% .[,x], human.layer.modules.df %>% filter(layer != y) %>% .[,x], alternative = "two.sided", conf.int = TRUE)
    tmp.df <- data.frame(region = y, cluster_module = x, statistic = tmp$statistic, p_value = tmp$p.value, alternative = tmp$alternative, method = tmp$method,
                         conf_int_lower = tmp$conf.int[1],
                         conf_int_upper = tmp$conf.int[2],
                         estimate = tmp$estimate,
                         cluster_average = human.layer.modules.df %>% filter(layer == y) %>% .[,x] %>% mean(.), 
                         other_average = human.layer.modules.df %>% filter(layer != y) %>% .[,x] %>% mean(.))
    wilcoxon.test.results.human <- rbind(wilcoxon.test.results.human, tmp.df)
  }
}
wilcoxon.test.results.human
wilcoxon.test.results.human <- wilcoxon.test.results.human %>% mutate(adj_p_val = p_value*54) # apply Bonferroni correction
wilcoxon.test.results.human <- wilcoxon.test.results.human %>% mutate(Cluster = as.numeric(gsub("Sadick_human_gene_module", "", cluster_module))-1)
wilcoxon.test.results.human <- wilcoxon.test.results.human %>% mutate(enrich = ifelse(adj_p_val >= 0.05, "",
                                                                                      ifelse(estimate > 0, "+", "-")))
saveRDS(wilcoxon.test.results.human, "wilcoxon_test_results_human.rds")
write.csv(wilcoxon.test.results.human, "wilcoxon_test_results_human.csv")

wilcoxon.human.summary <- wilcoxon.test.results.human[,c("Cluster", "region", "enrich")]
saveRDS(wilcoxon.human.summary, "wilcoxon_human_summary.rds")
write.csv(wilcoxon.human.summary, "wilcoxon_human_summary.csv")

# create differential enrichment summary plots

# layer summary plot: dotplot of z-scored average gene module score (sized by % of spots with > 0 gene module score; labeled percent enriched)

# create data frame of mouse gene module scores, layer annotations
layer_modules_df <- as.data.frame(combo[[c("Sadick_mouse_gene_modules1", "Sadick_mouse_gene_modules2", "Sadick_mouse_gene_modules3", "Sadick_mouse_gene_modules4", 
                                           "Sadick_mouse_gene_modules5", "Sadick_mouse_gene_modules6", "Sadick_mouse_gene_modules7", "Sadick_mouse_gene_modules8", 
                                           "Sadick_mouse_gene_modules9", "final_layer_annotation")]])
colnames(layer_modules_df) <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8",
                                "Region")

# calculate average module scores per layer
layer_avg <- layer_modules_df %>% 
  group_by(Region) %>%
  summarise(across(everything(), mean))

layer_avg <- layer_avg[1:6,]

layer_avg[,2:10] <- scale(layer_avg[,2:10]) # z-score along each gene module
layer_avg <- reshape2::melt(layer_avg) # melt into long data frame
colnames(layer_avg) <- c("Region", "Module", "average_score") # rename columns

layer_modules_df <- reshape2::melt(layer_modules_df) # melt into long data frame
colnames(layer_modules_df) <- c("Region", "Module", "value") # rename columns 

# calculate "percent enriched"; % of spots in given layer with > 0 module score
layer_modules_df <- layer_modules_df %>%  
  group_by(Region, Module) %>% 
  summarise(count = n(),
            count_enriched = sum(value > 0),
            percent_enriched = count_enriched / count)

# merge two data frames
layer_module_stats <- merge(layer_avg, layer_modules_df, by = c("Region", "Module"))

# create data frame with spots and layer identities
human.layer.modules.df <- data.frame(layer = maynard.combo$layer)
# combine layers 2 and 3 to L2/3 to match mouse data
human.layer.modules.df <- human.layer.modules.df %>% mutate(layer = ifelse(layer == "Layer2" | layer == "Layer3", "L2/3", 
                                                                           ifelse(layer == "Layer1", "L1", 
                                                                                  ifelse(layer == "Layer4", "L4",
                                                                                         ifelse(layer == "Layer5", "L5",
                                                                                                ifelse(layer == "Layer6", "L6", layer))))))
# reorder column names
human.layer.modules.df[match(colnames(maynard.combo), rownames(human.layer.modules.df)),]

# assign updated layer labels to seurat object
maynard.combo$final_layer <- human.layer.modules.df$layer

# create data frame with module scores and layer identities
human_layer_modules_df <- as.data.frame(maynard.combo[[c("Sadick_human_gene_module1", "Sadick_human_gene_module2", "Sadick_human_gene_module3", "Sadick_human_gene_module4", 
                                                         "Sadick_human_gene_module5", "Sadick_human_gene_module6", "Sadick_human_gene_module7", "Sadick_human_gene_module8", 
                                                         "Sadick_human_gene_module9", "final_layer")]])
colnames(human_layer_modules_df) <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8",
                                      "Region")

# calculate average layer module scores
layer_avg <- human_layer_modules_df %>% 
  group_by(Region) %>%
  summarise(across(everything(), mean))

layer_avg <- layer_avg[1:6,]

# scale
layer_avg[,2:10] <- scale(layer_avg[,2:10])

# melt into long data frame
layer_avg <- reshape2::melt(layer_avg)
colnames(layer_avg) <- c("Region", "Module", "average_score")

# melt into long data frame
human_layer_modules_df <- reshape2::melt(human_layer_modules_df)
colnames(human_layer_modules_df) <- c("Region", "Module", "value")

# calculate "percentage enriched": the proportion of spots in each layer with > 0 module score
human_layer_modules_df <- human_layer_modules_df %>% 
  group_by(Region, Module) %>% 
  summarise(count = n(),
            count_enriched = sum(value > 0),
            percent_enriched = count_enriched / count)

# merge data frames
human_layer_module_stats <- merge(layer_avg, human_layer_modules_df, by = c("Region", "Module"))

# creating human region dendrogram; reorganize data
mat <- human_layer_module_stats %>% 
  select(-count, -count_enriched, -percent_enriched) %>%  
  pivot_wider(names_from = Region, values_from = average_score) %>% 
  data.frame() 
row.names(mat) <- mat$Module  
mat <- mat[,-1] 

# cluster
region_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
# create dendrogram
region_dend <- as.dendrogram(region_clust) 

# rename/rotate dendrogram leaves
labels(region_dend) <- gsub("\\.", "/", region_clust$labels[region_clust$order])
region_dend <- dendextend::rotate(region_dend, order = c("WM", "L1", "L2/3", "L4", "L5", "L6"))

# create ggplot dendrogrma
region_tree <- ggdendrogram(region_dend) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())


# creating mouse region dendrogram; reorganize data 
mat <- layer_module_stats %>% 
  select(-count, -count_enriched, -percent_enriched) %>% 
  pivot_wider(names_from = Region, values_from = average_score) %>% 
  data.frame() 
row.names(mat) <- mat$Module  
mat <- mat[,-1] 

# clustering
region_clust <- hclust(dist(mat %>% as.matrix() %>% t())) 
# create dendrogram
region_dend_mouse <- as.dendrogram(region_clust) 
# rename/reorder dendrogram leaves
labels(region_dend_mouse) <- gsub("\\.", "/", region_clust$labels[region_clust$order])
region_dend_mouse <- dendextend::rotate(region_dend_mouse, order = c("WM", "L1", "L2/3", "L4", "L5", "L6"))
# plot dendrogram
region_tree_mouse <- ggdendrogram(region_dend_mouse) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())


# reorganize data frames for clustering
mat1 <- human_layer_module_stats %>% 
  select(-count, -count_enriched, -percent_enriched) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = Region, values_from = average_score) %>% 
  data.frame() 
mat2 <- layer_module_stats %>% 
  select(-count, -count_enriched, -percent_enriched) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = Region, values_from = average_score) %>% 
  data.frame() 
row.names(mat1) <- mat1$Module  
mat1 <- mat1[,-1] 
row.names(mat2) <- mat2$Module  
mat2 <- mat2[,-1] 

# combine human and mouse data frames
mat <- cbind(mat1, mat2)

# cluster modules (across both human and mouse)
mod_clust <- hclust(dist(mat %>% as.matrix())) 
# create dendrogram
mod_dend <- as.dendrogram(mod_clust) 
# modify dendrogram
mod_tree <- ggdendrogram(mod_dend) + coord_flip() + scale_y_reverse() + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# create human dotplot 
dotplot <- human_layer_module_stats %>% mutate(Module = factor(Module, levels = mod_clust$labels[mod_clust$order]), Region = factor(Region, levels = labels(region_dend))) %>% 
  ggplot(aes(x=Region, y = Module, color = average_score, size = percent_enriched)) + 
  geom_point() + scale_size_area(labels = scales::percent, guide = guide_legend(
    direction = "horizontal",
    title.position = "top",
    order = 1,
    title.hjust = 0.5
  )) + 
  scale_color_viridis_c(name = 'Scaled Gene\nModule Score', option = "inferno", lim = range(mat),  guide = guide_colorbar(direction = "horizontal", title.position = "top", order = 2)) + 
  cowplot::theme_cowplot() + labs(x = "Cortical Region", y = NULL, size = "% of Spots Enriched") + 
  theme(axis.line  = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal")

# create human dotplot with module dendrogram
p1 <- plot_grid(mod_tree, NULL, dotplot + NoLegend(), nrow = 1, ncol = 3, rel_widths = c(0.3,-0.05, 2), align = 'h', scale = c(0.961, 1, 1))

# save human dotplot
postscript("human_cortical_region_dotpot.ps", width = 6, height = 4)
p1
dev.off()

# save human region dendrogram
postscript("human_cortical_region_dotpot_region_dendrogram.ps", width = 6, height = 5)
region_tree
dev.off()

# save legend
postscript("human_cortical_region_dotpot_legend.ps", width = 6, height = 5)
as_ggplot(get_legend(dotplot))
dev.off()

# create mouse dotplot
dotplot <- layer_module_stats %>% mutate(Module = factor(Module, levels = mod_clust$labels[mod_clust$order]), Region = factor(Region, levels = labels(region_dend_mouse))) %>% 
  ggplot(aes(x=Region, y = Module, color = average_score, size = percent_enriched)) + 
  geom_point() + scale_size_area(labels = scales::percent, guide = guide_legend(
    direction = "horizontal",
    title.position = "top",
    order = 1,
    title.hjust = 0.5
  )) + 
  scale_color_viridis_c(name = 'Scaled Gene\nModule Score', option = "inferno", lim = range(mat), guide = guide_colorbar(direction = "horizontal", title.position = "top", order = 2)) + 
  cowplot::theme_cowplot() + labs(x = "Cortical Region", y = NULL, size = "% of Spots Enriched") + #theme_pubr() + 
  theme(axis.line  = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal")

# plot dotplot + dendrogram
p1 <- plot_grid(mod_tree, NULL, dotplot + NoLegend(), nrow = 1, ncol = 3, rel_widths = c(0.3,-0.05, 2), align = 'h', scale = c(0.961, 1, 1))

# save mouse dotplot
postscript("mouse_cortical_region_dotpot.ps", width = 6, height = 4)
p1
dev.off()

# save mouse region dendrogram
postscript("mouse_cortical_region_dotpot_region_dendrogram.ps", width = 6, height = 5)
region_tree_mouse
dev.off()

# save legend
postscript("mouse_cortical_region_dotpot_legend.ps", width = 6, height = 5)
as_ggplot(get_legend(dotplot))
dev.off()

# merge layer module stats data frames for mouse and human together
comparison.df <- merge(layer_module_stats, human_layer_module_stats, by = c("Region", "Module")) # x = mouse, y = human

# plot correlation plot legend
postscript("cor_plot_legend.ps", width = 8, height = 4)
as_ggplot(get_legend(ggplot(comparison.df, aes(x = average_score.x, y = average_score.y, color = Module)) + 
                       geom_point(aes(shape = Region), size = 5) + 
                       scale_shape_manual(values=c(15, 16, 17, 18, 3, 8), guide = guide_legend(title.position = "top"), name = "Cortical Region") + 
                       scale_color_manual(values = c("#008ECE", "#59C7EB", "#E0607E", "#ECA0B2", "#077187", "#6AAAB7", "#8E2043", "#BC7A8F", "#0A9086"), guide = guide_legend(title.position = "top", ncol = 3), name = "Gene Module") + 
                       theme_pubr() + 
                       labs(x = "MOUSE\n(scaled gene module score)", y = "HUMAN\n(scaled gene module score)") + 
                       theme(aspect.ratio = 1)))
dev.off()

# plot correlation plot
postscript("layer_cor_plot.ps", width = 4, height = 4)
ggplot(comparison.df, aes(x = average_score.x, y = average_score.y, color = Module)) + 
  geom_point(aes(shape = Region), size = 3) + 
  scale_shape_manual(values=c(15, 16, 17, 18, 3, 8), guide = guide_legend(title.position = "top"), name = "Cortical Region") + 
  scale_color_manual(values = c("#008ECE", "#59C7EB", "#E0607E", "#ECA0B2", "#077187", "#6AAAB7", "#8E2043", "#BC7A8F", "#0A9086"), guide = guide_legend(title.position = "top", ncol = 3), name = "Gene Module") + 
  geom_smooth(method=lm, se=FALSE, color = "black") + 
  theme_pubr() + 
  labs(x = "MOUSE\n(scaled gene module score)", y = "HUMAN\n(scaled gene module score)") + 
  theme(aspect.ratio = 1) + NoLegend()
dev.off()

# calculate human-mouse correlation 
round(cor(comparison.df$average_score.x, comparison.df$average_score.y), digits = 2) # r = 0.75

# create a matrix of kruskal wallis test statistics for each cluster module, for human and mouse: 
mat <- as.matrix(data.frame("Human" = kruskal.test.results.human$statistic, "Mouse" = kruskal.test.results.mouse$statistic))
rownames(mat) <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8")

# scale
mat <- scale(mat)

# create heatmap
heat <- Heatmap(mat, name = "H", col = circlize::colorRamp2(breaks = c(-2, -1, 0, 1, 2), colors = viridis::viridis_pal(option = "viridis")(5)),
                show_column_dend = FALSE, cluster_columns = FALSE, row_names_side = "left") # very subtle but clear trend

# plot heatmap
postscript("kruskal_wallis_hstat_heatmap.ps", width = 4, height = 4)
heat
dev.off()

