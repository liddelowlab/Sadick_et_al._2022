#Analysis of published spatial transcriptomics datasets

In order to examine the spatial localization of astrocyte subtype gene signatures, we examined previously published human and mouse spatial transcriptomics datasets.

'Hasel_visium_preprocessing.R' and 'Maynard_visium_preprocessing.R' provide code for the generation of Seurat objects from CellRanger output files. 'Creating_gene_modules.R' provides code for the generation of human and mouse astrocyte cluster gene modules. 'Regional_annotation.R' provides code for the regional annotation of the mouse spatial transcriptomics data. 'Differential_Enrichment_Analysis.R' provides code for the statistical testing of differential enrichment of cluster gene signatures across cortical regions and conditions. 'Visualization.R' provides code for the visualization of astrocyte cluster gene signature enrichment in the human and mouse spatial transcriptomics data.   

This analysis was performed using R 4.1.1. For the versions of R packages used (and dependencies), see '\[script_name]\_packages.txt' files. 
