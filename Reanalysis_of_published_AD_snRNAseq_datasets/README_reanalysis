#Reanalysis of published AD snRNA-seq datasets

We obtained FASTQ files from previously published human AD snRNA-seq datasets: Mathys et al., 2019, Grubman et al., 2019, and Zhou et al., 2020. 
All snRNA-seq datasets were aligned to the premRNA-modified GRCh38 human reference genome using Cell Ranger software (version 4.0.0). 
Each dataset was processed for quality control, normalization, anchoring, integration, dimension reduction, clustering, annotation, and marker identification. 
For each dataset, astrocyte and oligodendrocyte nuclei were subsetted as unique Seurat objects for reanalysis in isolation. 
For Mathys astrocyte and oligodendrocyte subsetted analyses, multiple donors were removed because their nuclei yields were lower than the number of principle components used to evaluate the data. 
Please note that reference-based integration was used for all objects in Mathys and Zhou datasets due to memory constraints. 
Additionally, please note that two rounds of subsetting and analysis were required for Zhou astrocytes in order to remove contaminating, non-astrocytic nuclei/cells. 
