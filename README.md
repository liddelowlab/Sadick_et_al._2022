# Sadick-et-al., XXX, XXX
Human astrocytes and oligodendrocytes to the max!

Libraries (including version) used to analyze single nuclei RNASeq (snRNA-seq) astrocyte data included in: Sadick et al., XXX, XXX

Included is the github libraries and code necessary to replicate the Seurat objects used for analysis as well as all differential gene expression analysis and GO/pathway analysis completed. The following versions of individual software/packages were used for analysis:

- R 3.6.2 - 4.0.3
- biomaRt 2.46.0
- clusterProfiler 3.18.0
- ComplexHeatmap 2.6.2
- cowplot 1.1.0
- dplyr 1.0.2
- edgeR 3.32.0
- EnhancedVolcano 1.8.0
- ggplot2 3.3.2
- gplots 3.1.0
- gridExtra 2.3
- grid 4.0.3
- inauguration 0.0.0.9000
- magrittr 2.0.1
- Matrix 1.2-18
- matrixStats 0.57.0
- org.Hs.eg.db 3.12.0
- pheatmap 1.0.12
- plyr 1.8.6
- purrr 0.3.4 
- ReactomePA 1.34.0
- RColorBrewer 1.1-2
- scater 1.18.3
- scds 1.6.0
- scran 1.18.1
- scRNAseq 2.4.0
- Seurat 3.2.2
- SingleCellExperiment 1.12.0
- topGO 2.42.0
- unikn 0.3.0
- UpSetR 1.4.0
- zinbwave 1.12.0

# Code included in repository:
No custom code included in analysis pipelines.

# Data
Both FASTQ and Cell Ranger-generated matrix files generated in this project can be directly downloaded from GEO: [GSE167494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167494)

**Reanalysis of Astrocytes from Published Mouse Datasets:**
Primary sources for reanalyzed data is summarized below.
| Dataset | Seq Used | Data Repository | Publication |
| --- | --- | --- | --- |
| Mathys | Chromium Single Cell 3' Gene Expression (v2) | [syn18485175](https://www.synapse.org/#!Synapse:syn18485175) | [Mathys et al. (2019) *Nature* 570:332-337](https://pubmed.ncbi.nlm.nih.gov/31042697/) |
| Grubman | Chromium Single Cell 3' Gene Expression (v2) | [GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852) | [Grubman et al. (2019) *Nature Neuroscience* 22:2087-2097](https://pubmed.ncbi.nlm.nih.gov/31768052/) |
| Zhou | Chromium Single Cell 5' Gene Expression | [syn21670836](https://www.synapse.org/#!Synapse:syn21670836) | [Zhou et al. (2020) *Nature Medicine* 26:131-142](https://pubmed.ncbi.nlm.nih.gov/31932797/) |

# Acknowledgements
Data analysis on astrocytes and oligodendrocytes was completed by Jessica Sadick with significant input on analysis pipeline by Taitea Dykstra.
