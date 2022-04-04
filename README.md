# [Sadick, O'Dea et al. 2022, Neuron](https://doi.org/10.1016/j.neuron.2022.03.008)
Human astrocytes and oligodendrocytes to the max!

This repository includes the code used in the analysis of single nuclei RNASeq (snRNA-seq) astrocyte and oligodendrocyte data from [Sadick, O'Dea et al. 2022](https://doi.org/10.1016/j.neuron.2022.03.008).

All analysis was performed using R 3.6-4.1. See 'packages.txt' for versions of individual software/packages were used for analysis.

# Data
Both FASTQ and Cell Ranger-generated matrix files generated in this project can be directly downloaded from GEO: [GSE167494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167494)

**Reanalysis of Astrocytes from Published Mouse Datasets:**
Primary sources for reanalyzed data are summarized below.
| Dataset | Seq Used | Data Repository | Publication |
| --- | --- | --- | --- |
| Mathys | Chromium Single Cell 3' Gene Expression (v2) | [syn18485175](https://www.synapse.org/#!Synapse:syn18485175) | [Mathys et al. (2019) *Nature* 570:332-337](https://pubmed.ncbi.nlm.nih.gov/31042697/) |
| Grubman | Chromium Single Cell 3' Gene Expression (v2) | [GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852) | [Grubman et al. (2019) *Nature Neuroscience* 22:2087-2097](https://pubmed.ncbi.nlm.nih.gov/31768052/) |
| Zhou | Chromium Single Cell 5' Gene Expression | [syn21670836](https://www.synapse.org/#!Synapse:syn21670836) | [Zhou et al. (2020) *Nature Medicine* 26:131-142](https://pubmed.ncbi.nlm.nih.gov/31932797/) |

**Analysis of Published Spatial Transcriptomics Datasets:**
Primary sources for this data are summarized below.
| Dataset | Data Repository | Publication |
| --- | --- | --- |
| Maynard | [Globus: jhpce#HumanPilot10x](http://research.libd.org/globus/jhpce_HumanPilot10x/index.html) | [Maynard et al. (2021) *Nature Neuroscience* 24:425–436](https://doi.org/10.1038/s41593-020-00787-0) |
| Hasel | [GSE165098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165098) | [Hasel et al. (2021) *Nature Neuroscience* 24:1475–1487](https://doi.org/10.1038/s41593-021-00905-6) |

# Acknowledgements
Data analysis was completed by Jessica Sadick, Michael O'Dea, and Taitea Dykstra.
