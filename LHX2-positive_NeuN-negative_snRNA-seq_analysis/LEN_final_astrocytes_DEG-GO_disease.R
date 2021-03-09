#Cluster-specific differential gene expression (DGE) and gene ontology(GO)/pathway analysis for astrocytes and oligodendrocytes in LIM Homeobox 2 (LHX2)-positive/NeuN-negative *FINAL* donor snRNA-seq dataset
#Goal: To evaluate DGE and GO terms specific to astrocyte and oligodendrocyte clusters comparing APOE e2/3 Alzheimer's disease versus age-matched non-symptomatic patients.
#Pipeline prepared by Taitea Dykstra and Jessica S. Sadick
#Majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#---------------------------------------------------------------------------------------------------
#Load libraries
library(SingleCellExperiment)
library(Seurat)
library(edgeR)
library(zinbwave)
library(scRNAseq)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(biomaRt)
library(scater)
library(scran)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(Matrix)
library(dplyr)
library(topGO)
library(EnhancedVolcano)

#Increase working memory
options(future.globals.maxSize = 2048 * 1024 ^20)

#Enable parallel processing for faster computations
BiocParallel::register(BiocParallel::SerialParam())
#Set working directory where files are located and will be saved
setwd('file_path')

#Read in subsetted and clustered astrocytes as SeuratObject
so <- readRDS('./CR4_LEN_final_so_astro_r2_20PC.rds')
#Define number of clusters
so <- SetIdent(so, value = "integrated_snn_res.0.3")

######ASTROCYTE CLUSTER 0######
#Subset for only astrocyte cluster 0
so.sub <- subset(so, idents=c(0))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
#Identify genes that have a count of at least 6 in at least 6 cells
filter <- rowSums(assay(sce)>5)>5
#Generate information on how many genes meet these criteria
table(filter)
#Filter matrix based on genes above
sce.filt <- sce[filter,]
#View matrix stats after filtering
sce.filt

##Select highly variable genes on which to focus analysis
#Model gene variance
sce.var <- modelGeneVar(sce.filt)
#Using gene variance, identify the top 2000
keep <- getTopHVGs(sce.var, n=2000)
#Filter matrix for only these genes
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
#Create zero-inflated negative binomial regression model to the data (note: this step is VERY computationally intensive)
#(ref: https://urldefense.proofpoint.com/v2/url?u=https-3A__bioconductor.org_packages_release_bioc_vignettes_zinbwave_inst_doc_intro.html-23normalized-2Dvalues-2Dand-2Ddeviance-2Dresiduals&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=gAum_CR2zJzIcfsIb-90NM5JpXJv8NNIreuA_1pR0rY&m=119tcA-SeS0G6xYL2dXy0nnVrliAEj3dmOWe9fiWKdk&s=TzAjyF_5R6nAQXh9Jb_KviDZPxTUXeztF2XRJyZHf_k&e= )
#Convert SingleCellExperiment object counts to matrix format
#(oberservational weights can only be computed from SummarizedExperiment or matrix)
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
#Isolate weights to pass to DGE calculations
weights <- assay(sce.zinb, "weights")
#Save model with weights
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster0.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster0.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
#Create design matrix specifying conditions and sex as variables
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)
#Create edgeR object from model with weights
dge <- DGEList(assay(sce.zinb))
#Caluclate normalization factors with TMM
dge <- calcNormFactors(dge)
#Add observational weights into the edgeR object
dge$weights <- weights
#Estimate dispersion, passing defined design matrix to function
dge <- estimateDisp(dge, design)
#Fit model
fit <- glmFit(dge, design)
#Make comparison between AD and NS (by specifying AD first, NS becomes control)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
#Zero-inflation adjusted F test for assessing DGE
#This object will have all of the DGE results
dge_res <- glmWeightedF(fit, contrast=contrast)
#Save DGE so it doesn't need to be recalculated
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster0.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster0.rds')

#Get counts of DGEs based on benjamini-hochberg adjusted p-value of less than 0.05
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
#View counts of DGEs
summary(is.de)

##Summarize results
#Basic MD plot to visualize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
#Filter DGEs meeting statistical significance mentioned above
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
#Create dataframe with gene name, log2 fold change, and adjusted p-value
genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
#Save dataframe as csv file
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster0-all.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
#Save dataframe as csv file
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster0-0.25lfc.csv')

##Pathway analysis + visulaization
#Isolate gene names/ENSEMBL IDs
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster0_0.25lfc.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster0_0.25lfc.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(7,60,68,71,72,73,79,82,83,98,111,115,117,134,135,144,150,162,164,166,168,
             173,176,178,188,192,196,199,213,222,225,236,286,296,308,310,314,340,342,
             345,348,350,359,368,372,374,375,376,378,379,385,389,404,407,409,411,413,
             416,429,431,433,434,437,439,441,442,444,445,446,447)] <- 
  c("ENSG00000280441","ENSG00000287134","ENSG00000236197","ENSG0000024708",
    "ENSG00000274441","ENSG00000286757","ENSG00000275830","ENSG00000132475",
    "ENSG00000233067","ENSG00000198804","ENSG00000259539","ENSG00000287159",
    "ENSG00000279168","ENSG00000271904","ENSG00000237978","ENSG00000233420",
    "ENSG00000214955","ENSG00000235904","ENSG00000274265","ENSG00000264630",
    "ENSG00000284624","ENSG00000286387","ENSG00000242593","ENSG00000285971",
    "ENSG00000232931","ENSG00000257545","ENSG00000237356","ENSG00000280683",
    "ENSG00000258603","ENSG00000263508","ENSG00000228113","ENSG00000255580",
    "ENSG00000260517","ENSG00000179066","ENSG00000286044","ENSG00000249335",
    "ENSG00000287001","ENSG00000234997","ENSG00000283380","ENSG00000287232",
    "ENSG00000287708","ENSG00000250195","ENSG00000239920","ENSG00000287044",
    "ENSG00000284196","ENSG00000237153","ENSG00000259345","ENSG00000286836",
    "ENSG00000285534","ENSG00000237742","ENSG00000254656","ENSG00000250403",
    "ENSG00000249330","ENSG00000224063","ENSG00000254532","ENSG00000226022",
    "ENSG00000242880","ENSG00000254279","ENSG00000287876","ENSG00000255910",
    "ENSG00000251129","ENSG00000257986","ENSG00000267284","ENSG00000274461",
    "ENSG00000259995","ENSG00000283692","ENSG00000249061","ENSG00000254394",
    "ENSG00000285751","ENSG00000259692")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

#Format as dataframe
ensembl_id <- as.data.frame(ensembl_id)
#Add log2 fold change and adjusted p-values
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

#Import biomaRt data for gene name conversions
mart <- useMart("ensembl","hsapiens_gene_ensembl")
#Select specific datasets from the import that will be of interest - gene name, ENSEMBL ID + entrezID
#(entrez ID is used for clusterprofiler pathway analysis)
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)
#Create second dataframe from dataframe above including biomaRt conversion from ENSEMBL ID to entrezID
annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)

#Separate dataframe into up and down DGEs for separate pathway analysis
#Shown to be more accurate than when all considered together:
#Hong G, Zhang W, Li H, Shen X, Guo Z. 2014 Separate enrichment analysis of pathways for up- and downregulated genes. J. R. Soc. Interface 11: 20130950. https://urldefense.proofpoint.com/v2/url?u=http-3A__dx.doi.org_10.1098_rsif.2013.0950&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=gAum_CR2zJzIcfsIb-90NM5JpXJv8NNIreuA_1pR0rY&m=119tcA-SeS0G6xYL2dXy0nnVrliAEj3dmOWe9fiWKdk&s=LXWAKMPB46NPo1c5zAVMw_R3X9XN7kHNJkaDtZuo600&e=
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

#Convert dataframe to matrices with just log2 fold change and entrezID
up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id
down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Run pathway analysis on upregulated genes
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster0_upregulated_pathways.csv')

#View pathway analysis as barplot, looking at top 10 pathways
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster0) Upregulated GO Biological Pathways",
        font.size = 16)

#View basic volcano plot
EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Run pathway analysis on downregulated genes
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster0_downregulated_pathways.csv')

#View pathway analysis as barplot, looking at top 10 pathways
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster0) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 1######
#Subset for only astrocyte cluster 1
so.sub <- subset(so, idents=c(1))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster1.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster1.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster1.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster1.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster1.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster1-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster1.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster1.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(6,43,50,78,118,119,124,145,148,150,167,190,200,202,223,227,241,255,272,296,
             304,305,338,352,356,365,376,377,380,383,386,394,400,401,403,405,410)] <- 
  c("ENSG00000280441","ENSG00000236197","ENSG00000287134","ENSG00000274441","ENSG00000283913",
    "ENSG00000233067","ENSG00000198804","ENSG00000247626","ENSG00000242593","ENSG00000237978",
    "ENSG00000259481","ENSG00000274265","ENSG00000250166","ENSG00000275830","ENSG00000233420",
    "ENSG00000271904","ENSG0000024708","ENSG00000228113","ENSG00000286033","ENSG00000132475",
    "ENSG00000284624","ENSG00000237356","ENSG00000285971","ENSG00000238755","ENSG00000232855",
    "ENSG00000249328","ENSG00000287876","ENSG00000285941","ENSG00000250195","ENSG00000283380",
    "ENSG00000254656","ENSG00000259345","ENSG00000251129","ENSG00000225249","ENSG00000249061",
    "ENSG00000231424","ENSG00000224063")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster1_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster1) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster1_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster1) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 2######
#Subset for only astrocyte cluster 2
so.sub <- subset(so, idents=c(2))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster2.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster2.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster2.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster2.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster2.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster2-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster2.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster2.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(1,7,9,15,31,38,59,72,88,97,98,112,117,131,132,137,167,172,189,191,194,226,227,
             231,232,238,241,248,251,254,264,265,269,270,272,278,280,298,299,300,301,304,305,
             306)] <- 
  c("ENSG00000125462","ENSG00000132475","ENSG00000280441","ENSG00000236197","ENSG00000275830",
    "ENSG00000231079","ENSG00000271860","ENSG00000271904","ENSG00000264630","ENSG00000287292",
    "ENSG00000285971","ENSG00000257545","ENSG00000250166","ENSG00000284977","ENSG00000287862",
    "ENSG00000232931","ENSG00000261404","ENSG00000237356","ENSG00000274441","ENSG00000241231",
    "ENSG00000225339","ENSG00000233420","ENSG0000024708","ENSG00000286863","ENSG00000259692",
    "ENSG00000162078","ENSG00000283380","ENSG00000236886","ENSG00000228061","ENSG00000279082",
    "ENSG00000224063","ENSG00000258342","ENSG00000285987","ENSG00000237742","ENSG00000274461",
    "ENSG00000249330","ENSG00000228113","ENSG00000249061","ENSG00000213121","ENSG00000259345",
    "ENSG00000287001","ENSG00000254532","ENSG00000253553","ENSG00000140105")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster2_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster2) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster2_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster2) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 3######
#Subset for only astrocyte cluster 3
so.sub <- subset(so, idents=c(3))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster3.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster3.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster3.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster3.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster3.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster3-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster3.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster3.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(30,32,69,98,99,119,124,136,140,152,157,194,212,215,231,238,239,241,243,281,312,
             328,334,344,348,367,371,378,400,408,427,428,431,433,435,443,490,492,504,507,526,
             539,541,554,557,574,578,580,587,594,602,609,619,625,627,630,633,637,643,646,648,
             656,658,660,662,667,668,669)] <- 
  c("ENSG00000286757","ENSG00000275830","ENSG00000233067","ENSG0000024708","ENSG00000280441",
    "ENSG00000198938","ENSG00000274265","ENSG00000198804","ENSG00000287277","ENSG00000279168",
    "ENSG00000251136","ENSG00000198712","ENSG00000287134","ENSG00000273118","ENSG00000237356",
    "ENSG00000132475","ENSG00000287159","ENSG00000198727","ENSG00000235904","ENSG00000233420",
    "ENSG00000203739","ENSG00000236197","ENSG00000237978","ENSG00000248079","ENSG00000283913",
    "ENSG00000232931","ENSG00000248932","ENSG00000225339","ENSG00000264630","ENSG00000253764",
    "ENSG00000251442","ENSG00000258603","ENSG00000274441","ENSG00000198899","ENSG00000235538",
    "ENSG00000284977","ENSG00000249335","ENSG00000285219","ENSG00000198886","ENSG00000284624",
    "ENSG00000254656","ENSG00000231079","ENSG00000242880","ENSG00000242593","ENSG00000242242",
    "ENSG00000285941","ENSG00000286533","ENSG00000258631","ENSG00000243176","ENSG00000224063",
    "ENSG00000231672","ENSG00000232855","ENSG00000285534","ENSG00000226022","ENSG00000238755",
    "ENSG00000286044","ENSG00000228065","ENSG00000232555","ENSG00000274461","ENSG00000259345",
    "ENSG00000250195","ENSG00000245768","ENSG00000258028","ENSG00000287001","ENSG00000251363",
    "ENSG00000227240","ENSG00000226868","ENSG00000249061")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster3_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster3) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster3_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster3) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 4######
#Subset for only astrocyte cluster 4
so.sub <- subset(so, idents=c(4))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster4.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster4.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster4.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster4.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster4.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster4-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- row.names(genes)
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster4.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster4.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(1,6,30,33,41,43,48,72,94,96,128,139,146,160,166,180,198,207,221,226,232,240,247,
             248,254,258,271,284,291,296,304,319,324,358,361,363,367,390,404,409,425,431,441,
             466,473,478,479,488,490,500,502,505,521,558,566,568,570,595,599)] <- 
  c("ENSG00000280441","ENSG00000125462","ENSG00000180573","ENSG00000236197","ENSG00000250166",
    "ENSG00000275830","ENSG00000236975","ENSG00000274265","ENSG00000284977","ENSG00000225339",
    "ENSG00000285987","ENSG00000286288","ENSG00000271860","ENSG00000140105","ENSG00000287704",
    "ENSG00000179300","ENSG00000235538","ENSG00000248994","ENSG00000231079","ENSG00000273247",
    "ENSG00000259255","ENSG00000251442","ENSG00000259336","ENSG00000249740","ENSG00000253944",
    "ENSG00000248079","ENSG00000251034","ENSG00000239381","ENSG00000228113","ENSG00000250326",
    "ENSG00000273118","ENSG00000230084","ENSG00000274461","ENSG00000287706","ENSG00000287134",
    "ENSG00000286533","ENSG00000260971","ENSG00000249335","ENSG00000226994","ENSG00000229956",
    "ENSG00000237877","ENSG00000237978","ENSG00000286757","ENSG00000288075","ENSG00000236886",
    "ENSG00000231252","ENSG00000227486","ENSG00000259345","ENSG00000253320","ENSG00000251314",
    "ENSG00000285367","ENSG00000231918","ENSG00000279168","ENSG00000132475","ENSG00000261404",
    "ENSG00000287292","ENSG00000233797","ENSG00000257060","ENSG00000198727")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes$table$logFC
ensembl_id$Pval <- genes$table$FDR

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster4_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster4) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster4_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster4) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 5######
#Subset for only astrocyte cluster 5
so.sub <- subset(so, idents=c(5))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster5.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster5.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster5.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster5.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster5.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster5-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster5.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster5.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(3,4,6,40,45,65,80,86,96,107,116,122,127,160,161,173,181,186,187,208,214,222,224,
             236,270,274,283,286,291,305,312,319,320,334,337,339,344,378,388,390,392,402,403,
             414,416,426,435,436,446,452,454,460,462,463,473,477,478)] <- 
  c("ENSG00000125462","ENSG00000132475","ENSG00000280441","ENSG00000284977","ENSG00000285987",
    "ENSG00000287372","ENSG00000230084","ENSG00000275830","ENSG00000258342","ENSG00000274265",
    "ENSG00000271904","ENSG00000231079","ENSG00000237742","ENSG00000249740","ENSG00000285971",
    "ENSG00000236197","ENSG00000271860","ENSG00000251034","ENSG00000286533","ENSG00000198886",
    "ENSG00000283380","ENSG00000251442","ENSG00000198840","ENSG00000273118","ENSG00000239381",
    "ENSG00000273247","ENSG00000287292","ENSG00000198899","ENSG00000253764","ENSG00000259692",
    "ENSG00000236975","ENSG00000198786","ENSG00000228061","ENSG00000231246","ENSG00000286162",
    "ENSG00000288088","ENSG00000166444","ENSG00000225249","ENSG00000226994","ENSG00000212907",
    "ENSG00000250166","ENSG00000248138","ENSG00000198727","ENSG00000237356","ENSG00000224063",
    "ENSG00000280241","ENSG00000225339","ENSG00000284959","ENSG00000239920","ENSG00000287544",
    "ENSG00000260971","ENSG00000257322","ENSG00000236886","ENSG00000280870","ENSG00000242880",
    "ENSG00000286072","ENSG00000254532")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster5_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster5) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster5_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster5) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 6######
#Subset for only astrocyte cluster 6
so.sub <- subset(so, idents=c(6))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster6.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster6.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster6.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster6.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster6.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster6-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster6.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster6.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(1,6,14,34,74,78,84,86,92,118,145,162,170,182,199,203,210,213,224,228,229,242,247,
             248,266,276,288,294,299,300,316,317,319,320,322,327,338,340)] <- 
  c("ENSG00000280441","ENSG00000125462","ENSG00000236197","ENSG00000274441","ENSG00000250166",
    "ENSG00000180573","ENSG00000271904","ENSG00000287134","ENSG00000274265","ENSG00000271860",
    "ENSG00000242593","ENSG00000255052","ENSG00000286288","ENSG00000198840","ENSG00000279082",
    "ENSG00000239381","ENSG00000236975","ENSG00000228113","ENSG00000248079","ENSG00000287292",
    "ENSG00000250326","ENSG00000273118","ENSG00000198888","ENSG00000284977","ENSG00000162078",
    "ENSG00000285971","ENSG00000251513","ENSG00000264630","ENSG00000187837","ENSG00000248138",
    "ENSG00000229956","ENSG00000233797","ENSG00000283380","ENSG00000236886","ENSG00000285987",
    "ENSG00000251442","ENSG00000286533","ENSG00000272865")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster6_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster6) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster6_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster6) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 7######
#Subset for only astrocyte cluster 7
so.sub <- subset(so, idents=c(7))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster7.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster7.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster7.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster7.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster7.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster7-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster7.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster7.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(1,6,15,22,26,28,43,46,49,62,64,91,108,114,122,125)] <- 
  c("ENSG00000280441","ENSG00000125462","ENSG00000236197","ENSG00000237356","ENSG00000231079",
    "ENSG00000233067","ENSG00000198804","ENSG00000274441","ENSG00000198763","ENSG00000198888",
    "ENSG00000285971","ENSG00000198840","ENSG00000249740","ENSG00000279082","ENSG00000287134",
    "ENSG00000273118")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save df as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluter7_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster7) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster7_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster7) Downregulated GO Biological Pathways",
        font.size = 16)

######ASTROCYTE CLUSTER 8######
#Subset for only astrocyte cluster 8
so.sub <- subset(so, idents=c(8))
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so.sub, assay="RNA")

##Filter out lowly expressed genes
filter <- rowSums(assay(sce)>5)>5
table(filter)
sce.filt <- sce[filter,]
sce.filt

##Select highly variable genes on which to focus analysis
sce.var <- modelGeneVar(sce.filt)
keep <- getTopHVGs(sce.var, n=2000)
sce.filt <- sce.filt[keep,]

##Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
weights <- assay(sce.zinb, "weights")
saveRDS(sce.zinb, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster8.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_weights_cluster8.rds')
#Isolate weights to pass to DGE calculations--load into local environment when model ran on HPC
weights <- assay(sce.zinb, "weights")

##DGE analysis with edgeR
condition <- factor(sce.zinb$disease_id)
sex <- factor(sce.zinb$sex_id)
design <- model.matrix(~0+condition+sex)

dge <- DGEList(assay(sce.zinb))
dge <- calcNormFactors(dge)
dge$weights <- weights
dge <- estimateDisp(dge, design)

fit <- glmFit(dge, design)
contrast <- makeContrasts(conditionAD-conditionNS, levels=design)
dge_res <- glmWeightedF(fit, contrast=contrast)
saveRDS(dge_res, file='CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster8.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_astro_r2_20PC_zinb_DGE_result_cluster8.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster8.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_sig_genes_cluster8-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Astro_gene_id_cluster8.csv')
write.csv(ensembl_id, file='./Astro_ensembl_id_cluster8.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(2,9,13,27,40,49,62,74,85,90,93,106,113)] <- 
  c("ENSG00000280441","ENSG00000125462","ENSG00000274265","ENSG00000236197","ENSG00000198804",
    "ENSG00000250166","ENSG00000237978","ENSG00000204929","ENSG00000287134","ENSG0000024708",
    "ENSG00000271904","ENSG00000274441","ENSG00000242593")
#Verify updated ENSEMBL IDs for NAs
sum(is.na(ensembl_id))

ensembl_id <- as.data.frame(ensembl_id)
ensembl_id$logFC <- genes.sig.lfc$logFC
ensembl_id$Pval <- genes.sig.lfc$Pval

mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ensembl_id$ensembl_id,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  ensembl_id[match(annotLookup$ensembl_gene_id, ensembl_id$ensembl_id),],
  annotLookup)
annotLookup_up <- subset(annotLookup, annotLookup$logFC > 0)
annotLookup_down <- subset(annotLookup, annotLookup$logFC < 0)

up_matrix <- annotLookup_up$logFC
names(up_matrix) <- annotLookup_up$entrezgene_id

down_matrix <- annotLookup_down$logFC
names(down_matrix) <- annotLookup_down$entrezgene_id

#Upregulated assessment
go_enrich <- enrichGO(gene = names(up_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
upregulated <- as.data.frame(go_enrich)
write.csv(upregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster8_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster8) Upregulated GO Biological Pathways",
        font.size = 16)

EnhancedVolcano(genes.sig, lab = row.names(genes.sig), x = "logFC", y = "Pval",
                pCutoff = 10e-6, FCcutoff = 0)

#Downregulated assessment
go_enrich <- enrichGO(gene = names(down_matrix),
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
downregulated <- as.data.frame(go_enrich)
write.csv(downregulated, file='./CR4_LEN_final_so_astro_r2_20PC_zinb_cluster8_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Astrocytes (cluster8) Downregulated GO Biological Pathways",
        font.size = 16)
