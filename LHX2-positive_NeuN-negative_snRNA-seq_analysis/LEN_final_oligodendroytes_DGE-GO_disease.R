#Cluster-specific differential gene expression (DGE) and gene ontology(GO)/pathway analysis for oligodendrocytes in LIM Homeobox 2 (LHX2)-positive/NeuN-negative *FINAL* donor snRNA-seq dataset
#Goal: To evaluate DGE and GO terms specific to oligodendrocyte clusters comparing APOE e2/3 Alzheimer's disease versus age-matched non-symptomatic patients.
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

#Read in subsetted and clustered oligodendrocytes as SeuratObject
so <- readRDS('./CR4_LEN_final_so_oligo_15PC.rds')
#Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(so, assay="RNA")

######OLIGODENDROCYTE CLUSTER 0######
#Subset for only oligodendrocyte cluster 0
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
saveRDS(sce.zinb, file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster0.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster0.rds')
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
saveRDS(dge_res, file='CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster0.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster0.rds')
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
write.csv(genes.sig, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster0.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
#Save dataframe as csv file
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster0-0.25lfc.csv')

##Pathway analysis + visulaization
#Isolate gene names/ENSEMBL IDs
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Oligo_gene_id_cluster0_0.25lfc.csv')
write.csv(ensembl_id, file='./Oligo_ensembl_id_cluster0_0.25lfc.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(17,19,23,52,59,60,76,77,81,90,93,99,103,104,148,152,155,168,175,180,184,188,189,192,193,195,203,208,210,225,231,233,240,248,252,255,257,259,274,276,282,286,287,294,297,303,305,307,308,310,311,312,313,315,318,320,328,331,332,333,335,337,340,342)] <-
c("ENSG00000227088","ENSG00000280441","ENSG00000260788","ENSG00000198840","ENSG00000260664","ENSG00000207955","ENSG00000280683","ENSG00000237949","ENSG00000287306","ENSG00000198727","ENSG00000273118","ENSG00000229989","ENSG00000253807","ENSG00000198899","ENSG00000230084","ENSG00000224536","ENSG00000287299","ENSG00000267922","ENSG00000250195","ENSG00000250220","ENSG00000287292","ENSG00000227757","ENSG00000253944","ENSG00000250326","ENSG00000260392","ENSG00000285667","ENSG00000286044","ENSG00000286856","ENSG00000285783","ENSG00000284977","ENSG00000259336","ENSG00000237813","ENSG00000248752","ENSG00000240687","ENSG00000246225","ENSG00000258561","ENSG00000241231","ENSG00000259124","ENSG00000283003","ENSG00000260517","ENSG00000274461","ENSG00000286924","ENSG00000242242","ENSG00000287135","ENSG00000287256","ENSG00000282943","ENSG00000287862","ENSG00000260042","ENSG00000281641","ENSG00000257261","ENSG00000285095","ENSG00000285769","ENSG00000287801","ENSG00000247081","ENSG00000286321","ENSG00000225746","ENSG00000287683","ENSG00000287001","ENSG00000286438","ENSG00000229491","ENSG00000243276","ENSG00000286797","ENSG00000287097","ENSG00000224063")
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
write.csv(upregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster0_upregulated_pathways.csv')

#View pathway analysis as barplot, looking at top 10 pathways
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster0) Upregulated GO Biological Pathways",
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
write.csv(downregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster0_downregulated_pathways.csv')

#View pathway analysis as barplot, looking at top 10 pathways
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster0) Downregulated GO Biological Pathways",
        font.size = 16)

######OLIGODENDROCYTE CLUSTER 1######
#Subset for only oligodendrocyte cluster 1
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
saveRDS(sce.zinb, file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster1.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster1.rds')
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
saveRDS(dge_res, file='CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster1.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster1.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster1.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster1-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Oligo_gene_id_cluster1_0.25lfc.csv')
write.csv(ensembl_id, file='./Oligo_ensembl_id_cluster1_0.25lfc.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(1,5,11,14,15,36,39,63,66,73,82,86,89,96,103,107,108,117,129,142,168,169,173,183,189,191,195,212,213,219,225,226,228,253,255,260)] <-
  c("ENSG00000280441","ENSG00000198840","ENSG00000260788","ENSG00000198899","ENSG00000198727",
    "ENSG00000198888","ENSG00000249042","ENSG00000198763","ENSG00000255311","ENSG00000198712",
    "ENSG00000250159","ENSG00000198886","ENSG00000207955","ENSG00000198938","ENSG00000231079",
    "ENSG00000285783","ENSG00000227088","ENSG00000260664","ENSG00000262223","ENSG00000233559",
    "ENSG00000253154","ENSG00000230084","ENSG00000260971","ENSG00000253807","ENSG00000232855",
    "ENSG00000273118","ENSG00000237166","ENSG00000250220","ENSG00000250195","ENSG00000241231",
    "ENSG00000286162","ENSG00000248752","ENSG00000286856","ENSG00000287277","ENSG00000287862",
    "ENSG00000287092")
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
write.csv(upregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster1_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster1) Upregulated GO Biological Pathways",
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
write.csv(downregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster1_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster1) Downregulated GO Biological Pathways",
        font.size = 16)

######OLIGODENDROCYTE CLUSTER 2######
#Subset for only oligodendrocyte cluster 2
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
saveRDS(sce.zinb, file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster2.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster2.rds')
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
saveRDS(dge_res, file='CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster2.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster2.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster2.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster2-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Oligo_gene_id_cluster2_0.25lfc.csv')
write.csv(ensembl_id, file='./Oligo_ensembl_id_cluster2_0.25lfc.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(15,21,32,85,90,108,112,117,131,132,134,141,158,193,196,197,199,208,217)] <- 
  c("ENSG00000260788","ENSG00000227088","ENSG00000280441","ENSG00000198840","ENSG00000286438",
    "ENSG00000259678","ENSG00000250326","ENSG00000261329","ENSG00000286162","ENSG00000260664",
    "ENSG00000250195","ENSG00000287306","ENSG00000233559","ENSG00000248752","ENSG00000198727",
    "ENSG00000207955","ENSG00000198899","ENSG00000286279","ENSG00000229989")
#Verify updated ensembls for NAs
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
write.csv(upregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster2_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster2) Upregulated GO Biological Pathways",
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
write.csv(downregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster2_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster2) Downregulated GO Biological Pathways",
        font.size = 16)

######OLIGODENDROCYTE CLUSTER 3######
#Subset for only oligodendrocyte cluster 3
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
saveRDS(sce.zinb, file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster3.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster3.rds')
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
saveRDS(dge_res, file='CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster3.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster3.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster3.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster3-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Oligo_gene_id_cluster3_0.25lfc.csv')
write.csv(ensembl_id, file='./Oligo_ensembl_id_cluster3_0.25lfc.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(19,21)] <- 
  c("ENSG00000280441","ENSG00000286749")
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
write.csv(upregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster3_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster3) Upregulated GO Biological Pathways",
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
write.csv(downregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster3_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster3) Downregulated GO Biological Pathways",
        font.size = 16)

######OLIGODENDROCYTE CLUSTER 4######
#Subset for only oligodendrocyte cluster 4
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
saveRDS(sce.zinb, file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster4.rds')

#Load model with weights if not already in environment
sce.zinb <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_weights_cluster4.rds')
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
saveRDS(dge_res, file='CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster4.rds')

#Load DGE results if not already in environment
dge_res <- readRDS(file='./CR4_LEN_final_so_oligo_15PC_zinb_DGE_result_cluster4.rds')
is.de <- decideTestsDGE(dge_res, adjust.method = "BH", p.value = 0.05, lfc = 0.25)
summary(is.de)

##Summarize results
plotMD(dge_res, status=is.de, values=c(1,-1), col=c("red","blue"), main = "EdgeR")
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
write.csv(genes.sig, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster4.csv')

#Apply threshold for 0.25 log2 fold change cut off
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
write.csv(genes.sig.lfc, file='./CR4_LEN_final_so_oligo_15PC_zinb_sig_genes_cluster4-0.25lfc.csv')

##Pathway analysis + visulaization
gene_id <- genes.sig.lfc[,1]
#Convert gene names to ENSEMBL IDs
ensembl_id <- mapIds(org.Hs.eg.db, keys = gene_id, keytype = "SYMBOL", column= ("ENSEMBL"))
#Save gene_id and ensembl_id files
write.csv(gene_id, file='./Oligo_gene_id_cluster4_0.25lfc.csv')
write.csv(ensembl_id, file='./Oligo_ensembl_id_cluster4_0.25lfc.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
ensembl_id[c(6,10)] <- 
  c("ENSG00000198840","ENSG00000198888")
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
write.csv(upregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster4_upregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster4) Upregulated GO Biological Pathways",
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
write.csv(downregulated, file='./CR4_LEN_final_so_oligo_15PC_zinb_cluster4_downregulated_pathways.csv')

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "AD Oligodendrocytes (cluster4) Downregulated GO Biological Pathways",
        font.size = 16)
