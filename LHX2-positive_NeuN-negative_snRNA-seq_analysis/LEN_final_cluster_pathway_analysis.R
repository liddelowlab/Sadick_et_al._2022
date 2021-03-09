#Cluster-specific pathway analysis for astrocytes and oligodendrocytes in LIM Homeobox 2 (LHX2)-positive/NeuN-negative *FINAL* donor snRNA-seq dataset
#Goal: To evaluate gene ontology (GO)/pathway analysis terms specific to astrocyte and oligodendrocyte clusters across all donors in APOE e2/3 Alzheimer's disease and age-matched non-symptomatic patient final cohort. NOTE: *NOT* evaluating differences between disease states.
#Pipeline prepared by Jessica S. Sadick
#Majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#---------------------------------------------------------------------------------------------------
#Load libraries
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

###ASTROCYTE CLUSTERS PATHWAY ANALYSES###
#Read csv file for CR4_LEN_final_so_astro_r2_20PC cluster-specific/enriched transcripts
astro <- read.csv(file = "./LEN_final_so_astro_r2_20PC_res0.3_genes-RNA.csv")

#Subset by astrocyte cluster
astro_C0 <- subset(astro, select = c(gene), astro$cluster == "0")
astro_C1 <- subset(astro, select = c(gene), astro$cluster == "1")
astro_C2 <- subset(astro, select = c(gene), astro$cluster == "2")
astro_C3 <- subset(astro, select = c(gene), astro$cluster == "3")
astro_C4 <- subset(astro, select = c(gene), astro$cluster == "4")
astro_C5 <- subset(astro, select = c(gene), astro$cluster == "5")
astro_C6 <- subset(astro, select = c(gene), astro$cluster == "6")
astro_C7 <- subset(astro, select = c(gene), astro$cluster == "7")
astro_C8 <- subset(astro, select = c(gene), astro$cluster == "8")

##Astrocyte cluster 0##
#Convert gene names to ENSEMBL IDs
astro_C0_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C0$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
#Verify updated ENSEMBL IDs for NAs
sum(is.na(astro_C0_ensembl))
#Save gene_id and ensembl_id files
write.csv(astro_C0, file='./Astro_gene_id_C0.csv')
write.csv(astro_C0_ensembl, file='./Astro_ensembl_id_C0.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
astro_C0_ensembl[c(9,25,31,33,42,48)] <- 
  c("ENSG00000287544","ENSG00000271860","ENSG00000287134","ENSG00000204929","ENSG00000279082",
    "ENSG00000233067")
#Convert to dataframe
astro_C0_ensembl <- as.data.frame(astro_C0_ensembl)
#Convert ENSEMBL IDs to entrezgene names
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C0_ensembl$astro_C0_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C0_ensembl[match(annotLookup$ensembl_gene_id, astro_C0_ensembl$astro_C0_ensembl),],
  annotLookup)
#Run GO/pathway analysis
go_enrich_C0 <- enrichGO(gene = annotLookup$entrezgene_id,
                      OrgDb = 'org.Hs.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
#Save dataframe as csv file
C0_pathways <- as.data.frame(go_enrich_C0)
write.csv(C0_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C0_pathways.csv')

##Astrocyte cluster 1##
astro_C1_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C1$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C1_ensembl))
write.csv(astro_C1, file='./Astro_gene_id_C1.csv')
write.csv(astro_C1_ensembl, file='./Astro_ensembl_id_C1.csv')
astro_C1_ensembl[c(18,22,33,35,46,50,55,58,59,78,82,85,86,89,95,109,112,114,119,122,128,130,133)] <-
  c("ENSG00000259481","ENSG00000255580","ENSG00000198840","ENSG00000198938","ENSG00000287159",
    "ENSG00000198727","ENSG00000233559","ENSG00000274441","ENSG00000198712","ENSG00000198899",
    "ENSG00000257545","ENSG00000223653","ENSG00000198888","ENSG00000264630","ENSG00000249373",
    "ENSG00000237596","ENSG00000224731","ENSG00000232885","ENSG00000230910","ENSG00000198886",
    "ENSG00000240291","ENSG00000235904","ENSG00000198804")
astro_C1_ensembl <- as.data.frame(astro_C1_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C1_ensembl$astro_C1_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C1_ensembl[match(annotLookup$ensembl_gene_id, astro_C1_ensembl$astro_C1_ensembl),],
  annotLookup)
go_enrich_C1 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C1_pathways <- as.data.frame(go_enrich_C1)
write.csv(C1_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C1_pathways.csv')

##Astrocyte cluster 2##
astro_C2_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C2$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C2_ensembl))
write.csv(astro_C2, file='./Astro_gene_id_C2.csv')
write.csv(astro_C2_ensembl, file='./Astro_ensembl_id_C2.csv')
astro_C2_ensembl[c(2,5,20,22,30,33,47,50,51,54)] <-
  c("ENSG00000287704","ENSG00000259255","ENSG00000228065","ENSG00000280441","ENSG00000279168",
    "ENSG00000287277","ENSG00000179300","ENSG00000232931","ENSG00000126870","ENSG00000259972")
astro_C2_ensembl <- as.data.frame(astro_C2_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C2_ensembl$astro_C2_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C2_ensembl[match(annotLookup$ensembl_gene_id, astro_C2_ensembl$astro_C2_ensembl),],
  annotLookup)
go_enrich_C2 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C2_pathways <- as.data.frame(go_enrich_C2)
write.csv(C2_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C2_pathways.csv')

##Astrocyte cluster 3##
astro_C3_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C3$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C3_ensembl))
write.csv(astro_C3, file='./Astro_gene_id_C3.csv')
write.csv(astro_C3_ensembl, file='./Astro_ensembl_id_C3.csv')
astro_C3_ensembl[c(24,25,41,49,60,78,87,95,97)] <-
  c("ENSG00000225339","ENSG00000275830","ENSG00000249740","ENSG00000140105","ENSG00000247081",
    "ENSG00000251034","ENSG00000286533","ENSG00000179300","ENSG00000259616")
astro_C3_ensembl <- as.data.frame(astro_C3_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C3_ensembl$astro_C3_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C3_ensembl[match(annotLookup$ensembl_gene_id, astro_C3_ensembl$astro_C3_ensembl),],
  annotLookup)
go_enrich_C3 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C3_pathways <- as.data.frame(go_enrich_C3)
write.csv(C3_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C3_pathways.csv')

##Astrocyte cluster 4##
astro_C4_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C4$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C4_ensembl))
write.csv(astro_C4, file='./Astro_gene_id_C4.csv')
write.csv(astro_C4_ensembl, file='./Astro_ensembl_id_C4.csv')
astro_C4_ensembl[c(16,41,45,77,81,84,109,115,119,136,141,152,157,166,167,174,239,291,334,338,349,350,364,374,393)] <-
  c("ENSG00000251442","ENSG00000237356","ENSG00000236975","ENSG00000287372","ENSG00000228408",
    "ENSG00000287939","ENSG00000286776","ENSG00000287704","ENSG00000251136","ENSG00000228065",
    "ENSG00000258342","ENSG00000258844","ENSG00000231672","ENSG00000287706","ENSG00000275448",
    "ENSG00000257319","ENSG00000259255","ENSG00000126870","ENSG00000285987","ENSG00000230084",
    "ENSG00000280441","ENSG00000231615","ENSG00000288075","ENSG00000180573","ENSG00000284977")
astro_C4_ensembl <- as.data.frame(astro_C4_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C4_ensembl$astro_C4_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C4_ensembl[match(annotLookup$ensembl_gene_id, astro_C4_ensembl$astro_C4_ensembl),],
  annotLookup)
go_enrich_C4 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C4_pathways <- as.data.frame(go_enrich_C4)
write.csv(C4_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C4_pathways.csv')

##Astrocyte cluster 5##
astro_C5_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C5$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C5_ensembl))
write.csv(astro_C5, file='./Astro_gene_id_C5.csv')
write.csv(astro_C5_ensembl, file='./Astro_ensembl_id_C5.csv')
astro_C5_ensembl[c(5,11,25,46,53,54,68,71,85,95,101,112,127,146,148,156,171,197,201,206,217,271,278,286,316,329,343,350,351,354,361,368,376,383)] <-
  c("ENSG00000287704","ENSG00000259255","ENSG00000179300","ENSG00000287372","ENSG00000258342",
    "ENSG00000258844","ENSG00000231672","ENSG00000228065","ENSG00000274461","ENSG00000287893",
    "ENSG00000287939","ENSG00000225249","ENSG00000275830","ENSG00000288088","ENSG00000285987",
    "ENSG00000251442","ENSG00000285744","ENSG00000286695","ENSG00000283380","ENSG00000251136",
    "ENSG00000239920","ENSG00000230084","ENSG00000232931","ENSG00000287277","ENSG00000166444",
    "ENSG00000175170","ENSG00000257261","ENSG00000257060","ENSG00000248138","ENSG00000287706",
    "ENSG00000270207","ENSG00000261404","ENSG00000284977","ENSG00000132475")
astro_C5_ensembl <- as.data.frame(astro_C5_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C5_ensembl$astro_C5_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C5_ensembl[match(annotLookup$ensembl_gene_id, astro_C5_ensembl$astro_C5_ensembl),],
  annotLookup)
go_enrich_C5 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C5_pathways <- as.data.frame(go_enrich_C5)
write.csv(C5_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C5_pathways.csv')

##Astrocyte cluster 6##
astro_C6_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C6$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C6_ensembl))
write.csv(astro_C6, file='./Astro_gene_id_C6.csv')
write.csv(astro_C6_ensembl, file='./Astro_ensembl_id_C6.csv')
astro_C6_ensembl[c(6,34,80,91)] <-
  c("ENSG00000286776","ENSG00000237356","ENSG00000132475","ENSG00000241231")
astro_C6_ensembl <- as.data.frame(astro_C6_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C6_ensembl$astro_C6_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C6_ensembl[match(annotLookup$ensembl_gene_id, astro_C6_ensembl$astro_C6_ensembl),],
  annotLookup)
go_enrich_C6 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C6_pathways <- as.data.frame(go_enrich_C6)
write.csv(C6_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C6_pathways.csv')

##Astrocyte cluster 7##
astro_C7_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C7$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C7_ensembl))
write.csv(astro_C7, file='./Astro_gene_id_C7.csv')
write.csv(astro_C7_ensembl, file='./Astro_ensembl_id_C7.csv')
astro_C7_ensembl[c(3,10,13)] <-
  c("ENSG00000286811","ENSG00000237978","ENSG00000224592")
astro_C7_ensembl <- as.data.frame(astro_C7_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C7_ensembl$astro_C7_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C7_ensembl[match(annotLookup$ensembl_gene_id, astro_C7_ensembl$astro_C7_ensembl),],
  annotLookup)
go_enrich_C7 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C7_pathways <- as.data.frame(go_enrich_C7)
write.csv(C7_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C7_pathways.csv')

##Astrocyte cluster 8##
astro_C8_ensembl <- mapIds(org.Hs.eg.db, keys = astro_C8$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(astro_C8_ensembl))
write.csv(astro_C8, file='./Astro_gene_id_C8.csv')
write.csv(astro_C8_ensembl, file='./Astro_ensembl_id_C8.csv')
astro_C8_ensembl[c(6,27,52,53)] <-
  c("ENSG00000271860","ENSG00000228222","ENSG00000235450","ENSG00000271904")
astro_C8_ensembl <- as.data.frame(astro_C8_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=astro_C8_ensembl$astro_C8_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  astro_C8_ensembl[match(annotLookup$ensembl_gene_id, astro_C8_ensembl$astro_C8_ensembl),],
  annotLookup)
go_enrich_C8 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C8_pathways <- as.data.frame(go_enrich_C8)
write.csv(C8_pathways, file='./CR4_LEN_final_so_astro_r2_20PC_C8_pathways.csv')

#---------------------------------------------------------------------------------------------------
#Load libraries
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

###OLIGODENDROCYTE CLUSTERS PATHWAY ANALYSES###
#Read csv file for CR4_LEN_final_so_oligo_15PC cluster-specific/enriched transcripts
oligo <- read.csv(file = "./CR4_LEN_final_so_oligo_15PC_res0.1_genes-RNA.csv")

#Subset by oligodendrocyte cluster
oligo_C0 <- subset(oligo, select = c(gene), oligo$cluster == "0")
oligo_C1 <- subset(oligo, select = c(gene), oligo$cluster == "1")
oligo_C2 <- subset(oligo, select = c(gene), oligo$cluster == "2")
oligo_C3 <- subset(oligo, select = c(gene), oligo$cluster == "3")
oligo_C4 <- subset(oligo, select = c(gene), oligo$cluster == "4")

##Oligodendrocyte cluster 0##
#Convert gene names to ENSEMBL IDs
oligo_C0_ensembl <- mapIds(org.Hs.eg.db, keys = oligo_C0$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
#Verify updated ENSEMBL IDs for NAs
sum(is.na(oligo_C0_ensembl))
#Save gene_id and ensembl_id files
write.csv(oligo_C0, file='./Oligo_gene_id_C0.csv')
write.csv(oligo_C0_ensembl, file='./Oligo_ensembl_id_C0.csv')
#Update missing ENSEMBL IDs manually by checking against human database: http://grch37.ensembl.org/Homo_sapiens/Info/Index
oligo_C0_ensembl[c(5,15,20,34,35,37,51,61,76,78,80)] <- 
  c("ENSG00000280441","ENSG00000157654","ENSG00000279168","ENSG00000286749","ENSG00000274265",
    "ENSG00000251138","ENSG00000232931","ENSG00000226824","ENSG00000242593","ENSG00000259972",
    "ENSG00000250195")
#Convert to dataframe
oligo_C0_ensembl <- as.data.frame(oligo_C0_ensembl)
#Convert ENSEMBL IDs to entrezgene names
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=oligo_C0_ensembl$oligo_C0_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  oligo_C0_ensembl[match(annotLookup$ensembl_gene_id, oligo_C0_ensembl$oligo_C0_ensembl),],
  annotLookup)
#Run GO/pathway analysis
go_enrich_C0 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
#Save dataframe as csv file
C0_pathways <- as.data.frame(go_enrich_C0)
write.csv(C0_pathways, file='./CR4_LEN_final_so_oligo_15PC_C0_pathways.csv')

##Oligodendrocyte cluster 1##
oligo_C1_ensembl <- mapIds(org.Hs.eg.db, keys = oligo_C1$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(oligo_C1_ensembl))
write.csv(oligo_C1, file='./Oligo_gene_id_C1.csv')
write.csv(oligo_C1_ensembl, file='./Oligo_ensembl_id_C1.csv')
oligo_C1_ensembl[c(9,46,54,79,86,87,94,99,158,181,227,250,321)] <-
  c("ENSG00000132475","ENSG00000237949","ENSG00000255311","ENSG00000229944","ENSG00000285783",
    "ENSG00000229294","ENSG00000234527","ENSG00000254420","ENSG00000249042","ENSG00000233559",
    "ENSG00000286695","ENSG00000227544","ENSG00000207955")
oligo_C1_ensembl <- as.data.frame(oligo_C1_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=oligo_C1_ensembl$oligo_C1_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  oligo_C1_ensembl[match(annotLookup$ensembl_gene_id, oligo_C1_ensembl$oligo_C1_ensembl),],
  annotLookup)
go_enrich_C1 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C1_pathways <- as.data.frame(go_enrich_C1)
write.csv(C1_pathways, file='./CR4_LEN_final_so_oligo_15PC_C1_pathways.csv')

##Oligodendrocyte cluster 2##
oligo_C2_ensembl <- mapIds(org.Hs.eg.db, keys = oligo_C2$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(oligo_C2_ensembl))
write.csv(oligo_C2, file='./Oligo_gene_id_C2.csv')
write.csv(oligo_C2_ensembl, file='./Oligo_ensembl_id_C2.csv')
oligo_C2_ensembl[c(6,80,83,97,98,102,109,115,123,125,134,169,172,176,215,227,232,252,261,277)] <-
  c("ENSG00000286438","ENSG00000270207","ENSG00000231533","ENSG00000241345","ENSG00000287277",
    "ENSG00000258342","ENSG00000258844","ENSG00000271860","ENSG00000228793","ENSG00000251600",
    "ENSG00000224363","ENSG00000283098","ENSG00000227733","ENSG00000286749","ENSG00000258526",
    "ENSG00000251136","ENSG00000257139","ENSG00000280870","ENSG00000254180","ENSG00000237356")
oligo_C2_ensembl <- as.data.frame(oligo_C2_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=oligo_C2_ensembl$oligo_C2_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  oligo_C2_ensembl[match(annotLookup$ensembl_gene_id, oligo_C2_ensembl$oligo_C2_ensembl),],
  annotLookup)
go_enrich_C2 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C2_pathways <- as.data.frame(go_enrich_C2)
write.csv(C2_pathways, file='./CR4_LEN_final_so_oligo_15PC_C2_pathways.csv')

##Oligodendrocyte cluster 3##
oligo_C3_ensembl <- mapIds(org.Hs.eg.db, keys = oligo_C3$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(oligo_C3_ensembl))
write.csv(oligo_C3, file='./Oligo_gene_id_C3.csv')
write.csv(oligo_C3_ensembl, file='./Oligo_ensembl_id_C3.csv')
oligo_C3_ensembl[c(2,73,82,83,98,117,138,156,157,165,180,185,191,194,197,199,200,212)] <-
  c("ENSG00000286757","ENSG00000286288","ENSG00000228222","ENSG00000239268","ENSG00000287862",
    "ENSG00000125462","ENSG00000271860","ENSG00000198840","ENSG00000250938","ENSG00000198938",
    "ENSG00000198899","ENSG00000251136","ENSG00000198712","ENSG00000198727","ENSG00000126870",
    "ENSG00000198888","ENSG00000224078","ENSG00000198886")
oligo_C3_ensembl <- as.data.frame(oligo_C3_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=oligo_C3_ensembl$oligo_C3_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  oligo_C3_ensembl[match(annotLookup$ensembl_gene_id, oligo_C3_ensembl$oligo_C3_ensembl),],
  annotLookup)
go_enrich_C3 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C3_pathways <- as.data.frame(go_enrich_C3)
write.csv(C3_pathways, file='./CR4_LEN_final_so_oligo_15PC_C3_pathways.csv')

##Oligodendrocyte cluster 4##
oligo_C4_ensembl <- mapIds(org.Hs.eg.db, keys = oligo_C4$gene, keytype = "SYMBOL", column= ("ENSEMBL"))
sum(is.na(oligo_C4_ensembl))
write.csv(oligo_C4, file='./Oligo_gene_id_C4.csv')
write.csv(oligo_C4_ensembl, file='./Oligo_ensembl_id_C4.csv')
oligo_C4_ensembl[c(5,19,25,58,82,117,123,145,184,307,426,464,467)] <-
  c("ENSG00000225339","ENSG00000179300","ENSG00000140105","ENSG00000248587","ENSG00000106105",
    "ENSG00000134684","ENSG00000196305","ENSG00000166986","ENSG00000110619","ENSG00000269439",
    "ENSG00000133706","ENSG00000280441","ENSG00000240849")
oligo_C4_ensembl <- as.data.frame(oligo_C4_ensembl)
mart <- useMart("ensembl","hsapiens_gene_ensembl")
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=oligo_C4_ensembl$oligo_C4_ensembl,
  uniqueRows=TRUE)
annotLookup <- data.frame(
  oligo_C4_ensembl[match(annotLookup$ensembl_gene_id, oligo_C4_ensembl$oligo_C4_ensembl),],
  annotLookup)
go_enrich_C4 <- enrichGO(gene = annotLookup$entrezgene_id,
                         OrgDb = 'org.Hs.eg.db', 
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.10)
C4_pathways <- as.data.frame(go_enrich_C4)
write.csv(C4_pathways, file='./CR4_LEN_final_so_oligo_15PC_C4_pathways.csv')
