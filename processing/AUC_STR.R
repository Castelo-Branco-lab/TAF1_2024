library(Seurat)
library(Augur)
set.seed(1234)


CRL.seurat.clean.h   <-  readRDS("snRNA_wt_taf1ex38del_seuratv4.rds")




CRL.seurat.clean.h$STR_cond <- paste(CRL.seurat.clean.h$sample, CRL.seurat.clean.h$Celltypes_Striatum, sep = "_")
CRL.seurat.clean.h$SCT08_cond <- paste(CRL.seurat.clean.h$sample, CRL.seurat.clean.h$SCT_snn_res.0.8, sep = "_")
CRL.seurat.clean.h$SCT1_cond <- paste(CRL.seurat.clean.h$sample, CRL.seurat.clean.h$SCT_snn_res.1, sep = "_")
table(CRL.seurat.clean.h$SCT1_cond)
table(CRL.seurat.clean.h$SCT08_cond)
table(CRL.seurat.clean.h$STR_cond)

Idents(CRL.seurat.clean.h) <- "Batch"

CRL.seurat.clean.h2 <- subset( CRL.seurat.clean.h , ident=("B1") )

CRL.seurat.clean.h <- CRL.seurat.clean.h2

#same results when downsampling
#Idents(CRL.seurat.clean.h) <- "STR_cond"
#CRL_subset <- subset(CRL.seurat.qchf , idents=c( "KO_s3") , invert=T)

#CRL.object.STR <- subset( CRL.seurat.clean.h , downsample = 200 )

#Idents(CRL.seurat.clean.h) <- "SCT08_cond"
#CRL.object.08 <- subset( CRL.seurat.clean.h , downsample = 200 )

#Idents(CRL.seurat.clean.h) <- "SCT1_cond"
#CRL.object.1 <- subset( CRL.seurat.clean.h , downsample = 200 )




#rm(CRL.seurat.clean.h)
#rm(CRL_subset)

#sct 0.8
SCT08_STR_auc_subset = calculate_auc( CRL.seurat.clean.h2 , cell_type_col = "SCT_snn_res.0.8", label_col = "sample"  , n_threads = 2 )


saveRDS(SCT08_STR_auc_subset , file="SCT08_STR.rds")


#AMM STR
AMM_STR_auc_subset2 = calculate_auc(CRL.seurat.clean.h2 , cell_type_col = "Celltypes_Striatum", label_col = "sample"  , n_threads = 2 )


saveRDS(AMM_STR_auc_subset2 , file="AMM_STR.rds")



#sct 0.8
SCT08_STR_auc_subset3 = calculate_auc(CRL.seurat.clean.h2 , cell_type_col = "SCT_snn_res.1", label_col = "sample"  , n_threads = 2 )


saveRDS(SCT08_STR_auc_subset3 , file="SCT1_STR.rds")


####
auc_allrep_SCT1 <- readRDS("SCT1_STR.rds")
auc_allrep_AMM <- readRDS("AMM_STR.rds")
auc_allrep_SCT08 <- readRDS("SCT08_STR.rds")


(auc_allrep_AMM$AUC)
###here the used AUC microglia and others significant is resolution 0.8 and 1
DimPlot_scCustom( CRL.seurat.clean.h , group.by = "Celltypes_Striatum")
(auc_allrep_SCT08$AUC) # 
DimPlot_scCustom( CRL.seurat.clean.h , group.by = "SCT_snn_res.0.8")
(auc_allrep_SCT1$AUC) # 

write.csv(auc_allrep_SCT1$AUC , file = "SCT1_STR_AUC.csv")
