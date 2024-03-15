library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(scCustomize)
set.seed(1234)


CRL_KO1 <- Read10X(data.dir = "CRL_scRNAseq/mm10_P29160_1001_RNA_CR710/outs/filtered_feature_bc_matrix/")

CRL_WT1 <- Read10X(data.dir = "CRL_scRNAseq/mm10_P29160_1002_RNA_CR710/outs/filtered_feature_bc_matrix/")

CRL_WT2 <- Read10X(data.dir = "CRL_scRNAseq/mm10_P29160_1003_RNA_CR710/outs/filtered_feature_bc_matrix/")

CRL_KO2 <- Read10X(data.dir = "CRL_scRNAseq/mm10_P29160_1004_RNA_CR710/outs/filtered_feature_bc_matrix/")

CRL_WT3 <- Read10X(data.dir = "CRL_scRNAseq/mm10_P29160_1005_RNA_CR710/outs/filtered_feature_bc_matrix/")

CRL_KO3 <- Read10X(data.dir = "CRL_scRNAseq/mm10_P29160_1006_RNA_CR710/outs/filtered_feature_bc_matrix/")





KO3.seurat <- CreateSeuratObject(counts = CRL_KO3  , project = "CRL_KO3")
KO2.seurat <- CreateSeuratObject(counts = CRL_KO2  , project = "CRL_KO2")
KO1.seurat <- CreateSeuratObject(counts = CRL_KO1  , project = "CRL_KO1")
WT3.seurat <- CreateSeuratObject(counts = CRL_WT3  , project = "CRL_WT3")
WT2.seurat <- CreateSeuratObject(counts = CRL_WT2  , project = "CRL_WT2")
WT1.seurat <- CreateSeuratObject(counts = CRL_WT1  , project = "CRL_WT1")


KO3.seurat$sample <- "KO"
KO1.seurat$sample <- "KO"
KO2.seurat$sample <- "KO"

WT1.seurat$sample <- "WT"
WT2.seurat$sample <- "WT"
WT3.seurat$sample <- "WT"


KO1.seurat$replicate <- "s1"
KO2.seurat$replicate <- "s2"
KO3.seurat$replicate <- "s3"

WT1.seurat$replicate <- "s1"
WT2.seurat$replicate <- "s2"
WT3.seurat$replicate <- "s3"


KO1.seurat$rep_sample <- "KO_s1"
KO2.seurat$rep_sample <- "KO_s2"
KO3.seurat$rep_sample <- "KO_s3"

WT1.seurat$rep_sample <- "WT_s1"
WT2.seurat$rep_sample <- "WT_s2"
WT3.seurat$rep_sample <- "WT_s3"


KO1.seurat <- Add_Mito_Ribo_Seurat(seurat_object = KO1.seurat, species = "Mouse") 
KO2.seurat <- Add_Mito_Ribo_Seurat(seurat_object = KO2.seurat, species = "Mouse") 
KO3.seurat <- Add_Mito_Ribo_Seurat(seurat_object = KO3.seurat, species = "Mouse") 
WT1.seurat <- Add_Mito_Ribo_Seurat(seurat_object = WT1.seurat, species = "Mouse") 
WT2.seurat <- Add_Mito_Ribo_Seurat(seurat_object = WT2.seurat, species = "Mouse") 
WT3.seurat <- Add_Mito_Ribo_Seurat(seurat_object = WT3.seurat, species = "Mouse") 

KO1.seurat <- PercentageFeatureSet( KO1.seurat , pattern = "^mt\\-", col.name = "percent_mito")
KO2.seurat <- PercentageFeatureSet( KO2.seurat , pattern = "^mt\\-", col.name = "percent_mito")
KO3.seurat <- PercentageFeatureSet( KO3.seurat , pattern = "^mt\\-", col.name = "percent_mito")
WT1.seurat <- PercentageFeatureSet( WT1.seurat , pattern = "^mt\\-", col.name = "percent_mito")
WT2.seurat <- PercentageFeatureSet( WT2.seurat , pattern = "^mt\\-", col.name = "percent_mito")
WT3.seurat <- PercentageFeatureSet( WT3.seurat , pattern = "^mt\\-", col.name = "percent_mito")


KO1.seurat$Batch <- "B1"
KO2.seurat$Batch <- "B1"
KO3.seurat$Batch <- "B2"
WT1.seurat$Batch <- "B1"
WT2.seurat$Batch <- "B1"
WT3.seurat$Batch <- "B2"
CRL.seurat.qc <- merge( x= KO1.seurat , y =  c( KO2.seurat , KO3.seurat , WT1.seurat , WT2.seurat , WT3.seurat ) )

#median stats
median_stats_raw <- Median_Stats(seurat_object = CRL.seurat.qc, group_by_var = "rep_sample")
median_stats_raw



WT1.seurat.qc <- subset(WT1.seurat, subset =  percent_mito < 2 & nFeature_RNA > 50 & nFeature_RNA < 5000 & nCount_RNA < 10000 & nCount_RNA > 5 )
VlnPlot( WT1.seurat.qc , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3 , pt.size = 0)


Median_Stats(WT1.seurat.qc)


WT2.seurat.qc <- subset(WT2.seurat, subset =  percent_mito < 2 & nFeature_RNA > 50 & nFeature_RNA < 8000 & nCount_RNA < 10000 & nCount_RNA > 5 )
VlnPlot( WT2.seurat.qc , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3 , pt.size = 0)


Median_Stats(WT2.seurat.qc)


WT3.seurat.qc <- subset(WT3.seurat, subset =  percent_mito < 2 & nFeature_RNA > 50 & nFeature_RNA < 8000 & nCount_RNA < 20000)
VlnPlot( WT3.seurat.qc , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3 , pt.size = 0)


Median_Stats(WT3.seurat.qc)


KO1.seurat.qc <- subset(KO1.seurat, subset =  percent_mito < 2 & nFeature_RNA > 50 & nFeature_RNA < 10000 & nCount_RNA < 8000 & nCount_RNA > 5 )
VlnPlot(KO1.seurat.qc , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3 , pt.size = 0)


Median_Stats(KO1.seurat.qc)


KO2.seurat.qc <- subset(KO2.seurat, subset =  percent_mito < 2 & nFeature_RNA > 50 & nFeature_RNA < 6000 & nCount_RNA < 50000 & nCount_RNA > 5 )
VlnPlot( KO2.seurat.qc , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3 , pt.size = 0)


Median_Stats(KO2.seurat.qc)




# change to 400 nfeature
KO3.seurat.qc <- subset(KO3.seurat, subset =  percent_mito < 2 & nFeature_RNA > 50 & nFeature_RNA < 4000 & nCount_RNA < 20000)
VlnPlot( KO3.seurat.qc , features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3 , pt.size = 0)


Median_Stats(KO3.seurat.qc)

KO1.seurat.qc$Batch <- "B1"
KO2.seurat.qc$Batch <- "B1"
KO3.seurat.qc$Batch <- "B2"
WT1.seurat.qc$Batch <- "B1"
WT2.seurat.qc$Batch <- "B1"
WT3.seurat.qc$Batch <- "B2"
CRL.seurat.qc <- merge( x= KO1.seurat.qc , y =  c( KO2.seurat.qc , KO3.seurat.qc , WT1.seurat.qc , WT2.seurat.qc , WT3.seurat.qc ) )
#saveRDS( CRL.seurat.qc , file="/date/gcb/gcb_EA/scRNAseq_CR_WTKO_LK_sEAE_MK_multi_LB_WTKO_140623/Processed_data/CRL_scRNAseq/CRL_striatum_KOWT_RNA_seurat_021023.rds")
KO1.seurat$Batch <- "B1"
KO2.seurat$Batch <- "B1"
KO3.seurat$Batch <- "B2"
WT1.seurat$Batch <- "B1"
WT2.seurat$Batch <- "B1"
WT3.seurat$Batch <- "B2"
CRL.seurat.qc <- merge( x= KO1.seurat , y =  c( KO2.seurat , KO3.seurat , WT1.seurat , WT2.seurat , WT3.seurat ) )



###regress out hormonal genes
cd_features <- list(c("Prl" , "Gh" ,  "Pou1f1" , "Tshb"))
CRL.seurat.qc <- AddModuleScore(
  object = CRL.seurat.qc,
  features = cd_features,
  ctrl = 4,
  nbin = 20, ## the default nbin=24 does not work. not enough data
  name = 'Horm_score_all'
)

summary(CRL.seurat.qc$Horm_score_all1)

median_stats <- Median_Stats(seurat_object = CRL.seurat.qc , group_by_var = "rep_sample", median_var = "Horm_score_all1")
median_stats

Plot_Median_Other(seurat_object = CRL.seurat.qc , median_var = "Horm_score_all1", group_by = "rep_sample")

#remove cells with Horm module > 0
QC_Plot_GenevsFeature(seurat_object = CRL.seurat.qc, feature1 = "Horm_score_all1", low_cutoff_gene = 800,
                      high_cutoff_gene = 5500, high_cutoff_feature = 20 , group.by = "rep_sample")

QC_Plot_GenevsFeature(seurat_object = CRL.seurat.qc, feature1 = "nCount_RNA", low_cutoff_gene = 800,
                      high_cutoff_gene = 5500, high_cutoff_feature = 20 , group.by = "rep_sample")


#modules for hormonal genes

CRL.seurat.qcfilterout <- CRL.seurat.qc


Plot_Median_Other(seurat_object = CRL.seurat.qc , median_var = "Horm_score_all1", group_by = "rep_sample")

#remove cells with Horm module > 0
QC_Plot_GenevsFeature(seurat_object = CRL.seurat.qc, feature1 = "Horm_score_all1", low_cutoff_gene = 800,
                      high_cutoff_gene = 5500, high_cutoff_feature = 20 , group.by = "rep_sample")

CRL.seurat.qcfilterout$orig.ident[CRL.seurat.qcfilterout$Horm_score_all1 > 0   ] <- "CONT_Horm"



Idents(CRL.seurat.qcfilterout) <- "orig.ident"
CRL.seurat.filterout <- subset(CRL.seurat.qcfilterout, idents=("CONT_Horm") , invert=T)


#read the qc seurat object and filter out cells
CRL.seurat.qc$cellIDs <- rownames(CRL.seurat.qc@meta.data)
CRL.seurat.qc.filter <- subset(CRL.seurat.qc ,  subset = cellIDs %in% rownames(CRL.seurat.filterout@meta.data) )

QC_Plot_GenevsFeature(seurat_object = CRL.seurat.qc.filter, feature1 = "Horm_score_all1", low_cutoff_gene = 800,
                      high_cutoff_gene = 5500, high_cutoff_feature = 20 , group.by = "rep_sample")
rm(CRL.seurat.filterout)

###########
#dim reduction and lcustering for the fitlered object

# Perform log-normalization and feature selection, as well as SCT normalization on global object
CRL.seurat.qch <- CRL.seurat.qc.filter %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  SCTransform(vars.to.regress = c("percent_mito"))

# Calculate PCs using variable features determined by SCTransform (3000 by default)
CRL.seurat.qch <- RunPCA(CRL.seurat.qch, assay = "SCT", npcs = 50)
library(harmony)
CRL.seurat.qch <- RunHarmony(CRL.seurat.qch, 
                             group.by.vars = c("rep_sample" , "Batch" , "sample"), 
                             reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

CRL.seurat.qch <- RunUMAP(CRL.seurat.qch, reduction = "harmony", assay = "SCT", dims = 1:50)

CRL.seurat.qch <- FindNeighbors(object = CRL.seurat.qch, reduction = "harmony")
CRL.seurat.qch <- FindClusters(CRL.seurat.qch, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))

####Muóz-Manchado et al reference
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97478

AMM_STR_counts <-  read.table("MunozManchado_GSE97478_STR/GSE97478_Munoz-Manchado_et_al_molecule_count.txt" , row.names=NULL)

head(AMM_STR_counts)

#is first 6 rows
first_five_rows <- AMM_STR_counts %>% as.data.frame() %>%  slice(1:6) 
transposed_matrix <- t(first_five_rows)
transposed_df <- as.data.frame(transposed_matrix)
colnames(transposed_df) <- 1:ncol(transposed_df)

colnames(transposed_df)
library(janitor)
final_annot <-  transposed_df %>% as.data.frame() %>% row_to_names(row_number = 1)

dim(AMM_STR_counts)

AMM_matrix <- AMM_STR_counts %>% as.data.frame() %>%  slice(7:12942)  

colnames(AMM_matrix) <-  colnames(AMM_STR_counts)

row_amm <- AMM_matrix[,1]
AMM_matrix <- AMM_matrix[,-1]
rownames(AMM_matrix) <- make.unique(row_amm)

head(AMM_matrix)

#####

AMM_seur <- CreateSeuratObject(counts = AMM_matrix , project = 'AMM_STR', meta.data = final_annot)


AMM_seur <- AMM_seur %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  SCTransform() %>%
  RunPCA( assay="SCT" , npcs=30) %>%
  RunUMAP( reduction = "pca", assay = "SCT", dims = 1:30) %>%
  FindNeighbors( reduction = "pca") %>% 
  FindClusters( resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))

##label transfer

Idents(AMM_seur) <- "cluster"
DefaultAssay(AMM_seur) <- "SCT"

transfer.anchors <- FindTransferAnchors(
  reference = AMM_seur,
  query = CRL.seurat.qch,
  reference.assay = 'SCT',
  query.assay = 'SCT',
  reduction = 'cca',
  
  k.filter = NA
)

predicted.id  <- TransferData(
  anchorset = transfer.anchors, 
  refdata = AMM_seur$cluster,
  dims = 1:30, weight.reduction = "cca")

CRL.seurat.qch  <- AddMetaData(
  object = CRL.seurat.qch,
  metadata = predicted.id
)
CRL.seurat.qch$Celltypes_Striatum <- CRL.seurat.qch$predicted.id
options(ggrepel.max.overlaps=Inf)

DimPlot(  CRL.seurat.qch , group.by = 'predicted.id', label = TRUE , repel = T) + NoLegend() + ggplot2::ggtitle("")


table(  CRL.seurat.qch$SCT_snn_res.1 , CRL.seurat.qch$predicted.id)




unique(CRL.seurat.qch$Celltypes_Striatum)
table(  CRL.seurat.qch$SCT_snn_res.1 , CRL.seurat.qch$Celltypes_Striatum)
###################################

#save object
#saveRDS(CRL.seurat.qch , file="snRNA_wt_taf1ex38del_seuratv4.rds")