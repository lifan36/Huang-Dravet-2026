
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/DVPV/DF_2ndRound")

DVPV_CTX_311 <- readRDS(file = "DVPV_CTX_311_singlets_PCA.rds")
DVPV_CTX_330 <- readRDS(file = "DVPV_CTX_330_singlets_PCA.rds")
DVPV_CTX_333 <- readRDS(file = "DVPV_CTX_333_singlets_PCA.rds")

DVPV_CTX_324 <- readRDS(file = "DVPV_CTX_324_singlets_PCA.rds")
DVPV_CTX_327 <- readRDS(file = "DVPV_CTX_327_singlets_PCA.rds")
DVPV_CTX_329 <- readRDS(file = "DVPV_CTX_329_singlets_PCA.rds")

DVPV_CTX_311[["Condition"]] = c('CTX_DV_WT')
DVPV_CTX_330[["Condition"]] = c('CTX_DV_WT')
DVPV_CTX_333[["Condition"]] = c('CTX_DV_WT')
DVPV_CTX_324[["Condition"]] = c('CTX_DV_KI')
DVPV_CTX_327[["Condition"]] = c('CTX_DV_KI')
DVPV_CTX_329[["Condition"]] = c('CTX_DV_KI')

DVPV_CTX_311[["Sample_Name"]] = c('CTX_DV_WT_1')
DVPV_CTX_330[["Sample_Name"]] = c('CTX_DV_WT_2')
DVPV_CTX_333[["Sample_Name"]] = c('CTX_DV_WT_3')
DVPV_CTX_324[["Sample_Name"]] = c('CTX_DV_KI_1')
DVPV_CTX_327[["Sample_Name"]] = c('CTX_DV_KI_2')
DVPV_CTX_329[["Sample_Name"]] = c('CTX_DV_KI_3')

DVPV_Hipp_311 <- readRDS(file = "DVPV_Hipp_311_singlets_PCA.rds")
DVPV_Hipp_330 <- readRDS(file = "DVPV_Hipp_330_singlets_PCA.rds")
DVPV_Hipp_333 <- readRDS(file = "DVPV_Hipp_333_singlets_PCA.rds")

DVPV_Hipp_324 <- readRDS(file = "DVPV_Hipp_324_singlets_PCA.rds")
DVPV_Hipp_327 <- readRDS(file = "DVPV_Hipp_327_singlets_PCA.rds")
DVPV_Hipp_329 <- readRDS(file = "DVPV_Hipp_329_singlets_PCA.rds")

DVPV_Hipp_311[["Condition"]] = c('Hipp_DV_WT')
DVPV_Hipp_330[["Condition"]] = c('Hipp_DV_WT')
DVPV_Hipp_333[["Condition"]] = c('Hipp_DV_WT')
DVPV_Hipp_324[["Condition"]] = c('Hipp_DV_KI')
DVPV_Hipp_327[["Condition"]] = c('Hipp_DV_KI')
DVPV_Hipp_329[["Condition"]] = c('Hipp_DV_KI')

DVPV_Hipp_311[["Sample_Name"]] = c('Hipp_DV_WT_1')
DVPV_Hipp_330[["Sample_Name"]] = c('Hipp_DV_WT_2')
DVPV_Hipp_333[["Sample_Name"]] = c('Hipp_DV_WT_3')
DVPV_Hipp_324[["Sample_Name"]] = c('Hipp_DV_KI_1')
DVPV_Hipp_327[["Sample_Name"]] = c('Hipp_DV_KI_2')
DVPV_Hipp_329[["Sample_Name"]] = c('Hipp_DV_KI_3')

setwd("/athena/ganlab/scratch/lif4001/DVPV/integration_both")

DV_WT_CTX <- c(DVPV_CTX_311, DVPV_CTX_330, DVPV_CTX_333)
anchors_DV_WT_CTX <- FindIntegrationAnchors(object.list = DV_WT_CTX, dims = 1:30)
DV_WT_CTX_integrated <- IntegrateData(anchorset = anchors_DV_WT_CTX, dims = 1:30)
rm(DVPV_CTX_311, DVPV_CTX_330, DVPV_CTX_333, DV_WT_CTX)

DV_KI_CTX <- c(DVPV_CTX_324, DVPV_CTX_327, DVPV_CTX_329)
anchors_DV_KI_CTX <- FindIntegrationAnchors(object.list = DV_KI_CTX, dims = 1:30)
DV_KI_CTX_integrated <- IntegrateData(anchorset = anchors_DV_KI_CTX, dims = 1:30)
rm(DVPV_CTX_324, DVPV_CTX_327, DVPV_CTX_329, DV_KI_CTX)

DV_WT_Hipp <- c(DVPV_Hipp_311, DVPV_Hipp_330, DVPV_Hipp_333)
anchors_DV_WT_Hipp <- FindIntegrationAnchors(object.list = DV_WT_Hipp, dims = 1:30)
DV_WT_Hipp_integrated <- IntegrateData(anchorset = anchors_DV_WT_Hipp, dims = 1:30)
rm(DVPV_Hipp_311, DVPV_Hipp_330, DVPV_Hipp_333, DV_WT_Hipp)

DV_KI_Hipp <- c(DVPV_Hipp_324, DVPV_Hipp_327, DVPV_Hipp_329)
anchors_DV_KI_Hipp <- FindIntegrationAnchors(object.list = DV_KI_Hipp, dims = 1:30)
DV_KI_Hipp_integrated <- IntegrateData(anchorset = anchors_DV_KI_Hipp, dims = 1:30)
rm(DVPV_Hipp_324, DVPV_Hipp_327, DVPV_Hipp_329, DV_KI_Hipp)

DVPV <- c(DV_WT_CTX_integrated,DV_KI_CTX_integrated, DV_WT_Hipp_integrated, DV_KI_Hipp_integrated)
anchors_DVPV <- FindIntegrationAnchors(object.list = DVPV, dims = 1:30)
DVPV_integrated <- IntegrateData(anchorset = anchors_DVPV, dims = 1:30)
rm(DV_WT_CTX_integrated,DV_KI_CTX_integrated, DV_WT_Hipp_integrated, DV_KI_Hipp_integrated, DVPV)

#saveRDS(DVPV_integrated, file = "DVPV_integrated.rds")

DefaultAssay(DVPV_integrated) <- 'integrated'

# DVPV_integrated <- NormalizeData(DVPV_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# DVPV_integrated <- FindVariableFeatures(DVPV_integrated, selection.method = "vst", nfeatures = 3000)

DVPV_integrated <- ScaleData(DVPV_integrated, verbose = FALSE)
DVPV_integrated <- RunPCA(DVPV_integrated, features = VariableFeatures(object = DVPV_integrated), verbose = FALSE)

DVPV_integrated <- FindNeighbors(DVPV_integrated, dims = 1:15)
DVPV_integrated <- FindClusters(DVPV_integrated, resolution = 0.1)
DVPV_integrated <- RunUMAP(DVPV_integrated, dims = 1: 15)

DefaultAssay(DVPV_integrated) <- 'RNA'
DVPV_integrated <- NormalizeData(DVPV_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
DVPV_integrated <- ScaleData(DVPV_integrated, features = rownames(DVPV_integrated))

#saveRDS(DVPV_integrated, file = 'DVPV_integrated_PCA_0.1.rds')
#DVPV_integrated <- readRDS(file = "DVPV_integrated_PCA_0.1.rds")

DVPV_integrated$Condition <- factor(x = DVPV_integrated$Condition, levels = c("CTX_DV_WT","CTX_DV_KI","Hipp_DV_WT","Hipp_DV_KI"))
DVPV_integrated$Sample_Name <- factor(x = DVPV_integrated$Sample_Name, levels = c("CTX_DV_WT_1","CTX_DV_WT_2","CTX_DV_WT_3",
                                                                                  "CTX_DV_KI_1","CTX_DV_KI_2","CTX_DV_KI_3",
                                                                                  "Hipp_DV_WT_1","Hipp_DV_WT_2","Hipp_DV_WT_3",
                                                                                  "Hipp_DV_KI_1","Hipp_DV_KI_2","Hipp_DV_KI_3"))

pdf("DVPV_QC.pdf", width=9, height=4)
Idents(DVPV_integrated) <- "Condition"
VlnPlot(object = DVPV_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(DVPV_integrated) <- "Sample_Name"
pdf("DVPV_QC_Sample.pdf", width=12, height=4)

VlnPlot(object = DVPV_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(DVPV_integrated) <- "seurat_clusters"
pdf("DVPV_integrated_umap.pdf", width=5, height=4)
DimPlot(DVPV_integrated, reduction = 'umap', label = T)
dev.off()
pdf("DVPV_integrated_umap_split_individual.pdf", width=9, height=12)
DimPlot(DVPV_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 3)
dev.off()
pdf("DVPV_integrated_umap_split_Condition.pdf", width=6, height=6)
DimPlot(DVPV_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

write.csv(table(DVPV_integrated$seurat_clusters, DVPV_integrated$Sample_Name), "DVPV_cell_counts_cluster_by_sample.csv")

DefaultAssay(DVPV_integrated) <- 'RNA'

DVPV_markers <- FindAllMarkers(DVPV_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(DVPV_markers, "DVPV_markers.csv")


saveRDS(DVPV_integrated, file = 'DVPV_integrated_PCA_0.1.rds')

DVPV_integrated <- readRDS(file = "DVPV_integrated_PCA_0.1.rds")
#DVPV_markers <- read.csv(file = "DVPV_markers.csv", header=T,row.names =1)
#top5 <- DVPV_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#top5$gene <- as.character(top5$gene)
#pdf("DVPV_HeatMapTop5_0.1_new.pdf", width=24, height=16)
#DoHeatmap(DVPV_integrated, features = top5$gene) + NoLegend()
#dev.off()

#Add marker genes

sig_EN<-c("Snap25","Vsnl1","Dnm1","Pde1a","Cplx1","Nap1l5","Erbb4","Slc6a1","Hpca","Grin2b","Rasgrp1","Ppp3r1","Camk4","Brinp1",
          "Nrgn","C1ql3","Cplx2","Rbfox1","Prox1","Slc17a7", "Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", 
          "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r","Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("DVPV_annotation_combine.pdf", width=20, height=5)
DotPlot(object = DVPV_integrated, features = markers.to.plot) + RotatedAxis()
dev.off()

