
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
setwd("/athena/ganlab/scratch/lif4001/LG815C/DF_2ndRound_mt5")

LG815_C_164 <- readRDS(file = "LG815_C_164_singlets_PCA.rds")
LG815_C_29 <- readRDS(file = "LG815_C_29_singlets_PCA.rds")
LG815_C_35 <- readRDS(file = "LG815_C_35_singlets_PCA.rds")
LG815_C_64 <- readRDS(file = "LG815_C_64_singlets_PCA.rds")

LG815_C_164[["Condition"]] = c('Scn1a: +/+; cGAS: +/+')
LG815_C_29[["Condition"]] = c('Scn1a: +/+; cGAS: +/+')
LG815_C_35[["Condition"]] = c('Scn1a: +/+; cGAS: +/+')
LG815_C_64[["Condition"]] = c('Scn1a: +/+; cGAS: +/+')

LG815_C_164[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/+_1')
LG815_C_29[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/+_2')
LG815_C_35[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/+_3')
LG815_C_64[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/+_4')

LG815_C_172 <- readRDS(file = "LG815_C_172_singlets_PCA.rds")
LG815_C_30 <- readRDS(file = "LG815_C_30_singlets_PCA.rds")
LG815_C_31 <- readRDS(file = "LG815_C_31_singlets_PCA.rds")
LG815_C_33 <- readRDS(file = "LG815_C_33_singlets_PCA.rds")

LG815_C_172[["Condition"]] = c('Scn1a: +/+; cGAS: +/-')
LG815_C_30[["Condition"]] = c('Scn1a: +/+; cGAS: +/-')
LG815_C_31[["Condition"]] = c('Scn1a: +/+; cGAS: +/-')
LG815_C_33[["Condition"]] = c('Scn1a: +/+; cGAS: +/-')

LG815_C_172[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/-_1')
LG815_C_30[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/-_2')
LG815_C_31[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/-_3')
LG815_C_33[["Sample_Name"]] = c('Scn1a: +/+; cGAS: +/-_4')

LG815_C_106 <- readRDS(file = "LG815_C_106_singlets_PCA.rds")
LG815_C_121 <- readRDS(file = "LG815_C_121_singlets_PCA.rds")
LG815_C_132 <- readRDS(file = "LG815_C_132_singlets_PCA.rds")
LG815_C_95 <- readRDS(file = "LG815_C_95_singlets_PCA.rds")

LG815_C_106[["Condition"]] = c('Scn1a: +/-; cGAS: +/+')
LG815_C_121[["Condition"]] = c('Scn1a: +/-; cGAS: +/+')
LG815_C_132[["Condition"]] = c('Scn1a: +/-; cGAS: +/+')
LG815_C_95[["Condition"]] = c('Scn1a: +/-; cGAS: +/+')

LG815_C_106[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/+_1')
LG815_C_121[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/+_2')
LG815_C_132[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/+_3')
LG815_C_95[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/+_4')

LG815_C_122 <- readRDS(file = "LG815_C_122_singlets_PCA.rds")
LG815_C_139 <- readRDS(file = "LG815_C_139_singlets_PCA.rds")
LG815_C_153 <- readRDS(file = "LG815_C_153_singlets_PCA.rds")
LG815_C_159 <- readRDS(file = "LG815_C_159_singlets_PCA.rds")

LG815_C_122[["Condition"]] = c('Scn1a: +/-; cGAS: +/-')
LG815_C_139[["Condition"]] = c('Scn1a: +/-; cGAS: +/-')
LG815_C_153[["Condition"]] = c('Scn1a: +/-; cGAS: +/-')
LG815_C_159[["Condition"]] = c('Scn1a: +/-; cGAS: +/-')

LG815_C_122[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/-_1')
LG815_C_139[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/-_2')
LG815_C_153[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/-_3')
LG815_C_159[["Sample_Name"]] = c('Scn1a: +/-; cGAS: +/-_4')

setwd("/athena/ganlab/scratch/lif4001/LG815C/integration_mt5")

Condition_1 <- c(LG815_C_164, LG815_C_29, LG815_C_35, LG815_C_64)
anchors_Condition_1 <- FindIntegrationAnchors(object.list = Condition_1, dims = 1:30)
Condition_1_integrated <- IntegrateData(anchorset = anchors_Condition_1, dims = 1:30)
rm(LG815_C_164, LG815_C_29, LG815_C_35, LG815_C_64, Condition_1)

Condition_2 <- c(LG815_C_172, LG815_C_30, LG815_C_31, LG815_C_33)
anchors_Condition_2 <- FindIntegrationAnchors(object.list = Condition_2, dims = 1:30)
Condition_2_integrated <- IntegrateData(anchorset = anchors_Condition_2, dims = 1:30)
rm(LG815_C_172, LG815_C_30, LG815_C_31, LG815_C_33, Condition_2)

Condition_3 <- c(LG815_C_106, LG815_C_121, LG815_C_132, LG815_C_95)
anchors_Condition_3 <- FindIntegrationAnchors(object.list = Condition_3, dims = 1:30)
Condition_3_integrated <- IntegrateData(anchorset = anchors_Condition_3, dims = 1:30)
rm(LG815_C_106, LG815_C_121, LG815_C_132, LG815_C_95, Condition_3)

Condition_4 <- c(LG815_C_122, LG815_C_139, LG815_C_153, LG815_C_159)
anchors_Condition_4 <- FindIntegrationAnchors(object.list = Condition_4, dims = 1:30)
Condition_4_integrated <- IntegrateData(anchorset = anchors_Condition_4, dims = 1:30)
rm(LG815_C_122, LG815_C_139, LG815_C_153, LG815_C_159, Condition_4)



LG815C <- c(Condition_1_integrated,Condition_2_integrated,Condition_3_integrated,Condition_4_integrated)
anchors_LG815C <- FindIntegrationAnchors(object.list = LG815C, dims = 1:30)
LG815C_integrated <- IntegrateData(anchorset = anchors_LG815C, dims = 1:30)
rm(Condition_1_integrated,Condition_2_integrated,Condition_3_integrated,Condition_4_integrated, LG815C)

#saveRDS(LG815C_integrated, file = "LG815C_integrated.rds")

DefaultAssay(LG815C_integrated) <- 'integrated'

# LG815C_integrated <- NormalizeData(LG815C_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# LG815C_integrated <- FindVariableFeatures(LG815C_integrated, selection.method = "vst", nfeatures = 3000)

LG815C_integrated <- ScaleData(LG815C_integrated, verbose = FALSE)
LG815C_integrated <- RunPCA(LG815C_integrated, features = VariableFeatures(object = LG815C_integrated), verbose = FALSE)

LG815C_integrated <- FindNeighbors(LG815C_integrated, dims = 1:15)
LG815C_integrated <- FindClusters(LG815C_integrated, resolution = 0.1)
LG815C_integrated <- RunUMAP(LG815C_integrated, dims = 1: 15)

DefaultAssay(LG815C_integrated) <- 'RNA'
LG815C_integrated <- NormalizeData(LG815C_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
LG815C_integrated <- ScaleData(LG815C_integrated, features = rownames(LG815C_integrated))

#saveRDS(LG815C_integrated, file = 'LG815C_integrated_PCA_0.1.rds')
#LG815C_integrated <- readRDS(file = "LG815C_integrated_PCA_0.1.rds")

LG815C_integrated <- JoinLayers(LG815C_integrated)

LG815C_integrated$Condition <- factor(x = LG815C_integrated$Condition, levels = c("Scn1a: +/+; cGAS: +/+","Scn1a: +/+; cGAS: +/-","Scn1a: +/-; cGAS: +/+", "Scn1a: +/-; cGAS: +/-"))
LG815C_integrated$Sample_Name <- factor(x = LG815C_integrated$Sample_Name, levels = c("Scn1a: +/+; cGAS: +/+_1","Scn1a: +/+; cGAS: +/+_2","Scn1a: +/+; cGAS: +/+_3","Scn1a: +/+; cGAS: +/+_4",
                                                                                      "Scn1a: +/+; cGAS: +/-_1","Scn1a: +/+; cGAS: +/-_2","Scn1a: +/+; cGAS: +/-_3","Scn1a: +/+; cGAS: +/-_4",
                                                                                      "Scn1a: +/-; cGAS: +/+_1","Scn1a: +/-; cGAS: +/+_2","Scn1a: +/-; cGAS: +/+_3","Scn1a: +/-; cGAS: +/+_4",
                                                                                      "Scn1a: +/-; cGAS: +/-_1","Scn1a: +/-; cGAS: +/-_2","Scn1a: +/-; cGAS: +/-_3","Scn1a: +/-; cGAS: +/-_4"))

pdf("LG815C_QC.pdf", width=9, height=6)
Idents(LG815C_integrated) <- "Condition"
VlnPlot(object = LG815C_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(LG815C_integrated) <- "Sample_Name"
pdf("LG815C_QC_Sample.pdf", width=18, height=6)

VlnPlot(object = LG815C_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(LG815C_integrated) <- "seurat_clusters"
pdf("LG815C_integrated_umap.pdf", width=5, height=4)
DimPlot(LG815C_integrated, reduction = 'umap', label = T)
dev.off()
pdf("LG815C_integrated_umap_split_individual.pdf", width=12, height=12)
DimPlot(LG815C_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
pdf("LG815C_integrated_umap_split_Condition.pdf", width=6, height=6)
DimPlot(LG815C_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

write.csv(table(LG815C_integrated$seurat_clusters, LG815C_integrated$Sample_Name), "LG815C_cell_counts_cluster_by_sample.csv")

DefaultAssay(LG815C_integrated) <- 'RNA'

LG815C_markers <- FindAllMarkers(LG815C_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(LG815C_markers, "LG815C_markers.csv")


saveRDS(LG815C_integrated, file = 'LG815C_integrated_PCA_0.1.rds')

#LG815C_markers <- read.csv(file = "LG815C_markers.csv", header=T,row.names =1)
#top5 <- LG815C_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#top5$gene <- as.character(top5$gene)
#pdf("LG815C_HeatMapTop5_0.1_new.pdf", width=24, height=16)
#DoHeatmap(LG815C_integrated, features = top5$gene) + NoLegend()
#dev.off()

#Add marker genes

sig_EN<-c("Snap25","Vsnl1","Dnm1","Pde1a","Cplx1","Nap1l5","Erbb4","Slc6a1","Hpca","Grin2b","Rasgrp1","Ppp3r1","Camk4","Brinp1",
          "Nrgn","C1ql3","Cplx2","Rbfox1","Prox1","Slc17a7", "Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", 
          "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r","Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("LG815C_annotation_combine.pdf", width=20, height=5)
DotPlot(object = LG815C_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

