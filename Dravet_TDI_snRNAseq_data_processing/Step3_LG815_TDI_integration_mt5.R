
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
setwd("/athena/ganlab/scratch/lif4001/LG815_TDI/DF_2ndRound_mt5")

LG815_Ctrl_160 <- readRDS(file = "LG815_Ctrl_160_singlets_PCA.rds")
LG815_Ctrl_168 <- readRDS(file = "LG815_Ctrl_168_singlets_PCA.rds")
LG815_Ctrl_180 <- readRDS(file = "LG815_Ctrl_180_singlets_PCA.rds")
LG815_Ctrl_216 <- readRDS(file = "LG815_Ctrl_216_singlets_PCA.rds")

LG815_Ctrl_160[["Condition"]] = c('Scn1a: +/+')
LG815_Ctrl_168[["Condition"]] = c('Scn1a: +/+')
LG815_Ctrl_180[["Condition"]] = c('Scn1a: +/+')
LG815_Ctrl_216[["Condition"]] = c('Scn1a: +/+')

LG815_Ctrl_160[["Condition_Sex"]] = c('Scn1a: +/+_M')
LG815_Ctrl_168[["Condition_Sex"]] = c('Scn1a: +/+_F')
LG815_Ctrl_180[["Condition_Sex"]] = c('Scn1a: +/+_F')
LG815_Ctrl_216[["Condition_Sex"]] = c('Scn1a: +/+_M')

LG815_Ctrl_160[["Sample_Name"]] = c('Scn1a: +/+_1')
LG815_Ctrl_168[["Sample_Name"]] = c('Scn1a: +/+_2')
LG815_Ctrl_180[["Sample_Name"]] = c('Scn1a: +/+_3')
LG815_Ctrl_216[["Sample_Name"]] = c('Scn1a: +/+_4')

LG815_TDI_141 <- readRDS(file = "LG815_TDI_141_singlets_PCA.rds")
LG815_TDI_163 <- readRDS(file = "LG815_TDI_163_singlets_PCA.rds")
LG815_TDI_204 <- readRDS(file = "LG815_TDI_204_singlets_PCA.rds")
LG815_TDI_205 <- readRDS(file = "LG815_TDI_205_singlets_PCA.rds")

LG815_TDI_141[["Condition"]] = c('Scn1a: +/+; TDI')
LG815_TDI_163[["Condition"]] = c('Scn1a: +/+; TDI')
LG815_TDI_204[["Condition"]] = c('Scn1a: +/+; TDI')
LG815_TDI_205[["Condition"]] = c('Scn1a: +/+; TDI')

LG815_TDI_141[["Condition_Sex"]] = c('Scn1a: +/+; TDI_M')
LG815_TDI_163[["Condition_Sex"]] = c('Scn1a: +/+; TDI_F')
LG815_TDI_204[["Condition_Sex"]] = c('Scn1a: +/+; TDI_M')
LG815_TDI_205[["Condition_Sex"]] = c('Scn1a: +/+; TDI_F')

LG815_TDI_141[["Sample_Name"]] = c('Scn1a: +/+; TDI_1')
LG815_TDI_163[["Sample_Name"]] = c('Scn1a: +/+; TDI_2')
LG815_TDI_204[["Sample_Name"]] = c('Scn1a: +/+; TDI_3')
LG815_TDI_205[["Sample_Name"]] = c('Scn1a: +/+; TDI_4')

LG815_Ctrl_162 <- readRDS(file = "LG815_Ctrl_162_singlets_PCA.rds")
LG815_Ctrl_169 <- readRDS(file = "LG815_Ctrl_169_singlets_PCA.rds")
LG815_Ctrl_174 <- readRDS(file = "LG815_Ctrl_174_singlets_PCA.rds")
LG815_Ctrl_195 <- readRDS(file = "LG815_Ctrl_195_singlets_PCA.rds")

LG815_Ctrl_162[["Condition"]] = c('Scn1a: +/-')
LG815_Ctrl_169[["Condition"]] = c('Scn1a: +/-')
LG815_Ctrl_174[["Condition"]] = c('Scn1a: +/-')
LG815_Ctrl_195[["Condition"]] = c('Scn1a: +/-')

LG815_Ctrl_162[["Condition_Sex"]] = c('Scn1a: +/-_F')
LG815_Ctrl_169[["Condition_Sex"]] = c('Scn1a: +/-_M')
LG815_Ctrl_174[["Condition_Sex"]] = c('Scn1a: +/-_M')
LG815_Ctrl_195[["Condition_Sex"]] = c('Scn1a: +/-_F')

LG815_Ctrl_162[["Sample_Name"]] = c('Scn1a: +/-_1')
LG815_Ctrl_169[["Sample_Name"]] = c('Scn1a: +/-_2')
LG815_Ctrl_174[["Sample_Name"]] = c('Scn1a: +/-_3')
LG815_Ctrl_195[["Sample_Name"]] = c('Scn1a: +/-_4')

LG815_TDI_149 <- readRDS(file = "LG815_TDI_149_singlets_PCA.rds")
LG815_TDI_152 <- readRDS(file = "LG815_TDI_152_singlets_PCA.rds")
LG815_TDI_160 <- readRDS(file = "LG815_TDI_160_singlets_PCA.rds")
LG815_TDI_199 <- readRDS(file = "LG815_TDI_199_singlets_PCA.rds")

LG815_TDI_149[["Condition"]] = c('Scn1a: +/-; TDI')
LG815_TDI_152[["Condition"]] = c('Scn1a: +/-; TDI')
LG815_TDI_160[["Condition"]] = c('Scn1a: +/-; TDI')
LG815_TDI_199[["Condition"]] = c('Scn1a: +/-; TDI')

LG815_TDI_149[["Condition_Sex"]] = c('Scn1a: +/-; TDI_F')
LG815_TDI_152[["Condition_Sex"]] = c('Scn1a: +/-; TDI_M')
LG815_TDI_160[["Condition_Sex"]] = c('Scn1a: +/-; TDI_F')
LG815_TDI_199[["Condition_Sex"]] = c('Scn1a: +/-; TDI_M')

LG815_TDI_149[["Sample_Name"]] = c('Scn1a: +/-; TDI_1')
LG815_TDI_152[["Sample_Name"]] = c('Scn1a: +/-; TDI_2')
LG815_TDI_160[["Sample_Name"]] = c('Scn1a: +/-; TDI_3')
LG815_TDI_199[["Sample_Name"]] = c('Scn1a: +/-; TDI_4')

setwd("/athena/ganlab/scratch/lif4001/LG815_TDI/integration_mt5")

Condition_1 <- c(LG815_Ctrl_160, LG815_Ctrl_168, LG815_Ctrl_180, LG815_Ctrl_216)
anchors_Condition_1 <- FindIntegrationAnchors(object.list = Condition_1, dims = 1:30)
Condition_1_integrated <- IntegrateData(anchorset = anchors_Condition_1, dims = 1:30)
rm(LG815_Ctrl_160, LG815_Ctrl_168, LG815_Ctrl_180, LG815_Ctrl_216, Condition_1)

Condition_2 <- c(LG815_TDI_141, LG815_TDI_163, LG815_TDI_204, LG815_TDI_205)
anchors_Condition_2 <- FindIntegrationAnchors(object.list = Condition_2, dims = 1:30)
Condition_2_integrated <- IntegrateData(anchorset = anchors_Condition_2, dims = 1:30)
rm(LG815_TDI_141, LG815_TDI_163, LG815_TDI_204, LG815_TDI_205, Condition_2)

Condition_3 <- c(LG815_Ctrl_162, LG815_Ctrl_169, LG815_Ctrl_174, LG815_Ctrl_195)
anchors_Condition_3 <- FindIntegrationAnchors(object.list = Condition_3, dims = 1:30)
Condition_3_integrated <- IntegrateData(anchorset = anchors_Condition_3, dims = 1:30)
rm(LG815_Ctrl_162, LG815_Ctrl_169, LG815_Ctrl_174, LG815_Ctrl_195, Condition_3)

Condition_4 <- c(LG815_TDI_149, LG815_TDI_152, LG815_TDI_160, LG815_TDI_199)
anchors_Condition_4 <- FindIntegrationAnchors(object.list = Condition_4, dims = 1:30)
Condition_4_integrated <- IntegrateData(anchorset = anchors_Condition_4, dims = 1:30)
rm(LG815_TDI_149, LG815_TDI_152, LG815_TDI_160, LG815_TDI_199, Condition_4)



LG815_TDI <- c(Condition_1_integrated,Condition_2_integrated,Condition_3_integrated,Condition_4_integrated)
anchors_LG815_TDI <- FindIntegrationAnchors(object.list = LG815_TDI, dims = 1:30)
LG815_TDI_integrated <- IntegrateData(anchorset = anchors_LG815_TDI, dims = 1:30)
rm(Condition_1_integrated,Condition_2_integrated,Condition_3_integrated,Condition_4_integrated, LG815_TDI)

#saveRDS(LG815_TDI_integrated, file = "LG815_TDI_integrated.rds")

DefaultAssay(LG815_TDI_integrated) <- 'integrated'

# LG815_TDI_integrated <- NormalizeData(LG815_TDI_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# LG815_TDI_integrated <- FindVariableFeatures(LG815_TDI_integrated, selection.method = "vst", nfeatures = 3000)

LG815_TDI_integrated <- ScaleData(LG815_TDI_integrated, verbose = FALSE)
LG815_TDI_integrated <- RunPCA(LG815_TDI_integrated, features = VariableFeatures(object = LG815_TDI_integrated), verbose = FALSE)

LG815_TDI_integrated <- FindNeighbors(LG815_TDI_integrated, dims = 1:15)
LG815_TDI_integrated <- FindClusters(LG815_TDI_integrated, resolution = 0.1)
LG815_TDI_integrated <- RunUMAP(LG815_TDI_integrated, dims = 1: 15)

DefaultAssay(LG815_TDI_integrated) <- 'RNA'
LG815_TDI_integrated <- NormalizeData(LG815_TDI_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
LG815_TDI_integrated <- ScaleData(LG815_TDI_integrated, features = rownames(LG815_TDI_integrated))

#saveRDS(LG815_TDI_integrated, file = 'LG815_TDI_integrated_PCA_0.1.rds')
#LG815_TDI_integrated <- readRDS(file = "LG815_TDI_integrated_PCA_0.1.rds")

LG815_TDI_integrated <- JoinLayers(LG815_TDI_integrated)

LG815_TDI_integrated$Condition <- factor(x = LG815_TDI_integrated$Condition, levels = c("Scn1a: +/+","Scn1a: +/+; TDI","Scn1a: +/-", "Scn1a: +/-; TDI"))
LG815_TDI_integrated$Sample_Name <- factor(x = LG815_TDI_integrated$Sample_Name, levels = c("Scn1a: +/+_1","Scn1a: +/+_2","Scn1a: +/+_3","Scn1a: +/+_4",
                                                                                      "Scn1a: +/+; TDI_1","Scn1a: +/+; TDI_2","Scn1a: +/+; TDI_3","Scn1a: +/+; TDI_4",
                                                                                      "Scn1a: +/-_1","Scn1a: +/-_2","Scn1a: +/-_3","Scn1a: +/-_4",
                                                                                      "Scn1a: +/-; TDI_1","Scn1a: +/-; TDI_2","Scn1a: +/-; TDI_3","Scn1a: +/-; TDI_4"))

pdf("LG815_TDI_QC.pdf", width=9, height=6)
Idents(LG815_TDI_integrated) <- "Condition"
VlnPlot(object = LG815_TDI_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(LG815_TDI_integrated) <- "Sample_Name"
pdf("LG815_TDI_QC_Sample.pdf", width=18, height=6)

VlnPlot(object = LG815_TDI_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(LG815_TDI_integrated) <- "seurat_clusters"
pdf("LG815_TDI_integrated_umap.pdf", width=5, height=4)
DimPlot(LG815_TDI_integrated, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_integrated_umap_split_individual.pdf", width=12, height=12)
DimPlot(LG815_TDI_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
pdf("LG815_TDI_integrated_umap_split_Condition.pdf", width=6, height=6)
DimPlot(LG815_TDI_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

write.csv(table(LG815_TDI_integrated$seurat_clusters, LG815_TDI_integrated$Sample_Name), "LG815_TDI_cell_counts_cluster_by_sample.csv")

DefaultAssay(LG815_TDI_integrated) <- 'RNA'

LG815_TDI_markers <- FindAllMarkers(LG815_TDI_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(LG815_TDI_markers, "LG815_TDI_markers.csv")


saveRDS(LG815_TDI_integrated, file = 'LG815_TDI_integrated_PCA_0.1.rds')

#LG815_TDI_markers <- read.csv(file = "LG815_TDI_markers.csv", header=T,row.names =1)
#top5 <- LG815_TDI_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#top5$gene <- as.character(top5$gene)
#pdf("LG815_TDI_HeatMapTop5_0.1_new.pdf", width=24, height=16)
#DoHeatmap(LG815_TDI_integrated, features = top5$gene) + NoLegend()
#dev.off()

#Add marker genes

sig_EN<-c("Snap25","Erbb4","Slc6a1","Hpca","C1ql3","Cplx2","Prox1","Slc17a7", "Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", 
          "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r","Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("LG815_TDI_annotation_combine.pdf", width=12, height=5)
DotPlot(object = LG815_TDI_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

