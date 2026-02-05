library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MIN)
library(EnhancedVolcano)


MG <- subset(DV,idents = c('5','26') )
saveRDS(MG,"DVPV_MG_subset.rds")
AST <- subset(DV,idents = c('1','25') )
saveRDS(MG,"DVPV_AST_subset.rds")

DefaultAssay(MG) <- 'integrated'
all.genes <- rownames(MG)
MG<- ScaleData(MG, features = all.genes)
MG<- FindVariableFeatures(object = MG)
MG<- RunPCA(MG, features = VariableFeatures(object = MG))

ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.3)
MG <- RunUMAP(MG, dims = 1: 15)
DimPlot(MG, reduction = 'umap', label = T, split.by = "Condition")

MG2 <- subset(MG, idents = c("0",'1','2','3'))

DefaultAssay(MG2) <- 'integrated'
all.genes <- rownames(MG2)
MG2<- ScaleData(MG2, features = all.genes)
MG2<- FindVariableFeatures(object = MG2)
MG2<- RunPCA(MG2, features = VariableFeatures(object = MG2))

MG2 <- FindNeighbors(MG2, dims = 1:15)
MG2 <- FindClusters(MG2, resolution = 0.3)
MG2 <- RunUMAP(MG2, dims = 1: 15)
DimPlot(MG2, reduction = 'umap', label = T, split.by = "Condition")

MG2 <- subset(MG2, idents = c("0",'1','2','3'))
DefaultAssay(MG2) <- 'integrated'
all.genes <- rownames(MG2)
MG2<- ScaleData(MG2, features = all.genes)
MG2<- FindVariableFeatures(object = MG2)
MG2<- RunPCA(MG2, features = VariableFeatures(object = MG2))

MG2 <- FindNeighbors(MG2, dims = 1:15)
MG2 <- FindClusters(MG2, resolution = 0.3)
MG2 <- RunUMAP(MG2, dims = 1: 15)
DimPlot(MG2, reduction = 'umap', label = T, split.by = "Condition")

# rename cluster
n <- dim(table(MG2@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG2@active.ident <- plyr::mapvalues(x = MG2@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG2@active.ident <- factor(MG2@active.ident, levels=1:n)
saveRDS(MG2, file = 'DVPV_MG_reclusted_res0.3.rds')
