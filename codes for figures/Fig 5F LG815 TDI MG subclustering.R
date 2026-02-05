library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MIN)
library(EnhancedVolcano)

Idents(LG815_TDI) <- "celltype"
levels(LG815_TDI$celltype)
MG <- subset(LG815_TDI,idents = c('MG'))

DefaultAssay(MG) <- 'integrated'
all.genes <- rownames(MG)
MG<- ScaleData(MG, features = all.genes)
MG<- FindVariableFeatures(object = MG)
MG<- RunPCA(MG, features = VariableFeatures(object = MG))

ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:10)
MG <- FindClusters(MG, resolution = 0.3)
MG <- RunUMAP(MG, dims = 1:10)
dev.off()

DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()


DefaultAssay(MG) <- 'RNA'
VlnPlot(MG, features = c('Csf1r','P2ry12','Cx3cr1','Mrc1','Skap1'),pt.size=0,ncol = 5)
Idents(MG) <-'seurat_clusters'
MG <- subset(MG, idents = c("0",'1','2','3','5','7')) #### remove macrophage and T cells

DefaultAssay(MG) <- 'integrated'
all.genes <- rownames(MG)
MG<- ScaleData(MG, features = all.genes)
MG<- FindVariableFeatures(object = MG)
MG<- RunPCA(MG, features = VariableFeatures(object = MG))

ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.15)
MG <- RunUMAP(MG, dims = 1: 15)

pdf("LG815_TDI_MG_umap_Condition_res0.15.pdf", width=5.5, height=4)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()


saveRDS(MG, "LG815_TDI_MG_subcluster_res0.15.rds")

# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'LG815_TDI_MG_subcluster_res0.15.rds')

#### cluster marker and DEGs
DefaultAssay(MG) <- 'RNA'
Idents(MG) <-'seurat_clusters'
LG815_TDI_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.25, test.use = "MAST",min.pct = 0.1, only.pos = T)
write.csv(LG815_TDI_MG_markers, "LG815C_MG_markers_res0.15_RNAassay_logFC0.25.csv")

Idents(MG) <-'Condition'
DVKIvsCtrl <- FindMarkers(MG, ident.1 = 'Scn1a: +/-', ident.2 = "Scn1a: +/+", logfc.threshold = 0.1,min.pct = 0.1,
                          test.use = "MAST", assay ='RNA')
write.csv(DVKIvsCtrl, "Dravet_TDIvsDravet_DE_MIN_RNAassay_PV_pct0.1.csv")

write.csv(table(MG$seurat_clusters, MG$Sample_Name), "LG815_TDI_MG_subcluster_cell_counts_res0.15.csv")

                                                                                                