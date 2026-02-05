
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/LG815_TDI/integration_mt5")
LG815_TDI_integrated <- readRDS("LG815_TDI_integrated_PCA_0.1.rds")
Idents(LG815_TDI_integrated) <- "seurat_clusters"
#Remove 9 - mixed neuron, 13-Granule neuroblasts, dentate gyrus, 19, a very small cluster 509/153k, 0.3%.
LG815_TDI_integrated <- subset(LG815_TDI_integrated, idents=c("9","13","19"), invert=T)
Idents(LG815_TDI_integrated) <- "seurat_clusters"
LG815_TDI_integrated <- RenameIdents(LG815_TDI_integrated,
                                `0` = "EN", `1`="OL", `2`="AST", `3`="EN",
                                `4`="EN", `5`="MG",`6`="OPC", `7`="EN",
                                `8`="IN", `10`="VC", `11`="IN",
                                `12`="CHOR", `14`="EN",`15`="IN",
                                `16`="EN", `17`="EN", `18`="VC"
)

pdf("LG815_TDI_integrated_umap_annotation_test.pdf", width=7, height=5)
DimPlot(LG815_TDI_integrated, reduction = 'umap', label = T)
dev.off()

#LG815_TDI_integrated$celltype.orig.ident <- paste(Idents(LG815_TDI_integrated), LG815_TDI_integrated$orig.ident, sep = "_")
LG815_TDI_integrated$celltype <- Idents(LG815_TDI_integrated)
LG815_TDI_integrated$celltype <- factor(LG815_TDI_integrated$celltype, levels = c("EN","IN","OL","AST","MG","OPC","VC","CHOR"))
Idents(LG815_TDI_integrated) <- "celltype"

saveRDS(LG815_TDI_integrated, file = "LG815_TDI_integrated_Annotation.rds")

# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(LG815_TDI_integrated$Condition,LG815_TDI_integrated$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(), width=4,height=4,units="in")


# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(LG815_TDI_integrated$Sample_Name,LG815_TDI_integrated$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(), width=6,height=4,units="in")

Idents(LG815_TDI_integrated) <- "celltype"
DefaultAssay(LG815_TDI_integrated) <- 'RNA'
pdf("LG815_TDI_integrated_umap_annotation_noLabel.pdf", width=7, height=5)
DimPlot(LG815_TDI_integrated, reduction = 'umap', label = F)
dev.off()



#markers for annotation
pdf("annotation_1.pdf", width=10.5, height=2.7)
DotPlot(data, features = c("Snap25","Slc17a7", "Nrgn", "Gad1", "Gad2", "Adarb1","Epb41l4b","Rspo3","Syt6","Syt9","Lypd6b","Cbln4","Cbln1", "Plp1", "Mbp", "Mobp","Clu", "Plpp3",
                           "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Vcan", "Pdgfra", "Bmp6", "Adam12",
                           "Cped1")) + RotatedAxis()
dev.off()


data <- LG815_TDI_integrated
# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(data$Sample_Name,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/LG815_TDI/integration_mt5",
       width=6,height=4,units="in")

Idents(LG815_TDI_integrated) <- "celltype"
DefaultAssay(LG815_TDI_integrated) <- 'RNA'
pdf("LG815_TDI_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(LG815_TDI_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(LG815_TDI_integrated, idents = "EN")
Cluster_IN <- subset(LG815_TDI_integrated, idents = "IN")
Cluster_MG <- subset(LG815_TDI_integrated, idents = "MG")
Cluster_AST <- subset(LG815_TDI_integrated, idents = "AST")
Cluster_OL <- subset(LG815_TDI_integrated, idents = "OL")
Cluster_OPC <- subset(LG815_TDI_integrated, idents = "OPC")
Cluster_VC <- subset(LG815_TDI_integrated, idents = "VC")


saveRDS(Cluster_EN, file = "LG815_TDI_EN_subset.rds")
saveRDS(Cluster_IN, file = "LG815_TDI_IN_subset.rds")
saveRDS(Cluster_MG, file = "LG815_TDI_MG_subset.rds")
saveRDS(Cluster_AST, file = "LG815_TDI_AST_subset.rds")
saveRDS(Cluster_OL, file = "LG815_TDI_OL_subset.rds")
saveRDS(Cluster_OPC, file = "LG815_TDI_OPC_subset.rds")
saveRDS(Cluster_VC, file = "LG815_TDI_VC_subset.rds")

#######################################################################
setwd("/athena/ganlab/scratch/lif4001/LG815_TDI/integration_mt5/subclustering")
MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.15)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'LG815_TDI_MG_reclusted_res0.15.rds')
pdf("LG815_TDI_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_MG_umap_Condition_res0.15.pdf", width=8, height=2.5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_MG_umap_Sample_res0.15.pdf", width=9.5, height=6)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(MG) <- 'RNA'
LG815_TDI_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_MG_markers, "LG815_TDI_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "LG815_TDI_MG_subcluster_cell_counts_res0.15.csv")
#######################################################################
MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.2)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'LG815_TDI_MG_reclusted_res0.2.rds')
pdf("LG815_TDI_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_MG_umap_Condition_res0.2.pdf", width=8, height=2.5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_MG_umap_Sample_res0.2.pdf", width=9.5, height=6)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(MG) <- 'RNA'
LG815_TDI_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_MG_markers, "LG815_TDI_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "LG815_TDI_MG_subcluster_cell_counts_res0.2.csv")
#######################################################################

MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.3)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'LG815_TDI_MG_reclusted_res0.3.rds')
pdf("LG815_TDI_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_MG_umap_Condition_res0.3.pdf", width=8, height=2.5)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_MG_umap_Sample_res0.3.pdf", width=9.5, height=6)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(MG) <- 'RNA'
LG815_TDI_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_MG_markers, "LG815_TDI_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "LG815_TDI_MG_subcluster_cell_counts_res0.3.csv")
#######################################################################
AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.15)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'LG815_TDI_AST_reclusted_res0.15.rds')
pdf("LG815_TDI_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_AST_umap_Condition.pdf", width=8, height=2.5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_AST_umap_Sample.pdf", width=9.5, height=6)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(AST) <- 'RNA'
LG815_TDI_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_AST_markers, "LG815_TDI_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "LG815_TDI_AST_subcluster_cell_counts.csv")
##################################
AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.2)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'LG815_TDI_AST_reclusted_res0.2.rds')
pdf("LG815_TDI_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_AST_umap_Condition_res0.2.pdf", width=8, height=2.5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_AST_umap_Sample_res0.2.pdf", width=9.5, height=6)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(AST) <- 'RNA'
LG815_TDI_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_AST_markers, "LG815_TDI_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "LG815_TDI_AST_subcluster_cell_counts_res0.2.csv")
#######################################################################

AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.3)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'LG815_TDI_AST_reclusted_res0.3.rds')
pdf("LG815_TDI_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_AST_umap_Condition_res0.3.pdf", width=8, height=2.5)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_AST_umap_Sample_res0.3.pdf", width=9.5, height=6)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(AST) <- 'RNA'
LG815_TDI_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_AST_markers, "LG815_TDI_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "LG815_TDI_AST_subcluster_cell_counts_res0.3.csv")
#######################################################################
#######################################################################
OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.15)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'LG815_TDI_OL_reclusted_res0.15.rds')
pdf("LG815_TDI_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_OL_umap_Condition.pdf", width=8, height=2.5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_OL_umap_Sample.pdf", width=9.5, height=6)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OL) <- 'RNA'
LG815_TDI_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_OL_markers, "LG815_TDI_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "LG815_TDI_OL_subcluster_cell_counts.csv")
#######################################################################

OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.2)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'LG815_TDI_OL_reclusted_res0.2.rds')
pdf("LG815_TDI_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_OL_umap_Condition_res0.2.pdf", width=8, height=2.5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_OL_umap_Sample_res0.2.pdf", width=9.5, height=6)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OL) <- 'RNA'
LG815_TDI_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_OL_markers, "LG815_TDI_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "LG815_TDI_OL_subcluster_cell_counts_res0.2.csv")
#######################################################################

OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.3)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'LG815_TDI_OL_reclusted_res0.3.rds')
pdf("LG815_TDI_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_OL_umap_Condition_res0.3.pdf", width=8, height=2.5)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_OL_umap_Sample_res0.3.pdf", width=9.5, height=6)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OL) <- 'RNA'
LG815_TDI_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_OL_markers, "LG815_TDI_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "LG815_TDI_OL_subcluster_cell_counts_res0.3.csv")
#######################################################################
OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.15)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'LG815_TDI_OPC_reclusted_res0.15.rds')
pdf("LG815_TDI_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_OPC_umap_Condition.pdf", width=8, height=2.5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_OPC_umap_Sample.pdf", width=9.5, height=6)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OPC) <- 'RNA'
LG815_TDI_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_OPC_markers, "LG815_TDI_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "LG815_TDI_OPC_subcluster_cell_counts.csv")
#######################################################################

OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.2)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'LG815_TDI_OPC_reclusted_res0.2.rds')
pdf("LG815_TDI_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_OPC_umap_Condition_res0.2.pdf", width=8, height=2.5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_OPC_umap_Sample_res0.2.pdf", width=9.5, height=6)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OPC) <- 'RNA'
LG815_TDI_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_OPC_markers, "LG815_TDI_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "LG815_TDI_OPC_subcluster_cell_counts_res0.2.csv")
#######################################################################
OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.3)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'LG815_TDI_OPC_reclusted_res0.3.rds')
pdf("LG815_TDI_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_OPC_umap_Condition_res0.3.pdf", width=8, height=2.5)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_OPC_umap_Sample_res0.3.pdf", width=9.5, height=6)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OPC) <- 'RNA'
LG815_TDI_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
write.csv(LG815_TDI_OPC_markers, "LG815_TDI_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "LG815_TDI_OPC_subcluster_cell_counts_res0.3.csv")
#######################################################################
#######################################################################
#######################################################################

EN <- Cluster_EN
DefaultAssay(EN) <- 'integrated'
EN <- ScaleData(EN, verbose = FALSE)
EN <- RunPCA(EN, features = VariableFeatures(object = EN), verbose = FALSE)
ElbowPlot(EN)
EN <- FindNeighbors(EN, dims = 1:15)
EN <- FindClusters(EN, resolution = 0.1)
EN <- RunUMAP(EN, dims = 1: 15)
# rename cluster
n <- dim(table(EN@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
EN@active.ident <- plyr::mapvalues(x = EN@active.ident, from = current.cluster.ids, to = new.cluster.ids)
EN@active.ident <- factor(EN@active.ident, levels=1:n)
saveRDS(EN, file = 'LG815_TDI_EN_reclusted_res0.1.rds')
pdf("LG815_TDI_EN_umap_res0.1.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_EN_umap_Condition_res0.1.pdf", width=8, height=2.5)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_EN_umap_Sample_res0.1.pdf", width=9.5, height=6)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
#DefaultAssay(EN) <- 'RNA'
#LG815_TDI_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
#write.csv(LG815_TDI_EN_markers, "LG815_TDI_EN_markers_res0.1.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "LG815_TDI_EN_subcluster_cell_counts_res0.1.csv")
#######################################################################
#######################################################################

IN <- Cluster_IN
DefaultAssay(IN) <- 'integrated'
IN <- ScaleData(IN, verbose = FALSE)
IN <- RunPCA(IN, features = VariableFeatures(object = IN), verbose = FALSE)
ElbowPlot(IN)
IN <- FindNeighbors(IN, dims = 1:15)
IN <- FindClusters(IN, resolution = 0.1)
IN <- RunUMAP(IN, dims = 1: 15)
# rename cluster
n <- dim(table(IN@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
IN@active.ident <- plyr::mapvalues(x = IN@active.ident, from = current.cluster.ids, to = new.cluster.ids)
IN@active.ident <- factor(IN@active.ident, levels=1:n)
saveRDS(IN, file = 'LG815_TDI_IN_reclusted_res0.1.rds')
pdf("LG815_TDI_IN_umap_res0.1.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("LG815_TDI_IN_umap_Condition_res0.1.pdf", width=8, height=2.5)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()
pdf("LG815_TDI_IN_umap_Sample_res0.1.pdf", width=9.5, height=6)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
#DefaultAssay(IN) <- 'RNA'
#LG815_TDI_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0, only.pos = T)
#write.csv(LG815_TDI_IN_markers, "LG815_TDI_IN_markers_res0.1.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "LG815_TDI_IN_subcluster_cell_counts_res0.1.csv")
#######################################################################
