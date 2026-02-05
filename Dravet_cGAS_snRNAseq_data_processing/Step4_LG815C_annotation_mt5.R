
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/LG815C/integration_mt5")
LG815C_integrated <- readRDS("LG815C_integrated_PCA_0.1.rds")
Idents(LG815C_integrated) <- "seurat_clusters"
#Remove 5 - mixed neuron
LG815C_integrated <- subset(LG815C_integrated, idents=c("5"), invert=T)
Idents(LG815C_integrated) <- "seurat_clusters"
LG815C_integrated <- RenameIdents(LG815C_integrated,
                                `0` = "EN", `1`="OL", `2`="AST", `3`="EN",
                                `4`="EN",`6`="MG", `7`="EN",
                                `8`="EN", `9`="IN", `10`="OPC", `11`="IN",
                                `12`="EN", `13`="VC", `14`="CHOR",`15`="EN",
                                `16`="EN", `17`="EN"
)

pdf("LG815C_integrated_umap_annotation_test.pdf", width=7, height=5)
DimPlot(LG815C_integrated, reduction = 'umap', label = T)
dev.off()

#LG815C_integrated$celltype.orig.ident <- paste(Idents(LG815C_integrated), LG815C_integrated$orig.ident, sep = "_")
LG815C_integrated$celltype <- Idents(LG815C_integrated)
LG815C_integrated$celltype <- factor(LG815C_integrated$celltype, levels = c("EN","IN","OL","AST","MG","OPC","VC","CHOR"))
Idents(LG815C_integrated) <- "celltype"
#markers for annotation
pdf("annotation_1.pdf", width=13, height=4)
DotPlot(LG815C_integrated, features = c("Snap25","Syt1","Slc17a7", "Nrgn","Cux2","Pdzrn3","Mlip","Foxp2","Hs3st4","Grik3",
                                      "Cpne4","Tshz2","Vwc2l","Ntng2","Nr4a2","Col11a1","Gad1","Gad2","Plp1", "Mbp", "Mobp", "Aqp4", "Clu",
                                      "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Vcan", "Pdgfra",
                                      "Cped1","Flt1","Ebf1","Ttr","Htr2c")) + RotatedAxis()
dev.off()


saveRDS(LG815C_integrated, file = "LG815C_integrated_Annotation.rds")

# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(LG815C_integrated$Condition,LG815C_integrated$celltype))
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
a<-as.data.frame(table(LG815C_integrated$Sample_Name,LG815C_integrated$celltype))
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

Idents(LG815C_integrated) <- "celltype"
DefaultAssay(LG815C_integrated) <- 'RNA'
pdf("LG815C_integrated_umap_annotation_noLabel.pdf", width=7, height=5)
DimPlot(LG815C_integrated, reduction = 'umap', label = F)
dev.off()
