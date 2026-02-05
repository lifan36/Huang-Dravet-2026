
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/DVPV/integration_both")
DVPV_integrated <- readRDS("DVPV_integrated_PCA_0.1.rds")
Idents(DVPV_integrated) <- "seurat_clusters"
#Remove cluster 14 - Choroid plexus epithelial cells, cluster 10-Reln+ CR cells
DVPV_integrated <- subset(DVPV_integrated, idents=c("13","18"), invert=T)
Idents(DVPV_integrated) <- "seurat_clusters"
DVPV_integrated <- RenameIdents(DVPV_integrated,
                                 `0` = "OL_Mbp", `1`="CTX_Cux2", `2`="AST_Atp1a2", `3`="DG_Prox1",
                                 `4`="IN_Gad1", `5`="CTX_Foxp2", `6`="CA1_CA2_Cpne7", `7`="MG_Cx3cr1",
                                 `8`="IN_Gad1", `9`="OPC_Vcan", `10`="CTX_Il1rapl2", `11`="IN_Gad1",
                                 `12`="Reln", `14`="VC_Cped1", `15`="CA3_Cpne4", `16`="CHOR_Ttr", `17`="VC_Flt1"
)

pdf("DVPV_integrated_umap_annotation_test.pdf", width=7, height=5)
DimPlot(DVPV_integrated, reduction = 'umap', label = T)
dev.off()

#DVPV_integrated$celltype.orig.ident <- paste(Idents(DVPV_integrated), DVPV_integrated$orig.ident, sep = "_")
DVPV_integrated$celltype <- Idents(DVPV_integrated)
DVPV_integrated$celltype <- factor(DVPV_integrated$celltype, levels = c("CTX_Cux2","CTX_Foxp2","CTX_Il1rapl2","DG_Prox1","CA1_CA2_Cpne7","CA3_Cpne4",
                                                                        "IN_Gad1","OL_Mbp","AST_Atp1a2","MG_Cx3cr1","OPC_Vcan","VC_Cped1","VC_Flt1","Reln","CHOR_Ttr"))
Idents(DVPV_integrated) <- "celltype"
#markers for annotation
pdf("annotation_1.pdf", width=13, height=4)
DotPlot(DVPV_integrated, features = c("Snap25","Syt1","Slc17a7", "Nrgn","Cux2","Pdzrn3","Mlip","Foxp2","Hs3st4","Grik3","Il1rapl2","Satb2","Hs3st2","Prox1","Htr4","Glis3","Cpne7","Prkcg","Cntnap5c",
                                      "Cpne4","Tshz2","Vwc2l","Gad1","Gad2","Plp1", "Mbp", "Mobp", "Aqp4", "Clu",
                                      "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Vcan", "Pdgfra", "Bmp6", "Adam12",
                                      "Cped1","Flt1","Ebf1","Reln","Unc13c","Ttr","Htr2c")) + RotatedAxis()
dev.off()


saveRDS(DVPV_integrated, file = "DVPV_integrated_Annotation.rds")

# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(DVPV_integrated$Condition,DVPV_integrated$celltype))
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
a<-as.data.frame(table(DVPV_integrated$Sample_Name,DVPV_integrated$celltype))
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

Idents(DVPV_integrated) <- "celltype"
DefaultAssay(DVPV_integrated) <- 'RNA'
pdf("DVPV_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(DVPV_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(DVPV_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(DVPV_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(DVPV_integrated, idents = "MG_Cx3cr1")
Cluster_AST <- subset(DVPV_integrated, idents = "astrocytes")
Cluster_OL <- subset(DVPV_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(DVPV_integrated, idents = "OPCs")
Cluster_VC <- subset(DVPV_integrated, idents = "vascular cells")


saveRDS(Cluster_EN, file = "DVPV_EN_subset.rds")
saveRDS(Cluster_IN, file = "DVPV_IN_subset.rds")
saveRDS(Cluster_MG, file = "DVPV_MG_subset.rds")
saveRDS(Cluster_AST, file = "DVPV_AST_subset.rds")
saveRDS(Cluster_OL, file = "DVPV_OL_subset.rds")
saveRDS(Cluster_OPC, file = "DVPV_OPC_subset.rds")
saveRDS(Cluster_VC, file = "DVPV_VC_subset.rds")


