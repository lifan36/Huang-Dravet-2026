library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(EnhancedVolcano)
library(SCP)
library(scCustomize)
library(BiocParallel)
library(RColorBrewer)
register(MulticoreParam(workers = 8, progressbar = TRUE))



##### LG815C MG clustering
DefaultAssay(MG) <- 'integrated'
all.genes <- rownames(MG)
MG<- ScaleData(MG, features = all.genes)
MG<- FindVariableFeatures(object = MG)
MG<- RunPCA(MG, features = VariableFeatures(object = MG))

ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:50)
MG <- FindClusters(MG, resolution = 0.4)
MG <- RunUMAP(MG, dims = 1:50)


DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 2)

#####remove very small clusters
MG <- subset(MG, idents = c("0",'1','2','3','4','5'))
MG <- FindNeighbors(MG, dims = 1:50)
MG <- FindClusters(MG, resolution = 0.4)
MG <- RunUMAP(MG, dims = 1:50)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)


MG2[["RNA3"]] <- as(object = MG2[["RNA"]], Class = "Assay")

DefaultAssay(MG2) <- "RNA3"
DefaultAssay(MG) <- 'RNA'
DotPlot(object = MG, features = c('Stat1','Parp14','Rnf213','Ddx60','Trim30a','Cgas',"H2-D1", "H2-K1",'Ifnar1','Cd68'),scale.min = 0) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ RotatedAxis()
Idents(MG) <- 'seurat_clusters'
MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.1, only.pos = T)
write.csv(IFNARKO_OL_markers, "LG815_TDI_MG_markers_res0.15_RNAassay_logFC0.1.csv")


Idents(MG) <-'Condition'
DVKIvsCtrl <- FindMarkers(MG, ident.1 = 'Scn1a: +/-; cGAS: +/-', ident.2 = "Scn1a: +/-; cGAS: +/+", logfc.threshold = 0.15,min.pct = 0.1,
                          test.use = "MAST", assay ='RNA')
write.csv(DVKIvsCtrl, "Dravet_TDIvsDravet_DE_MMG_RNAassay_PV_pct0.1.csv")


saveRDS(MG,"LG815C_MG_final_dim50_res0.4.rds")


#### Fig 3a UMAP with pie chart
CellDimPlot(
  srt = MG, group.by = c("seurat_clusters"),
  reduction = "UMAP", theme_use = "theme_classic"
)

stat.colors <- c("Scn1a: +/+; cGAS: +/+" =  "#488CCA","Scn1a: +/+; cGAS: +/-" =  "#7DD3F6","Scn1a: +/-; cGAS: +/+" ="#EE3425","Scn1a: +/-; cGAS: +/-" = "#F79420")

CellDimPlot(MG, group.by = "seurat_clusters", 
            reduction = "UMAP", stat.by = "Condition",
            theme_use = "theme_classic", legend.position="none",
            stat_palcolor = stat.colors,
            stat_plot_alpha = 3,
            stat_plot_label = FALSE,
            stat_plot_label_size = 1,)

##### Fig 3b
ht <- GroupHeatmap(
  srt = MG,
  #cell_annotation = c("Condition"), cell_annotation_palette = c("Dark2"),
  features = c(
    "Trem2","Ctsl","Cd9", # stage-2 DAM
    'Ifnar2',"Ifngr2", # Interferon
    "H2-K1", "H2-D1","Raet1e" # MHC-II
  ),
  group.by = c( "seurat_clusters","Condition"),
  heatmap_palette = "RdBu",
  show_row_names = FALSE, row_names_side = "left",
  add_dot = TRUE,dot_size=unit(8,"mm"),add_reticle = FALSE,
  add_bg = FALSE,flip = TRUE
  
)
print(ht$plot)
