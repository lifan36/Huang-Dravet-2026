library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MIN)
library(EnhancedVolcano)

DV <- readRDS("DVPV_integrated_PCA_0.1_both.rds")
DefaultAssay(DV) <- 'integrated'
DV <- ScaleData(DV, verbose = FALSE)
DV <- RunPCA(DV, features = VariableFeatures(object = DV), verbose = FALSE)
ElbowPlot(DV)
DV <- FindNeighbors(DV, dims = 1:15)
DV <- FindClusters(DV, resolution = 0.25)
DV <- RunUMAP(DV, dims = 1: 15)
DimPlot(DV, reduction = 'umap', label = T)

sig <-c("Flt1",'Vwf','Pecam1','Cdh5','Cd96','Skap1','Cd3g','Cd3e','Cd3d',
        'Plxdc2','Cybb','Apbb1ip','Mobp','Plp1','Mbp','Sox10','Pcdh15','Lhfpl3',
        'Vcan','Pdgfra','Cspg4','Dock8','Mrc1','Sfmbt2','Cx3cr1','C1qc','Erbb4','Nxph1',
        'Gad2','Gad1','Syn3','Rbfox3','Camk2a',"Slc17a7", "Nrgn",'Aldh1l1','Aqp4','Gja1','Fgfr3','Gfap','Ttr','Clic6')
markers.to.plot <- as.matrix(sig)
DefaultAssay(DV) <- 'RNA'
DotPlot(object = DV, features = rev(x = markers.to.plot)) + RotatedAxis()

### both 2,3,4,7,9,12,13,15,18,19 EN, 6,8,11,14,22 IN, 0 Oligo, 
#20 ependymal, 23 endothelial, 10 OPC, 1,25 AST, 5,26 MG, 17 ambiguous  