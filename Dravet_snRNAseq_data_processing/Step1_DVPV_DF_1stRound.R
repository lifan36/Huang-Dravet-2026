#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/DVPV/DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_CTX_311/outs')
sc = autoEstCont(sc)
DVPV_CTX_311.counts = adjustCounts(sc)
DVPV_CTX_311 <- CreateSeuratObject(counts = DVPV_CTX_311.counts, project = "DVPV_CTX_311", min.cells = 3, min.features = 200)
rm(DVPV_CTX_311.counts)
#vizualize QC metrics and filtering====
DVPV_CTX_311[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_CTX_311, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_CTX_311
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_CTX_311_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_CTX_311_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_CTX_311_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_CTX_311_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_CTX_324/outs')
sc = autoEstCont(sc)
DVPV_CTX_324.counts = adjustCounts(sc)
DVPV_CTX_324 <- CreateSeuratObject(counts = DVPV_CTX_324.counts, project = "DVPV_CTX_324", min.cells = 3, min.features = 200)
rm(DVPV_CTX_324.counts)
#vizualize QC metrics and filtering====
DVPV_CTX_324[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_CTX_324, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_CTX_324
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_CTX_324_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_CTX_324_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_CTX_324_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_CTX_324_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_CTX_327/outs')
sc = autoEstCont(sc)
DVPV_CTX_327.counts = adjustCounts(sc)
DVPV_CTX_327 <- CreateSeuratObject(counts = DVPV_CTX_327.counts, project = "DVPV_CTX_327", min.cells = 3, min.features = 200)
rm(DVPV_CTX_327.counts)
#vizualize QC metrics and filtering====
DVPV_CTX_327[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_CTX_327, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_CTX_327
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_CTX_327_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_CTX_327_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_CTX_327_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_CTX_327_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_CTX_329/outs')
sc = autoEstCont(sc)
DVPV_CTX_329.counts = adjustCounts(sc)
DVPV_CTX_329 <- CreateSeuratObject(counts = DVPV_CTX_329.counts, project = "DVPV_CTX_329", min.cells = 3, min.features = 200)
rm(DVPV_CTX_329.counts)
#vizualize QC metrics and filtering====
DVPV_CTX_329[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_CTX_329, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_CTX_329
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_CTX_329_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_CTX_329_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_CTX_329_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_CTX_329_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_CTX_330/outs')
sc = autoEstCont(sc)
DVPV_CTX_330.counts = adjustCounts(sc)
DVPV_CTX_330 <- CreateSeuratObject(counts = DVPV_CTX_330.counts, project = "DVPV_CTX_330", min.cells = 3, min.features = 200)
rm(DVPV_CTX_330.counts)
#vizualize QC metrics and filtering====
DVPV_CTX_330[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_CTX_330, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_CTX_330
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_CTX_330_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_CTX_330_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_CTX_330_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_CTX_330_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_CTX_333/outs')
sc = autoEstCont(sc)
DVPV_CTX_333.counts = adjustCounts(sc)
DVPV_CTX_333 <- CreateSeuratObject(counts = DVPV_CTX_333.counts, project = "DVPV_CTX_333", min.cells = 3, min.features = 200)
rm(DVPV_CTX_333.counts)
#vizualize QC metrics and filtering====
DVPV_CTX_333[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_CTX_333, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_CTX_333
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_CTX_333_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_CTX_333_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_CTX_333_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_CTX_333_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_Hipp_311/outs')
sc = autoEstCont(sc)
DVPV_Hipp_311.counts = adjustCounts(sc)
DVPV_Hipp_311 <- CreateSeuratObject(counts = DVPV_Hipp_311.counts, project = "DVPV_Hipp_311", min.cells = 3, min.features = 200)
rm(DVPV_Hipp_311.counts)
#vizualize QC metrics and filtering====
DVPV_Hipp_311[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_Hipp_311, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_Hipp_311
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_Hipp_311_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_Hipp_311_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_Hipp_311_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_Hipp_311_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_Hipp_324/outs')
sc = autoEstCont(sc)
DVPV_Hipp_324.counts = adjustCounts(sc)
DVPV_Hipp_324 <- CreateSeuratObject(counts = DVPV_Hipp_324.counts, project = "DVPV_Hipp_324", min.cells = 3, min.features = 200)
rm(DVPV_Hipp_324.counts)
#vizualize QC metrics and filtering====
DVPV_Hipp_324[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_Hipp_324, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_Hipp_324
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_Hipp_324_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_Hipp_324_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_Hipp_324_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_Hipp_324_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_Hipp_327/outs')
sc = autoEstCont(sc)
DVPV_Hipp_327.counts = adjustCounts(sc)
DVPV_Hipp_327 <- CreateSeuratObject(counts = DVPV_Hipp_327.counts, project = "DVPV_Hipp_327", min.cells = 3, min.features = 200)
rm(DVPV_Hipp_327.counts)
#vizualize QC metrics and filtering====
DVPV_Hipp_327[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_Hipp_327, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_Hipp_327
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_Hipp_327_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_Hipp_327_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_Hipp_327_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_Hipp_327_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_Hipp_329/outs')
sc = autoEstCont(sc)
DVPV_Hipp_329.counts = adjustCounts(sc)
DVPV_Hipp_329 <- CreateSeuratObject(counts = DVPV_Hipp_329.counts, project = "DVPV_Hipp_329", min.cells = 3, min.features = 200)
rm(DVPV_Hipp_329.counts)
#vizualize QC metrics and filtering====
DVPV_Hipp_329[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_Hipp_329, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_Hipp_329
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_Hipp_329_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_Hipp_329_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_Hipp_329_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_Hipp_329_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_Hipp_330/outs')
sc = autoEstCont(sc)
DVPV_Hipp_330.counts = adjustCounts(sc)
DVPV_Hipp_330 <- CreateSeuratObject(counts = DVPV_Hipp_330.counts, project = "DVPV_Hipp_330", min.cells = 3, min.features = 200)
rm(DVPV_Hipp_330.counts)
#vizualize QC metrics and filtering====
DVPV_Hipp_330[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_Hipp_330, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_Hipp_330
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_Hipp_330_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_Hipp_330_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_Hipp_330_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_Hipp_330_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/DVPV/cellranger/DVPV_Hipp_333/outs')
sc = autoEstCont(sc)
DVPV_Hipp_333.counts = adjustCounts(sc)
DVPV_Hipp_333 <- CreateSeuratObject(counts = DVPV_Hipp_333.counts, project = "DVPV_Hipp_333", min.cells = 3, min.features = 200)
rm(DVPV_Hipp_333.counts)
#vizualize QC metrics and filtering====
DVPV_Hipp_333[["percent.mt"]] <- PercentageFeatureSet(object = DVPV_Hipp_333, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DVPV_Hipp_333
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DVPV_Hipp_333_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DVPV_Hipp_333_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"DVPV_Hipp_333_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DVPV_Hipp_333_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
