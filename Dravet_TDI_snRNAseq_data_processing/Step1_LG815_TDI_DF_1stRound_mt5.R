#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/LG815_TDI/DF_1stRound_mt5")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_160/outs')
sc = autoEstCont(sc)
LG815_Ctrl_160.counts = adjustCounts(sc)
LG815_Ctrl_160 <- CreateSeuratObject(counts = LG815_Ctrl_160.counts, project = "LG815_Ctrl_160", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_160.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_160[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_160, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_160
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_160_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_160_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_160_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_160_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_168/outs')
sc = autoEstCont(sc)
LG815_Ctrl_168.counts = adjustCounts(sc)
LG815_Ctrl_168 <- CreateSeuratObject(counts = LG815_Ctrl_168.counts, project = "LG815_Ctrl_168", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_168.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_168[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_168, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_168
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_168_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_168_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_168_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_168_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_180/outs')
sc = autoEstCont(sc)
LG815_Ctrl_180.counts = adjustCounts(sc)
LG815_Ctrl_180 <- CreateSeuratObject(counts = LG815_Ctrl_180.counts, project = "LG815_Ctrl_180", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_180.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_180[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_180, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_180
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_180_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_180_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_180_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_180_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_216/outs')
sc = autoEstCont(sc)
LG815_Ctrl_216.counts = adjustCounts(sc)
LG815_Ctrl_216 <- CreateSeuratObject(counts = LG815_Ctrl_216.counts, project = "LG815_Ctrl_216", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_216.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_216[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_216, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_216
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_216_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_216_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_216_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_216_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_162/outs')
sc = autoEstCont(sc)
LG815_Ctrl_162.counts = adjustCounts(sc)
LG815_Ctrl_162 <- CreateSeuratObject(counts = LG815_Ctrl_162.counts, project = "LG815_Ctrl_162", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_162.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_162[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_162, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_162
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_162_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_162_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_162_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_162_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_169/outs')
sc = autoEstCont(sc)
LG815_Ctrl_169.counts = adjustCounts(sc)
LG815_Ctrl_169 <- CreateSeuratObject(counts = LG815_Ctrl_169.counts, project = "LG815_Ctrl_169", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_169.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_169[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_169, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_169
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_169_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_169_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_169_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_169_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_174/outs')
sc = autoEstCont(sc)
LG815_Ctrl_174.counts = adjustCounts(sc)
LG815_Ctrl_174 <- CreateSeuratObject(counts = LG815_Ctrl_174.counts, project = "LG815_Ctrl_174", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_174.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_174[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_174, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_174
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_174_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_174_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_174_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_174_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_Ctrl_195/outs')
sc = autoEstCont(sc)
LG815_Ctrl_195.counts = adjustCounts(sc)
LG815_Ctrl_195 <- CreateSeuratObject(counts = LG815_Ctrl_195.counts, project = "LG815_Ctrl_195", min.cells = 3, min.features = 200)
rm(LG815_Ctrl_195.counts)
#vizualize QC metrics and filtering====
LG815_Ctrl_195[["percent.mt"]] <- PercentageFeatureSet(object = LG815_Ctrl_195, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_Ctrl_195
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_Ctrl_195_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_Ctrl_195_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_Ctrl_195_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_Ctrl_195_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_141/outs')
sc = autoEstCont(sc)
LG815_TDI_141.counts = adjustCounts(sc)
LG815_TDI_141 <- CreateSeuratObject(counts = LG815_TDI_141.counts, project = "LG815_TDI_141", min.cells = 3, min.features = 200)
rm(LG815_TDI_141.counts)
#vizualize QC metrics and filtering====
LG815_TDI_141[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_141, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_141
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_141_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_141_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_141_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_141_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_163/outs')
sc = autoEstCont(sc)
LG815_TDI_163.counts = adjustCounts(sc)
LG815_TDI_163 <- CreateSeuratObject(counts = LG815_TDI_163.counts, project = "LG815_TDI_163", min.cells = 3, min.features = 200)
rm(LG815_TDI_163.counts)
#vizualize QC metrics and filtering====
LG815_TDI_163[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_163, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_163
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_163_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_163_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_163_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_163_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_204/outs')
sc = autoEstCont(sc)
LG815_TDI_204.counts = adjustCounts(sc)
LG815_TDI_204 <- CreateSeuratObject(counts = LG815_TDI_204.counts, project = "LG815_TDI_204", min.cells = 3, min.features = 200)
rm(LG815_TDI_204.counts)
#vizualize QC metrics and filtering====
LG815_TDI_204[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_204, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_204
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_204_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_204_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_204_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_204_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_205/outs')
sc = autoEstCont(sc)
LG815_TDI_205.counts = adjustCounts(sc)
LG815_TDI_205 <- CreateSeuratObject(counts = LG815_TDI_205.counts, project = "LG815_TDI_205", min.cells = 3, min.features = 200)
rm(LG815_TDI_205.counts)
#vizualize QC metrics and filtering====
LG815_TDI_205[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_205, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_205
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_205_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_205_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_205_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_205_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_149/outs')
sc = autoEstCont(sc)
LG815_TDI_149.counts = adjustCounts(sc)
LG815_TDI_149 <- CreateSeuratObject(counts = LG815_TDI_149.counts, project = "LG815_TDI_149", min.cells = 3, min.features = 200)
rm(LG815_TDI_149.counts)
#vizualize QC metrics and filtering====
LG815_TDI_149[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_149, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_149
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_149_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_149_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_149_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_149_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_152/outs')
sc = autoEstCont(sc)
LG815_TDI_152.counts = adjustCounts(sc)
LG815_TDI_152 <- CreateSeuratObject(counts = LG815_TDI_152.counts, project = "LG815_TDI_152", min.cells = 3, min.features = 200)
rm(LG815_TDI_152.counts)
#vizualize QC metrics and filtering====
LG815_TDI_152[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_152, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_152
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_152_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_152_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_152_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_152_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_160/outs')
sc = autoEstCont(sc)
LG815_TDI_160.counts = adjustCounts(sc)
LG815_TDI_160 <- CreateSeuratObject(counts = LG815_TDI_160.counts, project = "LG815_TDI_160", min.cells = 3, min.features = 200)
rm(LG815_TDI_160.counts)
#vizualize QC metrics and filtering====
LG815_TDI_160[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_160, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_160
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_160_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_160_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_160_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_160_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815_TDI/cellranger/LG815_TDI_199/outs')
sc = autoEstCont(sc)
LG815_TDI_199.counts = adjustCounts(sc)
LG815_TDI_199 <- CreateSeuratObject(counts = LG815_TDI_199.counts, project = "LG815_TDI_199", min.cells = 3, min.features = 200)
rm(LG815_TDI_199.counts)
#vizualize QC metrics and filtering====
LG815_TDI_199[["percent.mt"]] <- PercentageFeatureSet(object = LG815_TDI_199, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_TDI_199
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_TDI_199_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4500 & nCount_RNA < 15000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG815_TDI_199_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_TDI_199_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_TDI_199_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
