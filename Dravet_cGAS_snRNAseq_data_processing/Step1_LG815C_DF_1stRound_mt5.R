#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/LG815C/DF_1stRound_mt5")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_106/outs')
sc = autoEstCont(sc)
LG815_C_106.counts = adjustCounts(sc)
LG815_C_106 <- CreateSeuratObject(counts = LG815_C_106.counts, project = "LG815_C_106", min.cells = 3, min.features = 200)
rm(LG815_C_106.counts)
#vizualize QC metrics and filtering====
LG815_C_106[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_106, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_106
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_106_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_106_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_106_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_106_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_121/outs')
sc = autoEstCont(sc)
LG815_C_121.counts = adjustCounts(sc)
LG815_C_121 <- CreateSeuratObject(counts = LG815_C_121.counts, project = "LG815_C_121", min.cells = 3, min.features = 200)
rm(LG815_C_121.counts)
#vizualize QC metrics and filtering====
LG815_C_121[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_121, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_121
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_121_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_121_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_121_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_121_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_122/outs')
sc = autoEstCont(sc)
LG815_C_122.counts = adjustCounts(sc)
LG815_C_122 <- CreateSeuratObject(counts = LG815_C_122.counts, project = "LG815_C_122", min.cells = 3, min.features = 200)
rm(LG815_C_122.counts)
#vizualize QC metrics and filtering====
LG815_C_122[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_122, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_122
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_122_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_122_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_122_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_122_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_132/outs')
sc = autoEstCont(sc)
LG815_C_132.counts = adjustCounts(sc)
LG815_C_132 <- CreateSeuratObject(counts = LG815_C_132.counts, project = "LG815_C_132", min.cells = 3, min.features = 200)
rm(LG815_C_132.counts)
#vizualize QC metrics and filtering====
LG815_C_132[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_132, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_132
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_132_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_132_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_132_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_132_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_139/outs')
sc = autoEstCont(sc)
LG815_C_139.counts = adjustCounts(sc)
LG815_C_139 <- CreateSeuratObject(counts = LG815_C_139.counts, project = "LG815_C_139", min.cells = 3, min.features = 200)
rm(LG815_C_139.counts)
#vizualize QC metrics and filtering====
LG815_C_139[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_139, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_139
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_139_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_139_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_139_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_139_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_153/outs')
sc = autoEstCont(sc)
LG815_C_153.counts = adjustCounts(sc)
LG815_C_153 <- CreateSeuratObject(counts = LG815_C_153.counts, project = "LG815_C_153", min.cells = 3, min.features = 200)
rm(LG815_C_153.counts)
#vizualize QC metrics and filtering====
LG815_C_153[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_153, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_153
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_153_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_153_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_153_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_153_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_159/outs')
sc = autoEstCont(sc)
LG815_C_159.counts = adjustCounts(sc)
LG815_C_159 <- CreateSeuratObject(counts = LG815_C_159.counts, project = "LG815_C_159", min.cells = 3, min.features = 200)
rm(LG815_C_159.counts)
#vizualize QC metrics and filtering====
LG815_C_159[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_159, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_159
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_159_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_159_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_159_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_159_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_164/outs')
sc = autoEstCont(sc)
LG815_C_164.counts = adjustCounts(sc)
LG815_C_164 <- CreateSeuratObject(counts = LG815_C_164.counts, project = "LG815_C_164", min.cells = 3, min.features = 200)
rm(LG815_C_164.counts)
#vizualize QC metrics and filtering====
LG815_C_164[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_164, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_164
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_164_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_164_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_164_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_164_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_172/outs')
sc = autoEstCont(sc)
LG815_C_172.counts = adjustCounts(sc)
LG815_C_172 <- CreateSeuratObject(counts = LG815_C_172.counts, project = "LG815_C_172", min.cells = 3, min.features = 200)
rm(LG815_C_172.counts)
#vizualize QC metrics and filtering====
LG815_C_172[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_172, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_172
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_172_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_172_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_172_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_172_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_29/outs')
sc = autoEstCont(sc)
LG815_C_29.counts = adjustCounts(sc)
LG815_C_29 <- CreateSeuratObject(counts = LG815_C_29.counts, project = "LG815_C_29", min.cells = 3, min.features = 200)
rm(LG815_C_29.counts)
#vizualize QC metrics and filtering====
LG815_C_29[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_29, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_29
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_29_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_29_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_29_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_29_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_30/outs')
sc = autoEstCont(sc)
LG815_C_30.counts = adjustCounts(sc)
LG815_C_30 <- CreateSeuratObject(counts = LG815_C_30.counts, project = "LG815_C_30", min.cells = 3, min.features = 200)
rm(LG815_C_30.counts)
#vizualize QC metrics and filtering====
LG815_C_30[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_30, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_30
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_30_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_30_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_30_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_30_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_31/outs')
sc = autoEstCont(sc)
LG815_C_31.counts = adjustCounts(sc)
LG815_C_31 <- CreateSeuratObject(counts = LG815_C_31.counts, project = "LG815_C_31", min.cells = 3, min.features = 200)
rm(LG815_C_31.counts)
#vizualize QC metrics and filtering====
LG815_C_31[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_31, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_31
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_31_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_31_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_31_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_31_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_33/outs')
sc = autoEstCont(sc)
LG815_C_33.counts = adjustCounts(sc)
LG815_C_33 <- CreateSeuratObject(counts = LG815_C_33.counts, project = "LG815_C_33", min.cells = 3, min.features = 200)
rm(LG815_C_33.counts)
#vizualize QC metrics and filtering====
LG815_C_33[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_33, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_33
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_33_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_33_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_33_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_33_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_35/outs')
sc = autoEstCont(sc)
LG815_C_35.counts = adjustCounts(sc)
LG815_C_35 <- CreateSeuratObject(counts = LG815_C_35.counts, project = "LG815_C_35", min.cells = 3, min.features = 200)
rm(LG815_C_35.counts)
#vizualize QC metrics and filtering====
LG815_C_35[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_35, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_35
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_35_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_35_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_35_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_35_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_64/outs')
sc = autoEstCont(sc)
LG815_C_64.counts = adjustCounts(sc)
LG815_C_64 <- CreateSeuratObject(counts = LG815_C_64.counts, project = "LG815_C_64", min.cells = 3, min.features = 200)
rm(LG815_C_64.counts)
#vizualize QC metrics and filtering====
LG815_C_64[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_64, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_64
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_64_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_64_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_64_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_64_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/LG815C/cellranger/LG815_C_95/outs')
sc = autoEstCont(sc)
LG815_C_95.counts = adjustCounts(sc)
LG815_C_95 <- CreateSeuratObject(counts = LG815_C_95.counts, project = "LG815_C_95", min.cells = 3, min.features = 200)
rm(LG815_C_95.counts)
#vizualize QC metrics and filtering====
LG815_C_95[["percent.mt"]] <- PercentageFeatureSet(object = LG815_C_95, pattern = "^mt-") #recognize mitochondrial transcripts
all <- LG815_C_95
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG815_C_95_FeatureScatter.pdf", width=12, height=4)
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
pdf("LG815_C_95_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"LG815_C_95_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG815_C_95_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
