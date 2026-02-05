#single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading, for non-mac
# enableWGCNAThreads(nThreads = 8)
# for mac
allowWGCNAThreads(nThreads = 8)

DefaultAssay(DV_MG) <- 'RNA'
Idents(DV_MG) <- 'seurat_clusters'

# load the Zhou et al snRNA-seq dataset
#DV_MG <- readRDS('Zhou_2020.rds')
DV_MG <- MG2

DV_MG <- SetupForWGCNA(
  DV_MG,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.01, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "test" # the name of the hdWGCNA experiment
)


# construct metacells  in each group
DV_MG <- MetacellsByGroups(
  DV_MG,
  group.by = c("seurat_clusters"), # specify the columns in DV_MG@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
DV_MG <- NormalizeMetacells(DV_MG)

DV_MG <- SetDatExpr(
  DV_MG,
  group_name = c('0','1','2','3'), # the name of the group of interest in the group.by column
  group.by='seurat_clusters', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

DV_MG <- TestSoftPowers(
  DV_MG,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(DV_MG)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(DV_MG)
head(power_table)

# construct co-expression network:
DV_MG <- ConstructNetwork(
  DV_MG, soft_power=3,
  setDatExpr=FALSE,
  tom_name = 'MG', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE
)

PlotDendrogram(DV_MG, main='MG hdWGCNA Dendrogram')

TOM <- GetTOM(DV_MG)
# need to run ScaleData first or else harmony throws an error:
DV_MG <- ScaleData(DV_MG, features=VariableFeatures(DV_MG))

# compute all MEs in the full single-cell dataset
DV_MG <- ModuleEigengenes(
  DV_MG,
  group.by.vars='seurat_clusters'
)

# harmonized module eigengenes:
hMEs <- GetMEs(DV_MG)

# module eigengenes:
MEs <- GetMEs(DV_MG, harmonized=FALSE)
write.csv(MEs, 'Module_Eigengen_cell.csv')
# compute eigengene-based connectivity (kME):
DV_MG <- ModuleConnectivity(
  DV_MG,
  group.by = 'seurat_clusters', group_name = c('0','1','2','3'), sparse=FALSE
)

# rename the modules
#DV_MG <- ResetModuleNames(
#DV_MG,
#  new_name = "INH-M"
#)

# plot genes ranked by kME for each module
p <- PlotKMEs(DV_MG, ncol=5)

p

# get the module assignment table:
modules <- GetModules(DV_MG)
write.csv(modules, 'DVPV_MG_wgcna_module_marker.csv')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(DV_MG, n_hubs = 10)

head(hub_df)

saveRDS(DV_MG, file='DVPV_hdWGCNA_MG_allclusters.rds')

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method

DV_MG <- ModuleExprScore(
  DV_MG,
  n_genes = 25,
  method='Seurat'
)


# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  DV_MG,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

# stitch together with patchwork
pdf('DVPV_hdWGCNA_Module_Feature.pdf', height=4, width=8)
wrap_plots(plot_list, ncol=6)
dev.off()

# plot module correlagram
ModuleCorrelogram(DV_MG)

# get hMEs from seurat object
MEs <- GetMEs(DV_MG, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
#mods <- colnames(MEs); mods <- mods[mods == 'yellow'|mods == 'turquoise'|mods == 'blue']
# add hMEs to Seurat meta-data:
DV_MG@meta.data <- cbind(DV_MG@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(DV_MG, features=mods, group.by = 'seurat_clusters')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
pdf('DVPV_MG_WGCNA_module_cluster_dot.pdf', height = 3, width = 5)
p
dev.off()

#### Differential Module Expression Analysis
group1 <- DV_MG@meta.data %>% subset(seurat_clusters == c('0','1','2','3','4') & Condition == 'E3_P301S_R47H') %>% rownames
group2 <- DV_MG@meta.data %>% subset(seurat_clusters == c('0','1','2','3','4') & Condition == 'E3_P301S_WT') %>% rownames

head(group1)

DMEs <- FindDMEs(
  DV_MG,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='test'
)

head(DMEs)
write.csv(DMEs, 'DVPV_MG_WGCNA_bycondition.csv')
pdf('DVPV_MG_WGCNA_condition_lolipop_all.pdf', height=3, width = 3.5)
PlotDMEsLollipop(
  DV_MG, 
  DMEs, 
  wgcna_name='test', 
  pvalue = "p_val_adj"
)
dev.off()
