# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# network analysis & visualization package:
library(igraph)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# load the Zhou et al snRNA-seq dataset
seurat_obj <- readRDS('data/Zhou_control.rds')

ModuleNetworkPlot(DV_MG)

DV_MG <- RunModuleUMAP(
  DV_MG,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(DV_MG)
  
  # plot with ggplot
  ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    geom_point(
      color=umap_df$color, # color each point by WGCNA module
      size=umap_df$kME*2 # size of each point based on intramodular connectivity
    ) +
    umap_theme()


ModuleUMAPPlot(
  DV_MG,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.2, # proportion of edges to sample (20% here)
  label_hubs=5 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
