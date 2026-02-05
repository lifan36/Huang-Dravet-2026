# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')
dbs <- c('Reactome_2022')
dbs <- c('MSigDB_Hallmark_2020')
# perform enrichment tests
DV_MG <- RunEnrichr(
  DV_MG,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(DV_MG)

# enrichr dotplot
pdf('Gillian_E3_MG_WGCNA_enrichplot.pdf', height = 5.5, width=6.5)
EnrichrDotPlot(
  DV_MG,
  mods ="all", # use all modules (this is the default behavior)
  database = "MSigDB_Hallmark_2020", # this has to be one of the lists we used above!!!
  n_terms=5 # number of terms for each module
)
dev.off()
