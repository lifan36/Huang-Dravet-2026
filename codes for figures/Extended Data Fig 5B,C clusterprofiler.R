library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(ggplot2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

d = read.csv("Cgashet_DV_MG_Overlap_DEG.csv", header=T)
d <- subset(d, d$color==c('blue'))
d <- d %>% mutate(rank = rank(-(d$avg_log2FC),ties.method = "random")) %>%
  arrange(desc(rank))
## assume 1st column is ID
## 2nd column is FC

## feature 1: numeric vector
geneList = d[,4]

## feature 2: named vector
x=as.character(d[,2])
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
names(geneList) = as.character(eg[,2])
## feature 3: decreasing orde
##geneList = sort(geneList, decreasing = TRUE,ties.method = "random")


### Reactome Analysis
library(ReactomePA)
#data(geneList, package="DOSE")
de <- names(geneList)#[abs(geneList) > 1.5]
head(de)

x <- enrichPathway(gene=de, organism = "mouse", 
                   pvalueCutoff = 0.05,
                   
                   readable=TRUE)
head(x)
selected_pathways <-head(x@result$Description,30)
selected_pathways

### Hallmark Analysis
library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
em <- enricher(de, TERM2GENE=m_t2g,
               pAdjustMethod = "fdr",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.05)
head(em)
em <- setReadable(em, 'org.Hs.eg.db', 'ENTREZID')

####GO analysis
de <- names(geneList)#[abs(geneList) > 1.5]
head(de)
ego <- enrichGO(gene          = de,
                keyType = 'ENTREZID',
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
ego <- simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min)
ego <- setReadable(ego, 'org.Mm.eg.db', 'ENTREZID')
selected_pathways <-head(ego$Description, 30)
selected_pathways

### Extended Data Fig 5B
dotplot(ego, showCategory=8) + ggtitle("dotplot for GSEA")
### Extended Data Fig 5C
cnetplot(ego, showCategory = selected_pathways[c(1,3,4)], foldChange=geneList, circular = TRUE, colorEdge = TRUE,cex_label_gene=0.65,cex_label_category=0.65,cex_category=0.75) 
