library(ggplot2)

##### Fig 5H dotplot for GSEA

data <- read.csv("Fig 5H LG815_TDI_MG_C5_Reactome_curated.csv")
ggplot() + 
       geom_point(data=data, aes(x=as.numeric(p.value.of.overlap), y=reorder(Gene.Set.Name,-log(FDR.q.value)), colour=as.numeric(p.value.of.overlap),size=X..Genes.in.Overlap..k.),
                                   shape=19, alpha=1) +
       scale_size_continuous(range = c(3, 9))+
       scale_colour_gradient2(low = "#C4B5D7", mid = "white", high = "#CF89BB")+
        xlab("fdr q value") +
       ggtitle("MG Cluster 4 marker")+
       theme_bw()+theme(panel.grid.major.x  = element_blank(),
                                               panel.grid.major.y  = element_blank(),
                                               panel.grid.minor = element_blank())+
       scale_x_continuous()+xlim(5,12.5)
