library(ggplot2)
library(Seurat)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)

# volcano plot of markers of MG4, focusing on DAM genes (Figure 2E) ======
setwd()

dam <-read.csv("DVvsWT_DE_MAST_RNAassay_MG_pct0.15_dim15.csv")
DV_MG <- read.csv("Cgashet_DVvsDV_DE_MAST_RNAassay_MG_pct0.15_dim15.csv")
 

# filter significant DEGs 
dam <- subset(dam, p_val_adj<0.05&(avg_log2FC>=0.1|avg_log2FC<=-0.1))
mg<- subset(DV_MG, p_val_adj<0.05&(avg_log2FC>=0.1|avg_log2FC<=-0.1))
mg$damlogFC <- dam$logFC[match(mg$gene,dam$Gene)]
mg$damlogFC <- dam$avg_log2FC[match(mg$X,dam$X)]

#scatterplot for correlation analysis
mg$color <- "NC"
mg$color[mg$avg_log2FC >0 & mg$damlogFC <0]<-"red"
mg$color[mg$avg_log2FC <0 & mg$damlogFC >0]<-"blue"

write.csv(subset(mg,color=='red'|color=="blue"),'Cgashet_DV_MG_Overlap_DEG.csv')

pdf('DV_Cgashet_MG_cor.pdf', height=4, width = 6)
ggplot(data=subset(mg,color=='red'|color=="blue"), aes(x=avg_log2FC, y=damlogFC)) + 
  geom_point(aes(color= color), size=2,alpha = 1/5)+
  scale_color_manual(values = c( "red" = "salmon", "blue" = "dodgerblue"))+
  geom_text_repel(subset(mg, X==c('H2-D1')|X==c('Ctsl')|X==c('Trem2')|X==c('Il17ra')|X==c('Ifnar2')|X==c('Tbk1')|X==c('Tmem173')|X==c('Aqp4')|X==c('Ddr1')|X==c('Tlr3')|X==c('H2-K1')|X==c('Prdx6')), mapping = aes(label = X), 
                  max.overlaps = nrow(mg),size = 5,segment.color = "Black",min.segment.length = 0)+
  #geom_text_repel(subset(mg, X==c('Cd9')|X==c('Gfap')|X==c('C4b')|X==c('Id3')|X==c('Apoe')|X==c('Cst3')|X==c('Cd81')|X==c('Aqp4')|X==c('Ddr1')|X==c('Tlr3')|X==c('H2-K1')|X==c('Prdx6')), mapping = aes(label = X), 
                  #max.overlaps = nrow(mg),size = 5,segment.color = "Black",min.segment.length = 0)+
  geom_text_repel(subset(mg, X==c('Bdnf')|X==c('Stxbp3')|X==c('Syn3')|X==c('Syndig1')|X==c('Unc5d')|X==c('Cacna1g')|X==c('Satb2')|X==c('Ryr3')|X==c('Pde10a')|X==c('Camk4')|X==c('Robo1')|X==c('Pde7b')|X==c('Homer1')), mapping = aes(label = X), 
  max.overlaps = nrow(mg),size = 5,segment.color = "Black",min.segment.length = 0)+
  
  geom_smooth(method = "lm", se = T)+
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  geom_hline(yintercept = c(0), linetype = "dotted") +
  geom_vline(xintercept = c(0), linetype = "dotted") +
  ylab("logFC[Scn1a+/- vs Scn1a+/+]") + xlab("logFC[Scn1a+/-;Cgas+/- vs Scn1a+/-]") 
dev.off()
cor(mg$avg_log2FC, mg$damlogFC, use = "complete.obs") #0.7908248
cor.test(mg$avg_log2FC, mg$damlogFC, use = "complete.obs", method = "pearson")




