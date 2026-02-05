seizure <- c("NRP2,MAPK10,TRPC5,KCNJ6,HOMER1,CNTNAP2,KCNK1,ITPR1,PCDH15,MPC1,GRM1,NTRK2,NEUROD2,NF1,PLCL1,BDNF,SYN3,GABRG2,GAP43,NRN1,SPTBN2,KCNA4,SV2C,MALAT1,CACNA1H,GABRG3,PHACTR1,SCG2,CACNA1G,GABRA5,SV2B,CHGB,KCNQ3,KCNC2,SPNS2,HTR4,KCNT2")
seizure <- unlist(strsplit(seizure,","))
foo = function(x){
  paste(toupper(substring(x, 1, 1)),
        tolower(substring(x, 2, nchar(x))),
        sep = "")
}
seizure_gene <- foo(seizure)

###heatmap showing average expression of genes of interest
Idents(AST) <- 'Condition'
cluster.averages <- AggregateExpression(AST, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
seizure_gene <- seizure$X..2000.2021.QIAGEN..All.rights.reserved.[2:nrow(seizure)]
memory_gene <- memory_deficit$ID
dendrite_gene <-dendrites$ID
activity_gene <- Cgas_het_vs_DV_EN_LTD$`Entrez Gene/Gene Symbol - mouse (Entrez Gene)`

cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% activity_gene,]

cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% seizure,]
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% memory_gene,]
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% dendrite_gene,]
cluster.averages_mat <- as.matrix(cluster.averages_mat)
colnames(cluster.averages_mat ) <- c('Scn1a+/+;Cgas+/+','Scn1a+/+;Cgas+/-','Scn1a+/-;Cgas+/+','Scn1a+/-;Cgas+/-')
col <- as.data.frame(colnames(cluster.averages_mat ))
rownames(col) <- colnames(cluster.averages_mat )
colnames(col) <- 'Genotype'

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

pdf("LG815C_EN_heatmap_LTD_genes.pdf", width=2.5, height=4)
pheatmap::pheatmap(cluster.averages_mat,
                   scale="row", clustering_method="ward.D2", color = mypal,annotation_col = col,
                   angle_col = 45, fontsize = 11, main='LTP/LTD genes', legend_labels = 'Average expression',
                   show_colnames = F, annotation_names_col=F,border_color = NA,
                   treeheight_row=F,treeheight_col=F )
dev.off()

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
mat <- cluster.averages_mat[,c(1,4,3)]
#mat <- cluster.averages_mat
mat <- round(scale_rows(mat), 3)


library(circlize)
useT.mat <- t(mat)
dim(useT.mat)
mat_list = useT.mat
dend_list = as.dendrogram(hclust(dist(t(mat_list)))) #changed to calculate for 1 matrix
plot(dend_list,horiz = F)
col_fun = colorRamp2(breaks= c(-1.5, 0, 1.5), 
                     colors = c("#4646CE","white", "#ED2996")) 
circos.par("start.degree" = 90,cell.padding = c(0, 0, 0, 0), gap.degree = 30,gap.after = c(30)) 
circos.initialize("a", xlim =c(0,nrow(mat))) #changed to 1 heatmap setting
circos.track(ylim = c(0, 4), bg.border = NA, track.height = 0.1, 
             panel.fun = function(x, y) {
               for(i in seq_len(ncol(useT.mat))) {
                 circos.text(i-0.5, 0, colnames(useT.mat)[order.dendrogram(dend_list)][i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             cex = 0.5, colnames.side = "inside")                
               }
             })
circos.track(ylim = c(0, 10), bg.border = NA, panel.fun = function(x, y) {
  m = mat_list 
  dend = dend_list
  #changed variable for 1 heatmap setting
  m2 = m[, order.dendrogram(dend)]
  col_mat = col_fun(m2)
  nr = nrow(m2)
  nc = ncol(m2)
  for(i in 1:nr) {
    circos.rect(1:nc - 1, rep(nr - i, nc), 
                1:nc, rep(nr - i + 1, nc), 
                border = col_mat[i, ], col = col_mat[i, ])
  }
  #adding row label
  circos.text(rep(1, 4), 4:1, 
              rownames(useT.mat), 
              facing = "downward", adj = c(1.05,0.5), cex = 0.7) 
  breaks = seq(0, 85, by = 5)
  
})
max_height = attr(dend_list, "height") #changed for 1 dendrogram setting
circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) {
               dend = dend_list
               circos.dendrogram(dend, max_height = max_height)
             })
circos.clear()

library(ComplexHeatmap)
lgd_links = Legend(at = c(-1.5, 0, 1.5), col_fun = col_fun, 
                  title_position = "topleft", title = "Expression", title_gp = gpar(fontsize = 5, fontface = "bold"))

draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))




#trajectoryPlot
barplot_data_overlap_filter <- as.data.frame(as.table(mat))

lineWidth = 0.5
pointSize = 10
Upstream_regulator_trajectory <- ggplot(barplot_data_overlap_filter, aes(x=Var2, y = Freq, group = Var1), fill = Var1)+
  #geom_line(aes(color = PROTEIN), linetype="dashed", size=1) +
  geom_smooth(method = loess,aes(color = Var1), size=0.5) +
  geom_point(aes(color = Var1),size = 1) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size = pointSize),
        plot.background=element_rect(fill="transparent",colour=NA),
        axis.title.x = element_text(color="black", size = pointSize),
        axis.title.y = element_text(color="black", size = pointSize),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 90, vjust=0.1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        axis.ticks = element_line(size = lineWidth, colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = pointSize , colour = "black"),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.background = element_rect(fill="transparent",colour=NA),
        panel.grid.major = element_blank()
        )+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  ggtitle(paste0('Zscore of upstream regulators ',' across Conditions'))+
  scale_y_continuous(expand = c(0, 0.1, .05, 0)) 
# scale_color_manual(values = c(cols)) +
# scale_fill_manual(values = c(cols)) 
Upstream_regulator_trajectory

