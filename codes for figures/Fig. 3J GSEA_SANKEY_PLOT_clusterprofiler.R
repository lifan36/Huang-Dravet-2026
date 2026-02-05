library(ggsankey)
kegg_data <- ego@result
kegg_data2 <- kegg_data %>%
  separate_rows(geneID,sep = "/")
kegg_data2$logFC <- d$avg_log2FC[match(kegg_data2$geneID,d$X)]

sankey <- kegg_data2

sankey <- sankey %>%
  arrange(Description) %>%
  mutate(Description = factor(Description, levels = unique(Description)))
sankey2 <- subset(sankey, sankey$Description==c("gliogenesis")|sankey$Description==c("glial cell differentiation")|sankey$Description==c("regulation of protein catabolic process"))

color_palette <- colorRampPalette(c("#009eff", "#eeeeee", "#f22942"))
# define color of path
ycol2 <- c("#00a08f","#dbc251","#f8b251","#e88e8f","#9a99e1")


sankey2$logFC <- as.numeric(as.character(sankey2$logFC))
mycol3 <- sankey2 %>%
  mutate(color = color_palette(100)[as.numeric(cut(logFC, breaks = seq(-1.25, 1.25, length.out = 101)))]) %>%
  dplyr::select(geneID, color)

# add color to path nodes
pathway_colors <- data.frame(
  path = unique(sankey2$Description),
  color = ycol2[1:length(unique(sankey2$Description))] 
)

# dataframe
sankeyl <- sankey2 %>% make_long(Description,geneID)
# add colors of element and path to sankeyl
sankeyl <- sankeyl %>%
  mutate(color = ifelse(node %in% sankey2$geneID,
                        mycol3$color[match(node, mycol3$geneID)],
                        pathway_colors$color[match(node, pathway_colors$path)]))

p <- ggplot(sankeyl, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                         fill = color, label = node)) +
  scale_fill_identity() + # predefined color
  geom_sankey(flow.alpha = 0.2,
              smooth = 8, 
              width = 0.08) + 
  geom_sankey_text(size = 3.2, color = 'black') +
  theme_void() +
  theme(legend.position = '')
p
# save Sankey plot
ggsave("AST_GO_Sankey.pdf", width = 3, height = 4)



custom_colors <- colorRampPalette(c( "#BBDDF1","#81C9F5", "#009EFF"))(100)
Pathway = c("gliogenesis", "glial cell differentiation", "regulation of protein catabolic process")
GeneRatio = c(13/172, 15/172, 13/172)
counts <- c(13,15,13)
q_value <- c(2.707025e-05,2.698082e-04,2.698082e-04)
y_axis<- c(3,10.5,18.5)

KEGG_point <- data.frame(
  Pathway = c("gliogenesis", "glial cell differentiation", "regulation of protein catabolic process"),
  GeneRatio = round(GeneRatio, 2),
  q_value = round(q_value, 4),
  counts = counts,
  y_axis = y_axis
)

# bubble plot
p1 <- ggplot(KEGG_point, aes(x = y_axis, y = GeneRatio)) +
  coord_flip() +
  geom_point(aes(size = counts, color = q_value), shape = 16) +
  scale_size_continuous(range = c(6, 8)) +
  scale_color_gradientn(colors = custom_colors, name = "q_value") +
  labs(x = NULL, y = bquote("GeneRatio"), size = "Count") +
  scale_y_continuous(limits = c(min(KEGG_point$GeneRatio) - 0.01, max(KEGG_point$GeneRatio) + 0.01)) +
  scale_x_continuous(limits = c(0,20), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
p1
# save bubble plot
ggsave(filename = "bubbles_kegg.pdf", plot = p1, height = 5.5, width = 3.5)


#### add legend 
emptyPlot(1,1, main='Test plot', axes=FALSE)
box()
# legend on outside of plotregion:
gradientLegend(valRange=c(-14,14), pos=.5, side=1)
gradientLegend(valRange=c(-14,14), pos=.5, side=2)
gradientLegend(valRange=c(-14,14), pos=.5, side=3)
gradientLegend(valRange=c(-14,14), pos=.5, side=4)

# legend on inside of plotregion:
gradientLegend(valRange=c(-14,14), pos=.5, side=1, inside=TRUE)
gradientLegend(valRange=c(-14,14), pos=.5, side=2, inside=TRUE)
gradientLegend(valRange=c(-14,14), pos=.5, side=3, inside=TRUE)
gradientLegend(valRange=c(-14,14), pos=.5, side=4, inside=TRUE)

# empty plot:
emptyPlot(1,1, main='Test plot', axes=FALSE)
box()
# number of segments:
gradientLegend(valRange=c(-1.2308088,-0.2542865), color =custom_colors, n.seg=1, side=1)
gradientLegend(valRange=c(-14,14), n.seg=c(-3,5), pos=.5, side=1, 
               inside=TRUE)

# This produces a warning, as there is no space for labels here:    
gradientLegend(valRange=c(-14.235,14.2), pos=.5, 
n.seg = c(-7,0), side=4)# different solutions:
# 1. adjust range (make sure also to adjust the range in the plot, 
#    for example by changing zlim)
emptyPlot(1,1, main='Test plot')
gradientLegend(valRange=c(-14,14), n.seg = c(-7,0), side=4)
# 2. reduce number of decimals:
emptyPlot(1,1, main='Test plot')
gradientLegend(valRange=c(-14.235,14.2), n.seg = c(-7,0), dec=1, side=4)
# 3. change labels to inside plot window:
emptyPlot(1,1, main='Test plot')
gradientLegend(valRange=c(-14.235,14.2), n.seg = c(-7,0), 
               dec=1, side=4, inside=TRUE)
# 4. increase right margin:
oldmar <- par()$mar
par(mar=c(5.1,3.1,4.1,4.1))
emptyPlot(1,1, main='Test plot')
gradientLegend(valRange=c(-14.235,14.2), dec=2, 
               n.seg = c(-7,0), side=4)
par(mar=oldmar) # return old values
# 5. change label position:
emptyPlot(1,1, main='Test plot')
gradientLegend(valRange=c(-14.235,14.2), dec=2, 
               n.seg = c(-7,0), side=4, pos.num=2)
gradientLegend(valRange=c(-14.235,14.2), dec=2, 
               n.seg = c(-7,0), side=4, pos.num=1, pos=.5)
# 6. change legend position and length:
emptyPlot(1,1, main='Test plot')
gradientLegend(valRange=c(-14.235,14.2), dec=2, 
               n.seg = c(-7,0), side=3, length=.5, pos=.75)

# change border color (and font color too!)
gradientLegend(valRange=c(-14,14),pos=.75, length=.5,
               color=alphaPalette('white', f.seq=seq(0,1, by=.1)), 
               border.col=alpha('gray'))

# when defining custom points, it is still important to specify side:

gradientLegend(valRange=c(-14,14), pos=c(.5,.25,.7,-.05), coords=TRUE, 
               border.col='red', side=1)
gradientLegend(valRange=c(-14,14), pos=c(.5,.25,.7,-.05), coords=TRUE, 
               border.col='red', side=2)

