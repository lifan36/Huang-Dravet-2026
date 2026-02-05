
# Figure 3C ###################################################################

DVPV_MG = read.csv("DVvsWT_DE_MAST_RNAassay_MG_pct0.15_dim15.csv", header = T)
Cgas_MG = read.csv("Cgashet_DVvsDV_DE_MAST_RNAassay_MG_pct0.15_dim15.csv", header = T)

DVPV_MG_UP = filter(DVPV_MG, p_val_adj < 0.05&(avg_log2FC>=0.1))
DVPV_MG_DN = filter(DVPV_MG, p_val_adj < 0.05&(avg_log2FC<=-0.1))

Cgas_MG_UP = filter(Cgas_MG, p_val_adj < 0.05&(avg_log2FC>=0.1))
Cgas_MG_DN = filter(Cgas_MG, p_val_adj < 0.05&(avg_log2FC<=-0.1))

library(GeneOverlap)
go.obj_UP <- newGeneOverlap(DVPV_MG_UP$X,Cgas_MG_DN$X,genome.size=21988)
go.obj_DN <- newGeneOverlap(DVPV_MG_DN$X,Cgas_MG_UP$X,genome.size=21988)
go.obj_UP <- testGeneOverlap(go.obj_UP)
go.obj_DN <- testGeneOverlap(go.obj_DN)
print(go.obj_UP)
print(go.obj_DN)

#save DEGs that are common
write.csv(go.obj_UP@intersection, 'DV_Cgashet_MG_Overlap_DEG_UP.cvs')
write.csv(go.obj_DN@intersection, 'DV_Cgashet_MG_Overlap_DEG_DN.cvs')

