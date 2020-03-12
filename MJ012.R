library(dplyr)
library(Seurat)

genelist = c("CD79A", "CD74", "CD3E", "CD3D", "CD3G", "CD4", "CD8A", "TRDC", "CD69", "ITGAE", "ITGA1", "CXCR6", "RUNX3", "PRDM1", "ZNF683", "CCR7", "KLF2", "S1PR1", "SELL", "MKI67", "CD28","TMIGD2", "TIGIT","TCF7", "RORA", "RORC", "AHR", "IL23R", "IL23A", "IL22RA1", "IL22", "IL17RA","IL17A","CCL20", "TNF", "IFNG", "IFNGR1","TBX21", "EOMES", "IL4", "IL5", "IL13", "GATA3", "BCL6", "CXCR5", "PDCD1", "ICOS", "IL2RA", "CXCL13", "IL21", "IL6", "STAT3", "IRF4", "IL21R", "IL6R", "BATF", "STAT5A", "MAF", "GZMB", "GZMA","PRF1","GNLY", "KLRD1", "KLRC1", "KLRK1","FOXP3", "CTLA4","IL10","TGFB1", "IL2","NR4A1", "SATB1", "EGR1", "IL4I1", "NFATC1","JUN", "FOS","NFKBIA", "NFKBID","CD83", "PTPN22", "CXCR4", "BTG2", "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", "HLA-E", "HLA-F", "LYN", "CCL4")

mj012_data <- Read10X(data.dir = "~/Documents/workstation/data/jian/MJ012")
mj012_so <- CreateSeuratObject(counts = mj012_data, project = "mj012")
mj012_so[["percent.mt"]] <- PercentageFeatureSet(mj012_so, pattern = "^MT-")

mj012_Q1 = quantile(mj012_so$nFeature_RNA)[2]
mj012_Q3 = quantile(mj012_so$nFeature_RNA)[4]
mj012_interQ = (mj012_Q3-mj012_Q1)*1.5
mj012_upperb = mj012_Q3 + mj012_interQ
mj012_lowerb = mj012_Q1 - mj012_interQ

VlnPlot(mj012_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

mj012_so <- subset(mj012_so, subset = nFeature_RNA > mj012_lowerb  & nFeature_RNA < mj012_upperb  & percent.mt <15)
mj012_so <- NormalizeData(mj012_so, normalization.method = "LogNormalize", scale.factor = 10000)
mj012_so <- FindVariableFeatures(mj012_so, selection.method = "vst", nfeatures = 2000)
mj012_all.genes <- rownames(mj012_so)
mj012_so <- ScaleData(mj012_so, features = mj012_all.genes)

mj012_so <- RunPCA(mj012_so, features = VariableFeatures(object = mj012_so))
#mj012_so <- JackStraw(mj012_so, num.replicate = 100)
#mj012_so <- ScoreJackStraw(mj012_so, dims = 1:20)

#ElbowPlot(mj012_so)

mj012_so <- FindNeighbors(mj012_so, dims = 1:15)
mj012_so <- FindClusters(mj012_so, resolution = 0.4)

mj012_so <- RunUMAP(mj012_so, dims = 1:15)


mj012_so <- mj012_so
#"MJ012" <-'mj012'

colors <- c('#E91E63','#3F51B5','#00BCD4','#8BC34A','#FFC107','#795548',
            '#9C27B0','#2196F3','#010688','#CDDC39','#FF9800','#9E9E9E',
            '#673AB7','#03A9F4','#4CAF50','#FFEB3B','#FF5722','#607D8B',
            '#F44336','#000000',
            '#F48FB1','#9FA8DA','#80DEEA','#C5E1A5','#FFE082','#BCAAA4',
            '#CE93D8','#90CAF9','#80CBC4','#E6EE9C','#FFCC80','#EEEEEE'
            )


pdf(paste("~/Desktop/data/mj012_output/","MJ012","_UMAP.pdf",sep=''), width = 8, height = 4.5)
DimPlot(mj012_so, cols=colors, label =TRUE) + NoLegend()
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_UMAP_right.jpeg", width = 1600, height = 900, res=300)
DimPlot(mj012_so, cols=colors)
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_UMAP_free.jpeg", width = 1600, height = 900, res=300)
DimPlot(mj012_so, cols=colors)+ NoLegend()
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_UMAP_on.jpeg", width = 1600, height = 900, res=300)
DimPlot(mj012_so, cols=colors, label =TRUE)+ NoLegend()
dev.off()



jpeg(paste("~/Desktop/data/mj012_output/","MJ012","_HEAT.jpeg",sep=''), width = 1200, height = 1200)
DoHeatmap(mj012_so , features = genelist,group.colors =colors)
dev.off()

mj012_so.markers <- FindAllMarkers(mj012_so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

drop.index <- grep(pattern = "^TR[ABDVG]|^IG[HLK][A-Z][1-9]{0,1}", x = mj012_so.markers$gene, value = FALSE)
mj012_so.markers <- mj012_so.markers[-drop.index, ]

mj012_top30 <- mj012_so.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_logFC)


jpeg(paste("~/Desktop/data/mj012_output/","MJ012","_HEAT_makers.jpeg",sep=''), width = 1200, height = length(table(mj012_so$seurat_clusters))*120)
DoHeatmap(mj012_so , features = c('CD3G','CD4','CD8A','TRDC', mj012_top30$gene), group.colors =colors)
dev.off()


mj012_so_cluster.averages <- AverageExpression(mj012_so, return.seurat=TRUE) 
jpeg(paste("~/Desktop/data/mj012_output/","MJ012","_HEAT_cluster_makers.jpeg",sep=''), width = length(table(mj012_so$seurat_clusters))*240, height = length(table(mj012_so$seurat_clusters))*480,res=300)
DoHeatmap(mj012_so_cluster.averages , features = c(mj012_top30$gene), group.colors =colors,draw.lines = FALSE,angle = 270)
dev.off()

gamma = c("TRGV2", "TRGV3", "TRGV4", "TRGV5", "TRGV8", "TRGV9", "TRGV11", "TRGV1", "TRGV5P", "TRGV6", "TRGV7", "TRGV10", "TRGVA", "TRGVB, TRGJ1", "TRGJ2", "TRGJP", "TRGJP1", "TRGJP2")
delta = c("TRDV1", "TRDV2", "TRDV3", "TRAV14DV4", "TRAV29DV5", "TRAV23DV6", "TRAV36DV7", "TRAV38-2DV8 TRDJ1", "TRDJ2", "TRDJ3", "TRDJ4")

jpeg("~/Desktop/data/mj012_output/MJ012_DOT_gamma.jpeg", width = 2400, height = 1200, res=300)
DotPlot(mj012_so, features = rev(gamma) , group.by = 'seurat_clusters', cols = c('#DDDDDD','#3F51B5')) + RotatedAxis()
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_DOT_delta.jpeg", width = 1800, height = 1200, res=300)
DotPlot(mj012_so, features = rev(delta) , group.by = 'seurat_clusters', cols = c('#DDDDDD','#3F51B5')) + RotatedAxis()
dev.off()


jpeg(paste("~/Desktop/data/mj012_output/","MJ012","_HEAT_cluster_gd.jpeg",sep=''), width = length(table(mj012_so$seurat_clusters))*240, height = length(table(mj012_so$seurat_clusters))*480,res=300)
DoHeatmap(mj012_so_cluster.averages , features =c(gamma,delta), group.colors =colors,draw.lines = FALSE,angle = 270)
dev.off()


Idents(mj012_so) <- 'seurat_clusters'
jpeg("~/Desktop/data/mj012_output/MJ012_DOT.jpeg", width = 1600, height = 1200, res=300)
DotPlot(mj012_so, features = rev(c('CD3E', 'CD4', 'CD8A', 'TRAC','TRBC1', 'TRDC', 'CD79A', 'CD14','CXCR4', 'ITGA4')), group.by = 'seurat_clusters', cols = c('#DDDDDD','#3F51B5')) + RotatedAxis()
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_DOT_TF.jpeg", width = 1600, height = 1200, res=300)
DotPlot(mj012_so, features = rev(c('TBX21', 'EOMES', 'RORC', 'RUNX3','PRDM1', 'ZNF683', 'TCF7')), group.by = 'seurat_clusters', cols = c('#DDDDDD','#3F51B5')) + RotatedAxis()
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_DOT_CT.jpeg", width = 1600, height = 1200, res=300)
DotPlot(mj012_so, features = rev(c('IFNG', 'TNF', 'IL17A', 'GZMA','GZMB', 'PRF1', 'FASLG','TNFSF10')), group.by = 'seurat_clusters', cols = c('#DDDDDD','#3F51B5')) + RotatedAxis()
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_VLN_pt.jpeg", width = 2500, height = 1000, res=300)
VlnPlot(mj012_so, features = c('CD3E', 'CD4', 'CD8A', 'TRAC','TRBC1', 'TRDC', 'CD79A', 'CD14','CXCR4', 'ITGA4'), ncol = 5, cols = colors,pt.size = 0.05,combine = TRUE, group.by = 'seurat_clusters', same.y.lims = TRUE)
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_VLN.jpeg", width = 2500, height = 1000, res=300)
VlnPlot(mj012_so, features = c('CD3E', 'CD4', 'CD8A', 'TRAC','TRBC1', 'TRDC', 'CD79A', 'CD14','CXCR4', 'ITGA4'), ncol = 5, cols = colors,pt.size = 0,combine = TRUE, group.by = 'seurat_clusters', same.y.lims = TRUE)
dev.off()


Idents(mj012_so) <- 'seurat_clusters'
write.csv(mj012_so.markers,'~/Desktop/data/mj012_output/mj012_markers.csv')
write.csv(mj012_so$orig.ident, '~/Desktop/data/mj012_output/mj012_cell.csv')
write.csv(mj012_so$orig.ident, '~/Desktop/hvg/mj012_cell_post.csv')




info = read.csv('~/Desktop/hvg/mj012_post.csv')
mj012_so$pre_post <- info[,5]
mj012_so$clono <-info[,4]
clono <- c('clonotype1','clonotype2','clonotype3','clonotype4','clonotype5','clonotype6','clonotype7','clonotype8','clonotype9','clonotype10','Others', 'None')
'%ni%' <- Negate('%in%')
clonotype = as.vector(info[,4])

mj012_so$clonotype <- replace(clonotype, seq(length(clonotype))[clonotype %ni% clono], 'Others')
clono_colors <- c('#EEEEEE','#BBBBBB','#FF9800','#FFEB3B','#00BCD4','#2196F3','#3F51B5','#4CAF50','#009688','#673AB7','#9C27B0','#E91E63')

'''
mj012_so$productive <- replace(clonotype, seq(length(clonotype))[clonotype != 'None'], 'Productive')
jpeg("~/Desktop/data/mj012_output/MJ012_Productive.jpeg", width = 1600, height = 900, res=300)
DimPlot(mj012_so, group.by = 'productive' ,pt.size=1)
dev.off()
'''

sel_lab = c()
sel_col = c()
for (i in seq(12)){
    if (clono[i] %in% rownames(table(mj012_so$clonotype)))
        {
        sel_lab = append(sel_lab, clono[i])
        sel_col = append(sel_col, rev(clono_colors)[i])
        }
}

sel_col = rev(sel_col)


jpeg("~/Desktop/data/mj012_output/MJ012_CLONO.jpeg", width = 1600, height = 900, res=300)
DimPlot(mj012_so, group.by = 'clonotype', order = sel_lab, cols =sel_col ,pt.size=1)
dev.off()


info = read.csv('~/Desktop/hvg/mj012_post.csv')
mj012_so$pre_post <- info[,5]
mj012_so$clono <-info[,4]
clono <- c('clonotype1','clonotype2','clonotype3','clonotype4','clonotype5','clonotype6','clonotype7','clonotype8','clonotype9','clonotype10','Others', 'None')
'%ni%' <- Negate('%in%')
clonotype = as.vector(info[,4])
gvh = as.vector(info[,5])
gvh_sel = c('CD4_GvH','CD4_nonGvH','CD8_GvH','CD8_nonGvH','Other_TCR_Clones','Undetectable_for_TCR')
mix <- replace(gvh, seq(length(gvh))[gvh %ni% gvh_sel], 'Undetectable_for_TCR')
mix <- replace(mix, seq(length(gvh))[clonotype != 'None' & gvh %ni% gvh_sel], 'Other_TCR_Clones')

mj012_so$mix <- mix

color6 <- c('#EEEEEE','#BBBBBB',color5[2:5])
sel_lab = c()
sel_col = c()
for (i in seq(6)){
    if (gvh_sel[i] %in% rownames(table(mj012_so$mix)))
        {
        sel_lab = append(sel_lab, gvh_sel[i])
        sel_col = append(sel_col, rev(color6)[i])
        }
}



sel_col = rev(sel_col)
jpeg("~/Desktop/data/mj012_output/MJ012_MIX.jpeg", width = 1600, height = 900, res=300)
DimPlot(mj012_so, group.by = 'mix', order = sel_lab, cols =sel_col ,pt.size=1)
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_MIX_NoLegend.jpeg", width = 1600, height = 900, res=300)
DimPlot(mj012_so, group.by = 'mix', order = sel_lab, cols =sel_col ,pt.size=1)+NoLegend()
dev.off()



Idents(mj012_so) <- 'mix'
mj012_maGvH <- FindAllMarkers(mj012_so, only.pos = TRUE, min.cells.group = 1)
write.csv(mj012_maGvH,'~/Desktop/data/mj012_output/mj012_gvh_markers.csv')

mj012_top30 <- mj012_maGvH %>% group_by(cluster) %>% top_n(n = 12, wt = avg_logFC)



Idents(mj012_so) <- 'pre_post'
mj012_maGvH <- FindallMarkers(mj012_so, ident.1 = 'CD4_GvH',min.cells.group = 1)


mj012_top30 <- mj012_maGvH  %>% group_by(cluster) %>% top_n(n = 12, wt = avg_logFC)
mj012_maGvH[mj012_maGvH$p_val_adj <0.05]




lst <- sort(mj012_so@assays$RNA@data[mj012_so$pre_post == 'CD4_GvH'], index.return=TRUE, decreasing=TRUE)




jpeg("~/Desktop/data/mj012_output/mj012_FP_CD4.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'CD4')
dev.off()

jpeg("~/Desktop/data/mj012_output/mj012_FP_CD8A.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'CD8A')
dev.off()

jpeg("~/Desktop/data/mj012_output/mj012_FP_CD69.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'CD69')
dev.off()

jpeg("~/Desktop/data/mj012_output/mj012_FP_IL2.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'IL2')
dev.off()

jpeg("~/Desktop/data/mj012_output/mj012_FP_ITGAE.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'ITGAE')
dev.off()

jpeg("~/Desktop/data/mj012_output/mj012_FP_FOXP3.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'FOXP3')
dev.off()

jpeg("~/Desktop/data/mj012_output/mj012_FP_CCR7.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'CCR7')
dev.off()

jpeg("~/Desktop/data/mj012_output/mj012_FP_KLF2.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'KLF2')
dev.off()


jpeg("~/Desktop/data/mj012_output/mj012_FP_IL22.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'IL22')
dev.off()


jpeg("~/Desktop/data/mj012_output/MJ012_FP_CD3E.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'CD3E')
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_FP_TRDC.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'TRDC')
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_FP_CD79A.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'CD79A')
dev.off()

jpeg("~/Desktop/data/mj012_output/MJ012_FP_CD14.jpeg", width = 1600, height = 900, res=300)
FeaturePlot(mj012_so,features = 'CD14')
dev.off()








