library(Seurat)
library(dplyr)
library(gdata)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(SummarizedExperiment)
library(BiocParallel)
library(openxlsx)
library(magrittr)
#library(monocle3)
library(tidyr)
library(patchwork)
#library(SeuratWrappers)
library(viridis)
library(pheatmap)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

## Marker expression

genes <- Hmisc::Cs(PDGFRA, MEST, COL5A2, MPZ, S100B, SOX10, ACAN, S100A1, KCNMA1, KRT14, COL17A1, 
                   KRT15, MLANA, PMEL, KIT, MKI67, TOP2A, CDK1, TUBB3, UCHL1, SOX2, PPARG, LPL, PLIN4, 
                   CD93, PECAM1, LYVE1)
genes <-  c("PTGDS","MEF2C","MFAP5","ASPN","IGF1","ITGA8","MYLK")
#   xlim(-12,15) & ylim(-12,15)


#Figure_2a
tiff("All_UMAP.tiff",units="in", width=7, height=6, res=140)
DimPlot(object = obj_sc2,order = TRUE, pt.size = 0.6, 
        cols = c("lightskyblue3","deeppink","skyblue4","deepskyblue","violet","deeppink3","green4","violetred2", "orange3")) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=7, face = "bold")) & xlim(-12,15) & ylim(-12,15) &
  guides(color = guide_legend(override.aes = list(size=5,shape =18), keyheight = 0.4,keywidth = 0.7, ncol=1)) & 
  labs(title = "All samples")
dev.off()

## colors
c("#EB6154", "#1D91C0", "#238443", "#316675FF",  "#FFC0DA",
  "#D95F02","black", "#FF800E","#f6d746"  ,"#FF00FFFF","blue","#A6D854", "#898989")


## Figure_2b Markers dotplot
tiff("Markers4e.tiff", units="in", width=9, height=5, res=100)
DotPlot(obj_sc2, features = genes,scale =  FALSE,cluster.idents = FALSE) &
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 8),
        legend.title.position = "left",
        legend.title = element_text(angle = 90,size = 8)) &
  scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) &
  ylab("Clusters") 
dev.off() 

#Figure_2g
tiff("Fib_Markers1e.tiff", units="in", width=3.5, height=3.6, res=110)
ggplot(fib_dot, aes(x=features.plot, y= id)) +
  geom_point(aes(color = `Average Expression`,size = `Percent Expressed`), shape=19) +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", 
                                    fill=NA),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.key.size = unit(0.30, "cm"),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.title.position = "top",
        legend.title = element_text(hjust = 0.5,size = 7), 
        aspect.ratio = 1) +
  scale_color_gradientn(colors=paletteer_c("ggthemes::Classic Red-Blue", 30,direction = -1)) +
  ylab("Clusters") 
dev.off() 
#scale_color_gradientn(colors=paletteer_c("ggthemes::Classic Red-Blue", 30,direction = -1)) 



#Figure 2e and 2i
tiff("Avg_2.tiff", units="in", width=14, height=8, res=120)
FeaturePlot(obj_t, features = genes[c(i+2,i+3)],pt.size =1.2,ncol = 2,order = TRUE) & theme(text = element_text(size = 10),
                                                                                            axis.title = element_blank(),
                                                                                            axis.line = element_blank(),
                                                                                            axis.text.y = element_blank(),
                                                                                            axis.text.x = element_blank(),
                                                                                            axis.ticks.x=element_blank(),
                                                                                            axis.ticks.y = element_blank(),
                                                                                            aspect.ratio = 1,
                                                                                            legend.position = c(.15,0.15), #(.855,0.85)
                                                                                            legend.key.size = unit(0.30, "cm"),
                                                                                            legend.text = element_text(size = 8),
                                                                                            plot.background = element_rect(fill = "grey"))  & 
  xlim(-10,10) & ylim(-10,10) & scale_colour_gradientn(colours = paletteer_c("grDevices::Oslo", 100,direction = 1))
dev.off()  




#Figure_2F
tiff("Fib_UMAP.tiff",units="in", width=7, height=6, res=140)
DimPlot(obj_t,order = TRUE, pt.size = 0.6, cols = c("green4", "deepskyblue4","red3","grey35", "deepskyblue","orangered")) &
  theme(text = element_text(size = 5),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=7, face = "bold")) & xlim(-10,8) & ylim(-10,8) &
  guides(color = guide_legend(override.aes = list(size=5,shape =18), keyheight = 0.4,keywidth = 0.7, ncol=1)) & 
  labs(title = "All samples")
dev.off()




#ISG15
tiff("ISG15_120dpi.tiff",units="in", width=14, height=9, res=120)
FeaturePlot(subset(obj_t, stim == "STIM2"),ncol=2,features = c("ISG15"),order = TRUE,pt.size = 1.4) &
  theme(text = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1,
        legend.position = c(.15,0.15), #(.855,0.85)
        legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 8),
        plot.background = element_rect(fill = "grey")) & xlim(-10,8) & ylim(-10,8) &
  scale_colour_gradientn(colours = paletteer_c("grDevices::Oslo", 100,direction = 1),limit=c(0,5.1))
dev.off()

#cycling
tiff("cyc_immunity.tiff",units="in", width=14, height=9, res=70)
FeaturePlot(obj_sc2_cyc,split.by = "stim", ncol=2,features = c("immunity1"),order = TRUE,pt.size = 4.9) &
  theme(text = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1,
        legend.position = c(.95,.75), #(.855,0.85)
        legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 7),
        plot.background = element_rect(fill = "grey")) & xlim(-10,10) & ylim(-10,10) &
  #scale_colour_gradientn(colours =c("#f8f8f8","#f8baba","#f87c7c","#cc0000","#b00000","#960000"))
  scale_colour_gradientn(colours = paletteer_c("grDevices::Oslo", 100,direction = 1),limit=c(-0.4,2))
dev.off()


##counts
tiff("All_counts1.tiff", units="in", width=5, height=9, res=80)
ggplot(all_counts1, aes(x=variable, y = value,fill=cell_types)) +
  geom_bar(stat="identity",position = "fill") + 
  scale_fill_manual(values=c("lightskyblue3","deeppink","skyblue4","deepskyblue","violet","deeppink3","green4","violetred2", "orange3"))
dev.off()



"#FFD1BA"
ids <- c("mesench_fib","pap_fib","reticular_fib","pap_fib","pap_fib","pap_fib","pap_fib","reticular_fib","pap_fib","pap_fib","pap_fib","pap_fib","pap_fib","DP_DC","DS_fib")
Idents(obj_sc2) <- "seurat_clusters"
names(ids) <- levels(obj_sc2)
obj_t <- RenameIdents(obj_sc2, ids)



# Figure 3c  viral expression
tiff("Fib_viral2_120dpi.tiff",units="in", width=14, height=8, res=120)
FeaturePlot(subset(obj_t,stim == "STIM2"),keep.scale = "feature",features = c("MT-early","MT-late"),order = TRUE,pt.size = 1.2) &
  theme(text = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1,
        legend.position.inside = c(0.9,0.87),
        legend.key.size = unit(0.30, "cm"),
        legend.text=element_text(size=rel(1.2)),
        plot.background = element_rect(fill = "grey")) & xlim(-10,8) & ylim(-10,8) &
  scale_colour_gradientn(colours =c("#f8f8f8","#f8baba","#f87c7c","#cc0000","#b00000","#960000"),limit=c(0,7.5))
#scale_colour_gradientn(colours = paletteer_c("grDevices::Oslo", 100,direction = 1))
dev.off()

# Figure 4c Module score
tiff("Fib_immunity_120dpi_2.tiff",units="in", width=14, height=8, res=150)
FeaturePlot(obj_t,split.by = "stim",features = c("immunity1"),keep.scale = "all",order = TRUE,pt.size = 1.2) &
  theme(text = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.ticks.y = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.1,0.2),
        legend.key.size = unit(0.30, "cm"),
        plot.background = element_rect(fill = "grey")) & xlim(-10,8) & ylim(-10,8) & 
  scale_colour_gradientn(colours = paletteer_c("grDevices::Oslo", 100,direction = 1),limit=c(-0.4,2))
#scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),limit=c(-1,2)) 
dev.off()

####Figure  4d and 3d
coex_df <- fib_coex_df2

x1 <- coex_df[which(rowSums(coex_df[,c(4,6,7)])==0),]
x2 <- coex_df[which(rowSums(coex_df[,c(4,6,7)])!=0),]


#x1 <- coex_df[which(rowSums(coex_df[,c(5:7)])==0),]
#x2 <- coex_df[which(rowSums(coex_df[,c(5:7)])!=0),]

tiff("Fib_viral_log2earl-late_120dpi_2.tiff",units="in", width=10, height=8, res=120)

tiff("Fib_viral_log2earl-ISG15_120dpi_2.tiff",units="in", width=10, height=8, res=120)
#scale_colour_gradientn(colours = paletteer_c("grDevices::Zissou 1", 800,direction = -1), limit=c(-2,2))
ggplot() +    geom_point(data=x1,aes(x = umapharmony_1, y = umapharmony_2),colour = "#EAEBEB", size=1.9, show.legend = TRUE) + 
  geom_point(data=x2,aes(x = umapharmony_1, y = umapharmony_2,colour=fc1), size=1.9, show.legend = TRUE) +
  scale_colour_gradientn(colours = viridis::plasma(100, direction = -1),limit=c(-1.5,1.5),na.value = "#F0F921FF" ) + 
  theme(text = element_text(size = 9), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text.y = element_blank(),  
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),aspect.ratio = 1,panel.background = element_rect(fill = "grey"),
        legend.text = element_text(size=7, face = "bold")) & xlim(-10,8) & ylim(-10,8)
dev.off()

tiff("Cyc_viral_log2earl-late_120dpi.tiff",units="in", width=3, height=3, res=90)
ggplot() +    geom_point(data=x1,aes(x = umapharmony_1, y = umapharmony_2),colour = "#EAEBEB", size=4.9, show.legend = TRUE) + 
  geom_point(data=x2,aes(x = umapharmony_1, y = umapharmony_2,colour=log2fc), size=4.9, show.legend = TRUE) +
  scale_colour_gradientn(colours = viridis::plasma(100, direction = -1),limit=c(-1.5,1.5),na.value = "#F0F921FF" ) +
  #scale_colour_gradientn(colours = brewer.pal(n = 11, name = "Spectral") , limit=c(-1.5,1.5),na.value = "#9E0142" ) + 
  theme(text = element_text(size = 9), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text.y = element_blank(),  
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.ticks.y = element_blank(),aspect.ratio = 1,panel.background = element_rect(fill = "grey"),
        legend.text = element_text(size=7, face = "bold"))
dev.off()




####### RNA-seq ##############

# PCA
tiff("PCA_RNAseq4.tiff",units="in", width=6, height=6, res=120)
ggplot(PCA_data, aes(x=PC1,y=PC2, label = rownames(PCA_data), color=condition)) +  
  geom_point(size=5.1, stroke = 1) +   geom_text_repel(size=3.15,stat = "identity") +
  theme(axis.text.x  = element_text(size = 3.5,margin = margin(r = 0))) +
  theme(axis.text.y = element_text(size = 3.5,margin = margin(r = 0))) +
  xlab(paste0("PC1, VarExp:", round(percentVar[1],4))) + 
  ylab(paste0("PC2, VarExp:", round(percentVar[2],4))) + 
  theme(axis.title.y = element_text(size = 10.5))+
  theme(axis.title.x = element_text(size = 10.5))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),aspect.ratio = 1,
        panel.border = element_rect(colour = "black", fill=NA)) + #, size=0.2
  scale_color_manual(values=c("#e95462","#319dF6AD","#2700D1"))  + scale_shape_manual(values = c(19,19)) +
  theme(legend.position="none") + ylim(-40,70) + xlim(-40, 70) 
dev.off()


## heatmap
tiff("HM_RNAseq2.tiff",units="in", width=7, height=11, res=90)
#pheatmap(hm_vsn[sig_genes_names$ensembl_id,],
#         scale = "row",clustering_distance_rows="correlation", show_rownames = FALSE, cluster_cols = FALSE)
pheatmap(y,scale = "row",cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = FALSE, cellheight=1, cellwidth = 40) 
dev.off()

pheatmap(hm_vsn[genes$ensembl_id,],
         scale = "row",clustering_distance_rows="correlation", show_rownames = FALSE, cluster_cols = FALSE,treeheight_row = 0)

#Volcano
tiff("Volcano2_100dpi.tiff",units="in", width=6, height=5, res=120)
ggplot(DESeq2Res_plot, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=genes)) +
  geom_point(size=2.5,alpha=0.6) +  scale_color_manual(values=c("#1F74B1FF", "black", "#EC6857FF"))  +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2), aspect.ratio = 1) +
  labs(title = "FDR = 0.05") + xlim(-5,12) + geom_vline(xintercept = c(-1,1),linetype="dashed",color="darkgrey")  + 
  geom_hline(yintercept =  -log10(0.05),linetype="dashed",color="darkgrey") +   
  geom_text_repel(aes(label=genes),size=2.7, max.overlaps = Inf,colour="#FF410D",fontface = "bold.italic")
dev.off()

## GO

tiff("GO_clusterprofiler.tiff", units="cm", width=17, height=19, res=90)
clusterProfiler::dotplot(CP_res,showCategory= 30, font.size=9)
dev.off()

## Figur 3f Viral expression plot
#MT-early
viral_meta <- merge(obj_sc2@meta.data,t(obj_sc2@assays$RNA$data[c("MT-early","MT-late"),]),by='row.names')
#all_cells
tiff("All_Viral_expression_norm.tiff", units="in", width=12, height=6, res=150)
ggplot(subset(viral_meta,`MT-early` > 0), aes(x = seurat_annotations, y = `MT-early`,color=seurat_annotations)) +
  geom_beeswarm(size=0.6,  corral = "gutter") + 
  scale_color_manual(values = c("lightskyblue3","deeppink","skyblue4","deepskyblue","violet","deeppink3","green4","violetred2", "orange3")) +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),aspect.ratio = 0.5,axis.title.x = element_blank()) +
  theme(legend.position="none")  +
  ylim(0,8) + 
  ggplot(subset(viral_meta,`MT-late` > 0), aes(x = seurat_annotations, y = `MT-late`,color=seurat_annotations)) +
  geom_beeswarm(size=0.6,  corral = "gutter") + 
  scale_color_manual(values = c("lightskyblue3","deeppink","skyblue4","deepskyblue","violet","deeppink3","green4","violetred2", "orange3")) +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),aspect.ratio = 0.5,axis.title.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=5,shape =18), keyheight = 0.4,keywidth = 0.7, ncol=1)) +
  ylim(0,8) 
dev.off()


#fib 3f
tiff("Fib_viral_expression_norm.tiff", units="in", width=12, height=6, res=140)
ggplot(subset(viral_meta,MT.early > 0), aes(x = annotations, y = MT.early,color=annotations)) +
  geom_beeswarm(size=0.6,  corral = "gutter") + scale_color_manual(values =c("green4", "deepskyblue4","red3","grey35", "deepskyblue","orangered","deeppink3")) +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),aspect.ratio = 0.5,axis.title.x = element_blank()) +
  theme(legend.position="none")  +
  ylim(0,8) + 
  ggplot(subset(viral_meta, MT.late > 0), aes(x = annotations, y = MT.late,color=annotations)) +
  geom_beeswarm(size=0.6,  corral = "gutter") + scale_color_manual(values = c("green4", "deepskyblue4","red3","grey35", "deepskyblue","orangered","deeppink3")) +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 0.0),
        panel.grid.minor = element_line(color="grey",size = 0.1),
        panel.grid.major = element_line(color="grey",size = 0.1),
        axis.ticks = element_blank(),aspect.ratio = 0.5,axis.title.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=5,shape =18), keyheight = 0.4,keywidth = 0.7, ncol=1)) +
  ylim(0,8) 
dev.off()



















