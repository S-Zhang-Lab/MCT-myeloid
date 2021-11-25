##T cell subset analysis Figure 5

library(dplyr)
library(Seurat)
library(ggplot2)
library(umap)
library(ggplot2)
library(reshape2)
library(psych)
library(ggpubr)
library(wesanderson)
library(plotly)


#######################################################
######Current CD3 gating for Figure5#####Nov42021
###For newer T cell analysis
###Alternatively set HTO_classification_group as the Idents 

load(file= "MCT.singlet_small_analysis_manuscript.rda")



Idents(MCT.singlet) <- MCT.singlet@meta.data$HTO_classification_group


DefaultAssay(MCT.singlet) <- "ADT"


Immunecells_Tcells <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD45", feature2 = "CITE-CD3", cells = CD45hptprch)
CD45CD3Tcells <- CellSelector(plot = Immunecells_Tcells)


Idents(MCT.singlet, cells = CD45CD3Tcells) <- "CD45CD3Tcells"

TcellCD3subset <- subset(MCT.singlet, idents = "CD45CD3Tcells")

table(TcellCD3subset$HTO_classification_group)

##




####

l <- c(ls())
l
l <- l[-172]
rm(list = (l))

save.image("TcellCD3subset_small.rda")

Idents(TcellCD3subset) <- "HTO_classification_group"

TcellCD3subset_small <- subset(TcellCD3subset, downsample = 256)
table(TcellCD3subset_small$HTO_classification_group)

#Control MCT_treatment 
#256          256 


TcellCD3subset_small <- NormalizeData(TcellCD3subset_small, 
                                   assay = "RNA", features = rownames(TcellCD3subset_small))

TcellCD3subset_small <- FindVariableFeatures(TcellCD3subset_small,
                                          assay = "RNA", features = rownames(TcellCD3subset_small))



### standard scaling (no regression)


TcellCD3subset_small <- ScaleData(TcellCD3subset_small,
                               assay = "RNA", features = VariableFeatures(TcellCD3subset_small))

DefaultAssay(TcellCD3subset_small) <- "RNA"

TcellCD3subset_small <- RunPCA(TcellCD3subset_small, features = VariableFeatures(object = TcellCD3subset_small))

TcellCD3subset_small <- FindNeighbors(TcellCD3subset_small, dims = 1:10)
TcellCD3subset_small <- FindClusters(TcellCD3subset_small, resolution = 0.8)

TcellCD3subset_small <- RunUMAP(TcellCD3subset_small, dims = 1:10)  

Idents(TcellCD3subset_small) <- "seurat_clusters"
pdf("Dimplot_Seurat_umap_TcellCD3subset_small.pdf", width = 12, height = 6, paper = 'special')
DimPlot(TcellCD3subset_small, split.by = "HTO_classification_group", reduction = "umap")
dev.off()


DefaultAssay(TcellCD3subset_small) <- "RNA"

TcellCD3subset_small_markers <- FindAllMarkers(TcellCD3subset_small, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TcellCD3subset_small_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10Tcell

DoHeatmap(TcellCD3subset_small, features = top10Tcell$gene)





#pdf("Rplot-heatmap_TcellCD3subset_small.pdf", width = 12, height = 6, paper = 'special')
#DoHeatmap(TcellCD3subset_small, features = unique(TcellCD3subset_small_markers$gene))
#dev.off()




ggplot(TcellCD3subset_small@meta.data, aes(seurat_clusters, group = HTO_classification_group)) +
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes(label = scales::percent(..prop..),
                
                y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="seurat_clusters") +
  facet_grid(~HTO_classification_group) +
  scale_y_continuous(labels = scales::percent)


freq_table_T <- prop.table(x = table(TcellCD3subset_small@meta.data$seurat_clusters,
                                     TcellCD3subset_small@meta.data$HTO_classification_group),margin = 2)

library("scales")

coloridentities <- levels(TcellCD3subset_small@meta.data$seurat_clusters)

my_color_palette <- hue_pal()(length(coloridentities))

barplot(height = freq_table_T, col = my_color_palette)

Idents(TcellCD3subset_small) <- TcellCD3subset_small$seurat_clusters

DoHeatmap(TcellCD3subset_small, features = c("Cd3e", "Cd8a", "Cd4", "Tcf7", "IL7r", "Ccr7", "Il17r", "Cxcr6",
                                                 "Gzmk","Lag3", "Layn", "Havrc2", "Ctla4","Tox", "Pdcd1","Cx3cr1",
                                                 "Fcgr3a", "Prf1", "Klrg1", "Gitr", "Cd69", "Cxcl13", "Cd103", "Cd39", "Tcf1"), slot = "data") +
  scale_fill_gradientn(colors = c("white", "red"))




### overlay cluster 0 and 2 with CITE-markers

Idents(TcellCD3subset_small) <- "seurat_clusters"

Cluster0TcellCD3subset_small <- subset(TcellCD3subset_small, idents = c("0"))

table(Cluster0TcellCD3subset_small$HTO_classification_group)

#Control MCT_treatment 
#47          74

DefaultAssay(Cluster0TcellCD3subset_small) <- "ADT"

Immunecells_Cluster0Tcells <- FeatureScatter(Cluster0TcellCD3subset_small, feature1 = "CITE-CD4", 
                                              feature2 = "CITE-CD3")
CD45CD3CD4Tcells <- CellSelector(plot = Immunecells_Cluster0Tcells)

Idents(Cluster0TcellCD3subset_small, cells = CD45CD3CD4Tcells) <- "CD4CD3Cluster0Tcellsubsetsmall"
CD4CD3Cluster0Tcellsubsetsmall.subset <- subset(Cluster0TcellCD3subset_small, 
                                                 idents = "CD4CD3Cluster0Tcellsubsetsmall")


table(CD4CD3Cluster0Tcellsubsetsmall.subset$HTO_classification_group)
#    Control MCT_treatment 
#19           24 


DefaultAssay(Cluster0TcellCD3subset_small) <- "RNA"





Idents(TcellCD3subset_small) <- "seurat_clusters"
Cluster0vsothersTcellCD3subset_small.markers <- FindMarkers(TcellCD3subset_small, 
                                            ident.1 = "0", ident.2 = c("1", "2", "3","4","5","6","7","8","9"), only.pos = FALSE)
write.csv(Cluster0vsothersTcellCD3subset_small.markers, file = "DGE_Cluster0vsothersTcellCD3subset_small.markers.csv")

####Figure 
volcanodata.DGE_Cluster0vsothersTcellCD3subset_small.markers <- read.csv("DGE_Cluster0vsothersTcellCD3subset_small.markers.csv",
                                   header=TRUE, stringsAsFactors=FALSE, row.names=1)

volcanodata.DGE_Cluster0vsothersTcellCD3subset_small.markers$p.adj <- -1*log(
  volcanodata.DGE_Cluster0vsothersTcellCD3subset_small.markers$p_val_adj, 10)

#convert p value into log values

pThresh <- 20

#plot Figure Supplementary 3D
ggplot(volcanodata.DGE_Cluster0vsothersTcellCD3subset_small.markers, aes(x=avg_logFC, y=p.adj , cex=0.5)) + 
  
  geom_point(shape=19, size=3) 

#to know the names of genes
ggplotly(ggplot(volcanodata.DGE_Cluster0vsothersTcellCD3subset_small.markers, aes(x=avg_logFC, y=p.adj, cex=0.5, 
                                                         text = row.names(volcanodata.DGE_Cluster0vsothersTcellCD3subset_small.markers))) + 
           geom_point(shape=21, size=3) + 
           theme(panel.background = element_rect(fill=NA), panel.grid.major = element_line(colour = "grey80")))


#whole T cell susbet
VlnPlot(TcellCD3subset_small, features = "Il17a", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

DefaultAssay(TcellCD3subset_small) <- "ADT"


VlnPlot(TcellCD3subset_small, features = "CITE-CD4", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

VlnPlot(TcellCD3subset_small, features = "CITE-CD8a", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(TcellCD3subset_small, features = "CITE-CD8a", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(TcellCD3subset_small, features = "CITE-CD25", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()



VlnPlot(TcellCD3subset_small, features = "CITE-CD44", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(TcellCD3subset_small, features = "CITE-PD-1", group.by = "seurat_clusters", 
        split.by = "HTO_classification_group", split.plot = TRUE) +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()





DefaultAssay(TcellCD3subset_small) <- "RNA"

VlnPlot(TcellCD3subset_small, features = "Il10", group.by = "seurat_clusters", 
        split.by = "HTO_classification_group", split.plot = TRUE) +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()



VlnPlot(TcellCD3subset_small, features = "Il2", group.by = "seurat_clusters", 
        split.by = "HTO_classification_group", split.plot = TRUE) +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()



VlnPlot(TcellCD3subset_small, features = "Lag3", group.by = "seurat_clusters", 
        split.by = "HTO_classification_group", split.plot = TRUE) +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()




VlnPlot(TcellCD3subset_small, features = "Ctla4", group.by = "seurat_clusters", 
        split.by = "HTO_classification_group", split.plot = TRUE) +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()



VlnPlot(TcellCD3subset_small, features = "Pdcd1", group.by = "seurat_clusters", 
        split.by = "HTO_classification_group", split.plot = TRUE) +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()
Idents(TcellCD3subset_small) <- "HTO_classification_group"

VlnPlot(TcellCD3subset_small, features = "Pdpn", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

RidgePlot(TcellCD3subset_small, features = "Pdpn")


Idents(TcellCD3subset_small) <- "HTO_classification_group"

RidgePlot(TcellCD3subset_small, features = c("Tox", "Pdpn", "Lgals3", "Cd44", "Gzmk", "Cxcr6", "Ctla4", "Pdcd1", "Lag3"))

#Cluster0TcellCD3subset_small subset

DefaultAssay(Cluster0TcellCD3subset_small) <- "ADT"


VlnPlot(Cluster0TcellCD3subset_small, features = "CITE-CD4", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "CITE-CD8a", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "CITE-CD3", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "CITE-CD25", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()



VlnPlot(Cluster0TcellCD3subset_small, features = "CITE-CD44", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "CITE-PD-1", group.by = "HTO_classification_group", 
        split.by = "HTO_classification_group", split.plot = TRUE) +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


DefaultAssay(Cluster0TcellCD3subset_small) <- "RNA"

VlnPlot(Cluster0TcellCD3subset_small, features = "Il10", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()



VlnPlot(Cluster0TcellCD3subset_small, features = "Il2", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()



VlnPlot(Cluster0TcellCD3subset_small, features = "Lag3", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "Tox", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "Ccl7", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "Ccl12", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "Lgals3", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Cluster0TcellCD3subset_small, features = "Pdpn", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

Idents(Cluster0TcellCD3subset_small) <- "HTO_classification_group"
RidgePlot(Cluster0TcellCD3subset_small, features = c("Tox", "Pdpn", "Lgals3", 
                                             
                                             "Cd44", "Gzmk", "Cxcr6", "Ctla4", "Pdcd1", "Lag3"))


DefaultAssay(TcellCD3subset_small) <- "RNA"

RidgePlot(TcellCD3subset_small, features = c("Tox", "Pdpn", "Lgals3", 
                                             
                                             "Cd44", "Gzmk", "Cxcr6", "Ctla4", "Pdcd1", "Lag3"))



DefaultAssay(TcellCD3subset_small) <- "ADT"

RidgePlot(TcellCD3subset_small, features = c( "CITE-CD8a", "CITE-CD4", 
                                              "CITE-CD44",
                                             "CITE-CD25"), ncol = 2)


DefaultAssay(Cluster0TcellCD3subset_small) <- "ADT"
RidgePlot(Cluster0TcellCD3subset_small, features = c( "CITE-CD8a", "CITE-CD4", 
                                              "CITE-CD44","CITE-CD25"), ncol = 2)






#####################

Cluster1and2cellCD3subset_small <- subset(TcellCD3subset_small, idents = c("1", "2"))




DefaultAssay(Cluster1and2cellCD3subset_small) <- "RNA"

Cluster1and2cellCD3subset_small <- RunPCA(Cluster1and2cellCD3subset_small, 
                                          features = VariableFeatures(object = Cluster1and2cellCD3subset_small))

Cluster1and2cellCD3subset_small <- FindNeighbors(Cluster1and2cellCD3subset_small, dims = 1:10)
Cluster1and2cellCD3subset_small <- FindClusters(Cluster1and2cellCD3subset_small, resolution = 0.8)

Cluster1and2cellCD3subset_small <- RunUMAP(Cluster1and2cellCD3subset_small, dims = 1:10)  




Idents(Cluster1and2cellCD3subset_small) <- "seurat_clusters"

DimPlot(Cluster1and2cellCD3subset_small, split.by = "HTO_classification_group", reduction = "umap")



DefaultAssay(Cluster1and2cellCD3subset_small) <- "RNA"

Cluster1and2cellCD3subset_small_markers <- FindAllMarkers(Cluster1and2cellCD3subset_small, 
                                                          only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cluster1and2cellCD3subset_small_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10Tcellsubsetted

DoHeatmap(Cluster1and2cellCD3subset_small, features = top10Tcellsubsetted$gene)



Idents(Cluster1and2cellCD3subset_small) <- "HTO_classification_group"
Cluster1and2cellCD3subset_small.MCT.markers <- FindMarkers(Cluster1and2cellCD3subset_small, 
                                                            ident.1 = "MCT_treatment", 
                                                           ident.2 = "Control", only.pos = FALSE)
write.csv(Cluster1and2cellCD3subset_small.MCT.markers, file = "Cluster1and2cellCD3subset_small.MCT.markers.csv")

####Figure 
volcanodata.Cluster1and2cellCD3subset_small.MCT.markers <- read.csv("Cluster1and2cellCD3subset_small.MCT.markers.csv",
                                                                         header=TRUE, stringsAsFactors=FALSE, row.names=1)

volcanodata.Cluster1and2cellCD3subset_small.MCT.markers$p.adj <- -1*log(
  volcanodata.Cluster1and2cellCD3subset_small.MCT.markers$p_val_adj, 10)

#convert p value into log values

pThresh <- 20

#plot Figure Supplementary 
ggplot(volcanodata.Cluster1and2cellCD3subset_small.MCT.markers, aes(x=avg_logFC, y=p.adj , cex=0.5)) + 
  
  geom_point(shape=19, size=3) 

#to know the names of genes
ggplotly(ggplot(volcanodata.Cluster1and2cellCD3subset_small.MCT.markers, aes(x=avg_logFC, y=p.adj, cex=0.5, 
                                                                                  text = row.names(
                                                                                    volcanodata.Cluster1and2cellCD3subset_small.MCT.markers))) + 
           geom_point(shape=21, size=3) + 
           theme(panel.background = element_rect(fill=NA), panel.grid.major = element_line(colour = "grey80")))



VlnPlot(Cluster1and2cellCD3subset_small, features = "Pdcd1", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()

VlnPlot(Cluster1and2cellCD3subset_small, features = "Ctla4", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()


VlnPlot(Cluster1and2cellCD3subset_small, features = "Tox", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()

VlnPlot(Cluster1and2cellCD3subset_small, features = "Foxp3", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()

VlnPlot(Cluster1and2cellCD3subset_small, features = "Gzmb", group.by = "HTO_classification_group") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = "t.test")+
  theme_classic()


# Heatmap for subsetted CD3e high T cells plotted for marker genes and gene markers for functional status
Idents(Cluster1and2cellCD3subset_small) <- Cluster1and2cellCD3subset_small$HTO_classification_group
Idents(Cluster1and2cellCD3subset_small) <- Cluster1and2cellCD3subset_small$seurat_clusters

DoHeatmap(Cluster1and2cellCD3subset_small, features = c("Cd8a", "Cd4", "Tcf7", "IL7r", "Ccr7", "Il17r", "Cxcr6",
                                             "Gzmk","Lag3", "Layn", "Havrc2", "Ctla4","Tox", "Pdcd1","Cx3cr1",
                                             "Fcgr3a", "Prf1", "Klrg1", "Gitr", "Cd69", "Cxcl13", "Cd103", "Cd39", "Tcf1",
                                             "Nkg7", "Ctsw", "Cd8b1", "Dusp2", "Ccl5","Ly6c2", "Xcl1", "Plac8", "Gzmb",
                                             "Ramp3", "Tmem176b", "Tmem176a", "S1pr1", "Pxdc1", "Hif1a", "Krt83", "Tnfsf8",
                                             "Klf2", "Itgav", "Foxp3", "Tff1", "Tnfrsf4", "Klrg1", "Ikzf2", "Tnfrsf9", 
                                             "Bmyc", "Ikzf4", "Rel"), slot = "data") +
  scale_fill_gradientn(colors = c("white", "red"))

DotPlot(Cluster1and2cellCD3subset_small, features = c("Cd8a", "Cd4", "Tcf7", "Ccr7", "Il17r", "Cxcr6",
                                                      "Gzmk","Lag3", "Layn", "Havrc2", "Ctla4","Tox", "Pdcd1","Cx3cr1",
                                                      "Prf1", "Klrg1", "Cd69", "Cxcl13","Nkg7", "Ctsw", "Cd8b1", "Dusp2", "Ccl5","Ly6c2",
                                                      "Xcl1", "Plac8", "Gzmb",
                                                      "Ramp3", "Tmem176b", "Tmem176a", "S1pr1", "Pxdc1", "Hif1a", "Krt83", "Tnfsf8",
                                                      "Klf2", "Itgav", "Foxp3", "Tff1", "Tnfrsf4", "Ikzf2", "Tnfrsf9", 
                                                      "Bmyc", "Ikzf4", 
                                                      "Rel"),
        split.by = "HTO_classification_group") + RotatedAxis()

save.image("TcellCD3subset_small.rda")

###new-lines of code
Cluster1and2cellCD3subset_small$HTO_classification_group<- paste(Idents(Cluster1and2cellCD3subset_small),
                                               Cluster1and2cellCD3subset_small$HTO_classification_group, sep = "_")# redefine Identity which will include the                                                                                                                             #information you want to split the plot by


new.idents <- Cluster1and2cellCD3subset_small$HTO_classification_group# set new identity
Idents(object = Cluster1and2cellCD3subset_small) <- new.idents# call the seurat object by new identity

levels(Cluster1and2cellCD3subset_small) <- c('0_Control','0_MCT_treatment','1_Control','1_MCT_treatment','2_Control','2_MCT_treatment')

#levels(Cluster1and2cellCD3subset_small) <- c('Control_0','MCT_treatment_0','Control_1','MCT_treatment_1','Control_2','MCT_treatment_2')
#set new subgroups for your seurat object using levels function


DotPlot(Cluster1and2cellCD3subset_small, features = c("Cd8a", "Cd4", "Tcf7", "Ccr7", "Il17r", "Cxcr6",
                                                      "Gzmk","Lag3", "Layn", "Havrc2", "Ctla4","Tox", "Pdcd1","Cx3cr1",
                                                      "Prf1", "Klrg1", "Cd69", "Cxcl13","Nkg7", "Ctsw", "Cd8b1", "Dusp2", "Ccl5","Ly6c2",
                                                      "Xcl1", "Plac8", "Gzmb",
                                                      "Ramp3", "Tmem176b", "Tmem176a", "S1pr1", "Pxdc1", "Hif1a", "Krt83", "Tnfsf8",
                                                      "Klf2", "Itgav", "Foxp3", "Tff1", "Tnfrsf4", "Ikzf2", "Tnfrsf9", 
                                                      "Bmyc", "Ikzf4", 
                                                      "Rel"),
        ) + RotatedAxis()

#smaller gene set for final supplementary panel "Cd8a", "Cd4",

DotPlot(Cluster1and2cellCD3subset_small, features = c( "Ccr7", "Tcf7", "Ikzf4","Foxp3", "Klrg1", "Gzmb",
                                                       "Ctla4", "Klf2", "Tox", 
                                                      "Lag3", "Pdcd1"),) + RotatedAxis()

