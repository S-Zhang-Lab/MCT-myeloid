

####Figure 2 2021 Submission ####################




###loadlibraries_important

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



load(file= "MCT.singlet_small_analysis_manuscript.rda") ## load for ADT related 
load(file = "MCT.singlet_small.rda") ## for RNA assay



### load 10X data (NOTE: outs folder from Cellranger V3.0.3)


MCT_data <- Read10X(data = "filtered_feature_bc_matrix")


### examine the structure of the data. It should be a list of 2 objects

str(MCT_data) 


### Lets examine the data

head(MCT_data)


### Initialize the Seurat object with the raw (non-normalized data).

MCT <- CreateSeuratObject(counts = MCT_data$"Gene Expression")



### add ADT objects

MCT$CITE <- CreateAssayObject(counts = 
                                MCT_data$"Antibody Capture")



### read in TotalA.feature reference file

TotalA.features = read.csv(file = "TotalA_Ref.csv")

HTOs <- c("M-HTO-1","M-HTO-2","M-HTO-3","M-HTO-4","M-HTO-5","M-HTO-6")




### load antibody names and replace "_" with "-", since it has been changed when load data in Seurat

CITE_abs <- TotalA.features$name
CITE_abs <- gsub("_", "-", CITE_abs)
CITE_ADT <- CITE_abs[grep("^CITE",CITE_abs)]


###extract HTOs reads  # change to [1:6] for 6 sample demultiplexing


MCT.HTOs <- as.data.frame(MCT$"CITE"@counts)
HTOs <- t(as.data.frame(MCT.HTOs))
HTOs[1:6,1:6]


### check dim of HTOs table, HTO colunm should at very end

dim(HTOs)  
colnames(HTOs)
MCT.HTOs <- HTOs[,35:40]
MCT.HTOs[1:4,1:6]
MCT.ADT <- HTOs[,1:34]


### Add HTO data as a new assay independent from RNA

MCT[["HTO"]] <- CreateAssayObject(counts = t(MCT.HTOs))
MCT[["ADT"]] <- CreateAssayObject(counts = t(MCT.ADT))


### Normalize HTO data, here we use centered log-ratio (CLR) transformation

MCT <- NormalizeData(MCT, assay = "HTO", normalization.method = "CLR")
MCT <- NormalizeData(MCT, assay = "ADT", normalization.method = "CLR")


### Demultiplex cells based on HTO enrichment

MCT <- HTODemux(MCT, assay = "HTO", positive.quantile = 0.99)



### Global classification results_look at dimensions


table(MCT$HTO_classification.global)

### Visualize enrichment for selected HTOs with ridge plots
### Group cells based on the max HTO signal

Idents(MCT) <- "HTO_maxID"




### Visualize pairs of HTO signals to confirm mutual exclusivity in singlets

DefaultAssay(MCT) <- "HTO"




### Compare number of UMIs for singlets, doublets and negative cells

Idents(MCT) <- "HTO_classification.global"



###Generate a two dimensional tSNE embedding for HTOs. Here we are grouping cells by singlets and doublets ###for simplicity.
### First, we will remove negative cells from the object

MCT <- subset(MCT, idents = "Negative", invert = TRUE)




### Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = MCT, assay = "HTO"))))


### Calculate tSNE embeddings with a distance matrix
MCT <- RunTSNE(MCT, distance.matrix = hto.dist.mtx, perplexity = 100)


###Plot_HTO_classification

pdf("HTO_classification.DimPlot.pdf",width=10,height=10,paper='special')
DimPlot(MCT, reduction = "tsne")
dev.off()


### To increase the efficiency of plotting, you can subsample cells using the num.cells argument

pdf("HTO.Heatmap.pdf",width=10,height=10,paper='special')
HTOHeatmap(MCT, assay = "HTO")
dev.off()



### Cluster and visualize cells using the usualHTOHeatmap(MCT.singlet, assay = "HTO") scRNA-seq workflow, 
### and examine for the potential presence of batch effects.
### Extract the singlets

MCT.singlet <- subset(MCT, idents = "Singlet")

# The [[ operator can add columns to object metadata. 
#This is a great place to stash QC stats

MCT.singlet[["percent.mt"]] <- PercentageFeatureSet(MCT.singlet, pattern = "^mt-")



### removing unwanted cells from the dataset

MCT.singlet <- subset(MCT.singlet, subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 8000 & percent.mt < 20)
Idents(MCT.singlet) <- MCT$HTO_maxID





###Alternatively set HTO_classification_group as the Idents

MCT.singlet@meta.data$HTO_classification_num <- 
  recode(MCT.singlet@meta.data$HTO_classification, 
         "M-HTO-1" = 1, "M-HTO-2" = 2, "M-HTO-3" = 3, "M-HTO-4" = 4, 
         "M-HTO-5" = 5, "M-HTO-6" = 6)


MCT.singlet@meta.data$HTO_classification_group <- 
  recode(MCT.singlet@meta.data$HTO_classification_num, 
         "1" = "Control", "2" = "Control", "3" = "Control", 
         "4" = "MCT_treatment", "5" = "MCT_treatment", "6" = "MCT_treatment")

Idents(MCT.singlet) <- MCT.singlet@meta.data$HTO_classification_group
#Idents(MCT.singlet_small) <- MCT.singlet@meta.data$DCid
Idents(MCT.singlet) #test if it works
MCT.singlet_small  <- subset(MCT.singlet, downsample = 2356)

MCT.singlet_small[["htoclusterid"]] <- Idents(MCT.singlet_small)

#Idents(MCT.singlet_small) <- MCT.singlet_small@meta.data$seurat_clusters


DefaultAssay(MCT.singlet_small) <- "RNA"


### Normalization
### run sctransform - did not run sct


MCT.singlet_small <- SCTransform(MCT.singlet_small, vars.to.regress = "percent.mt", 
                                 verbose = TRUE)

##alternatively - yes



MCT.singlet_small <- NormalizeData(MCT.singlet_small, 
                                   assay = "RNA", features = rownames(MCT.singlet_small))




### choose ~1k variable features


MCT.singlet_small <- FindVariableFeatures(MCT.singlet_small,
                                          assay = "RNA", features = rownames(MCT.singlet_small))



### standard scaling (no regression)


MCT.singlet_small <- ScaleData(MCT.singlet_small,
                               assay = "RNA", features = VariableFeatures(MCT.singlet_small))





### Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(MCT.singlet_small), 10)



# plot variable features with and without labels

pdf("Variable.Genes.pdf",width=15,height=4,paper='special')
plot1 <- VariableFeaturePlot(MCT.singlet_small)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


# Dimensionality reduction-Fig.2A


DefaultAssay(MCT.singlet_small) <- "RNA"

MCT.singlet_small <- RunPCA(MCT.singlet_small, features = VariableFeatures(object = MCT.singlet_small))

MCT.singlet_small <- FindNeighbors(MCT.singlet_small, dims = 1:10)
MCT.singlet_small <- FindClusters(MCT.singlet_small, resolution = 0.8)

MCT.singlet_small <- RunUMAP(MCT.singlet_small, dims = 1:10)  

pdf("Dimplot_umap_percentmt.pdf", width = 12, height = 6, paper = 'special')
DimPlot(MCT.singlet_small, split.by = "HTO_classification_group", reduction = "umap")
dev.off()


##Heatmap for different RNA clusters-Fig.2B
Idents(MCT.singlet_small) <- "seurat_clusters"

DefaultAssay(MCT.singlet_small) <- "RNA"

MCT.markers <- FindAllMarkers(MCT.singlet_small, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MCT.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(MCT.markers, file = "RNA_cluster.markers.csv")
DoHeatmap(MCT.singlet_small, features = unique(MCT.markers$gene))





###### Figure 2C
DimPlot(MCT.singlet_small, reduction = "umap")

#subset for different assays
Idents(MCT.singlet_small) <- "seurat_clusters"
Dynosubclusters <- subset(MCT.singlet_small, idents = c("0", "1", "3", "4", "7", "8", "10", "12", "14"))

Velocyto_clusters <- subset(MCT.singlet_small, idents = c("0", "1", "3", "4", "7", "8", "10", "12", "14"))

save.image("MCT.singlet_small.rda")


####for GEOsubmission
write.csv(MCT.singlet_small@assays$RNA@data, file = "MCT.singlet_small_RNA.csv")
write.csv(MCT.singlet_small@assays$CITE@data, file = "MCT.singlet_small_CITE.csv")
write.csv(MCT.singlet_small@meta.data, file = "MCT.singlet_small_metadata.csv")




##########


DefaultAssay(MCT.singlet) <- "ADT"


plotimmunecells <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD45", 
                                  feature2 = "Ptprc")
ImmunecellsCD45 <- CellSelector(plot = plotimmunecells)

Idents(MCT.singlet, cells = ImmunecellsCD45) <- "CD45+"



# Subset the CD45+ immune cells
Immunecells_singlets <- subset(MCT.singlet, idents = "CD45+")

##Subset Broad immune cells from Immunecells_singlets subset
plotmyeloid <- FeatureScatter(Immunecells_singlets, feature1 = "CITE-CD45", feature2 = "CITE-CD11b")
CD11bofCD45 <- CellSelector(plot = plotmyeloid, cells = ImmunecellsCD45)


plotlymphoid <- FeatureScatter(Immunecells_singlets, feature1 = "CITE-CD45", feature2 = "CITE-CD3")
CD3ofCD45 <- CellSelector(plot = plotlymphoid, cells = ImmunecellsCD45)

plotlymphoid <- FeatureScatter(Immunecells_singlets, feature1 = "CITE-CD45", feature2 = "CITE-CD11c")
CD11cofCD45 <- CellSelector(plot = plotlymphoid, cells = ImmunecellsCD45)


Idents(Immunecells_singlets, cells = CD3ofCD45) <- "T cells"
Idents(Immunecells_singlets, cells = CD11bofCD45) <- "Myeloid"
Idents(Immunecells_singlets, cells = CD11cofCD45) <- "Dendritic cells"


Lymphoid_cellsubset <- subset(Immunecells_singlets_small, cells = CD3ofCD45)
Dendritc_cellsubset <- subset(Immunecells_singlets_small, cells = CD11cofCD45)
Myeloid_cellsubset <- subset(Immunecells_singlets_small, cells = CD11bofCD45)

##########################
### plot Figure 2 counts for different cell types ####

ggplot(Lymphoid_cellsubset@meta.data,
       aes(x=Lymphoid_cellsubset@meta.data$HTO_classification_group)) +
  geom_bar()


ggplot(Dendritc_cellsubset@meta.data,
       aes(x=Dendritc_cellsubset@meta.data$HTO_classification_group)) +
  geom_bar()



ggplot(Myeloid_cellsubset@meta.data,
       aes(x=Myeloid_cellsubset@meta.data$HTO_classification_group)) +
  geom_bar()


##subseting TADcs and TAMs

DefaultAssay(MCT.singlet_small) <- "ADT"

#Broz et.al.,

#remove monocytes

Immunecells_TADsTAMs <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD45", feature2 = "Ptprc")
CD45hptprch <- CellSelector(plot = Immunecells_TADsTAMs)


plottoremovemonocytes <- FeatureScatter(MCT.singlet, feature1 = "CITE-Ly6C", feature2 = "CITE-CD11b", cells = CD45hptprch)
CD11bhlyc6low <- CellSelector(plot = plottoremovemonocytes)


plotMHC2highDcs <- FeatureScatter(MCT.singlet, feature1 = "CITE-I-A-I-E", feature2 = "CITE-CD11c", cells = CD11bhlyc6low)
CD11bMHC2 <- CellSelector(plot = plotMHC2highDcs)

#distinguish dendritic cells from macrophages
plotCD24high <- FeatureScatter(MCT.singlet, feature1 = "CITE-F4-80", feature2 = "CITE-CD24", cells = CD11bMHC2)
CD24Dendritic <- CellSelector(plot = plotCD24high)

plot103DCs <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD103", feature2 = "CITE-CD11b", cells = CD24Dendritic)
CD103Dendritic <- CellSelector(plot = plot103DCs)

plot11bDCs <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD103", feature2 = "CITE-CD11b", cells = CD24Dendritic)
Dendritic11b <- CellSelector(plot = plot11bDCs)

#TAMs

plotTAM <- FeatureScatter(MCT.singlet, feature1 = "CITE-F4-80", feature2 = "CITE-CD24", cells = CD11bMHC2)
TAM <- CellSelector(plot = plotTAM)

plot11cTAM <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD11c", feature2 = "CITE-CD11b", cells = TAM)
TAM11c <- CellSelector(plot = plot11cTAM)

plot11bTAM <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD11c", feature2 = "CITE-CD11b", cells = TAM)
TAM11b <- CellSelector(plot = plot11bTAM)

# Laoui et al., monocyte derived DCs.
plotmonocytic <- FeatureScatter(MCT.singlet, feature1 = "CITE-Ly6C", feature2 = "CITE-CD64", cells = CD11bMHC2)
CD64hly6ch <- CellSelector(plot = plotmonocytic)


plotmonDC <- FeatureScatter(MCT.singlet, feature1 = "CITE-CD11b", feature2 = "CITE-CD24", cells = CD64hly6ch)#typo
monDC <- CellSelector(plot = plotmonDC)

Idents(MCT.singlet_small, cells = monDC) <- "mDC"
Idents(MCT.singlet_small, cells = Dendritic11b) <- "CD11bDcs"
Idents(MCT.singlet_small, cells = TAM11b) <- "TAM11b"
Idents(MCT.singlet_small, cells = TAM11c) <- "TAM11c"
Idents(MCT.singlet_small, cells = CD103Dendritic) <- "CD103Dcs"

MCT.singlet_small[["TADCs_TAMs"]] <- Idents(MCT.singlet_small)

Idents(MCT.singlet_small) <- "TADCs_TAMs"

#plot Dimplot for Figure2_PanelB right
pdf("Dimplot_umap_TADCs_TAMs.pdf", width = 7, height = 6, paper = 'special')
DimPlot(MCT.singlet_small, reduction = "umap")
dev.off()


MonocyticDC.subset <- subset(MCT.singlet_small, idents = "mDC")
CD11bDC.subset <- subset(MCT.singlet_small, idents = "CD11bDcs")
CD103DC.subset <- subset(MCT.singlet_small, idents = "CD103Dcs")
TAM11b.subset <- subset(MCT.singlet_small, idents = "TAM11b")
TAM11c.subset <- subset (MCT.singlet_small, idents = "TAM11c")


#plot bar plots for different TADCs and TAMs - Figure2 counts ###############

pdf("mDC_counts.pdf", width = 5, height = 5, paper = 'special')
ggplot(MonocyticDC.subset@meta.data,
       aes(x=MonocyticDC.subset@meta.data$HTO_classification_group)) +
  geom_bar()
dev.off()


pdf("CD11bDcs_counts.pdf", width = 5, height = 5, paper = 'special')
ggplot(CD11bDC.subset@meta.data,
       aes(x=CD11bDC.subset@meta.data$HTO_classification_group)) +
  geom_bar()
dev.off()


pdf("CD103Dcs_counts.pdf", width = 5, height = 5, paper = 'special')
ggplot(CD103DC.subset@meta.data,
       aes(x=CD103DC.subset@meta.data$HTO_classification_group)) +
  geom_bar()
dev.off()


pdf("TAM11b_counts.pdf", width = 5, height = 5, paper = 'special')
ggplot(TAM11b.subset@meta.data,
       aes(x=TAM11b.subset@meta.data$HTO_classification_group)) +
  geom_bar()
dev.off()

pdf("TAM11c_counts.pdf", width = 5, height = 5, paper = 'special')
ggplot(TAM11c.subset@meta.data,
       aes(x=TAM11c.subset@meta.data$HTO_classification_group)) +
  geom_bar()
dev.off()


###### TADS TAMS RNA - DEG ###################### Need the load(file= "MCT.singlet_small_analysis_manuscript.rda")

TADsTAMs <- subset(MCT.singlet_small, idents = c("mDC","CD103Dcs", "TAM11c", "TAM11b", "CD11bDcs"))
DefaultAssay(TADsTAMs) <- "RNA"


TADsTAMs.markers <- FindAllMarkers(TADsTAMs, assay = "RNA")
write.csv(TADsTAMs.markers , file = "TADsTAMs.markers.csv")

DoHeatmap(TADsTAMs, features = unique(TADsTAMs.markers$gene)) ##### Figure S2B #####

Idents(TADsTAMs) <- "HTO_classification_group"
TADsTAMs.markersMCTvsCtrl <- FindMarkers(TADsTAMs, assay = "RNA",ident.1 = "MCT_treatment", ident.2 ="Control")
write.csv(TADsTAMs.markersMCTvsCtrl, file = "TADsTAMs.markersMCTvsCtrl.csv")

volcanodata <- read.csv("TADsTAMs.markersMCTvsCtrl.csv",
                        header=TRUE, stringsAsFactors=FALSE, row.names=1)


#plot Figure Supplementary 2C
ggplot(volcanodata, aes(x=avg_logFC, y=p_val_adj, cex=0.5)) + 
  
  geom_point(shape=19, size=3, color=ifelse(volcanodata$p_val_adj>20,
                                            'red', 'black')) +
  ylim(1,-1)

#to know the names of genes
ggplotly(ggplot(volcanodata, aes(x=avg_logFC, y=p_val_adj, cex=0.5, text = row.names(volcanodata) )) + 
           geom_point(shape=21, size=3, color=ifelse(volcanodata$p_val_adj>pThresh,
                                                     'red', 'black')) + 
           ylim(1,-1) +
           theme(panel.background = element_rect(fill=NA), panel.grid.major = element_line(colour = "grey80")))


############Myeloid Clusters ############################

Idents(MCT.singlet_small) <- "seurat_clusters"
myeloid_clusters <- subset(MCT.singlet_small, idents = c("0", "1", "3", "4", "7", "8", "10", "12", "14"))


##Figure 2D#####
ggplot(myeloid_clusters@meta.data, aes(seurat_clusters, group = HTO_classification_group)) +
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes(label = scales::percent(..prop..),
                
                y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="seurat_clusters") +
  facet_grid(~HTO_classification_group) +
  scale_y_continuous(labels = scales::percent)


Idents(myeloid_clusters) <- "seurat_clusters"
DimPlot(myeloid_clusters, reduction = "umap", split.by = "HTO_classification_group") ##Figure 2C

Idents(myeloid_clusters) <- "HTO_classification_group"
MCT.Ctrlmyeloid_markers <- FindMarkers(myeloid_clusters, ident.1 = "MCT_treatment", ident.2 = "Control", only.pos = FALSE)
write.csv(MCT.Ctrlmyeloid_markers, file = "DGE_MCT.Ctrlmyeloid_markers.csv")

####Figure S2D
volcanodata.DEGMyeloid <- read.csv("DGE_MCT.Ctrlmyeloid_markers.csv",
                        header=TRUE, stringsAsFactors=FALSE, row.names=1)

volcanodata.DEGMyeloid$p.adj <- -1*log(volcanodata.DEGMyeloid$p_val_adj, 10)

#convert p value into log values

pThresh <- 20

#plot Figure Supplementary 2D
ggplot(volcanodata.DEGMyeloid, aes(x=avg_logFC, y=p.adj , cex=0.5)) + 
  
  geom_point(shape=19, size=3) 

#to know the names of genes
ggplotly(ggplot(volcanodata.DEGMyeloid, aes(x=avg_logFC, y=p.adj, cex=0.5, text = row.names(volcanodata.DEGMyeloid) )) + 
           geom_point(shape=21, size=3) + 
           theme(panel.background = element_rect(fill=NA), panel.grid.major = element_line(colour = "grey80")))



######
#subset for reclustering

load(file = "MCT.singlet_small.rda")
Idents(MCT.singlet_small) <- "seurat_clusters"
Zeroandfoursubclusters <- subset(MCT.singlet_small, idents = c("0", "4"))

Zeroandfoursubclusters <- RunPCA(Zeroandfoursubclusters, features = VariableFeatures(object = Zeroandfoursubclusters))

Zeroandfoursubclusters <- FindNeighbors(Zeroandfoursubclusters, dims = 1:10)
Zeroandfoursubclusters <- FindClusters(Zeroandfoursubclusters, resolution = 0.8)

Zeroandfoursubclusters <- RunUMAP(Zeroandfoursubclusters, dims = 1:10)

DimPlot(Zeroandfoursubclusters, split.by = "HTO_classification_group", reduction = "umap")#Fig2E



Idents(Zeroandfoursubclusters) <- "seurat_clusters"

DefaultAssay(Zeroandfoursubclusters) <- "RNA"

Zeroandfoursubclusters.markers <- FindAllMarkers(Zeroandfoursubclusters,
                                                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Zeroandfoursubclusters.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(Zeroandfoursubclusters.markers, file = "ZeroandfoursubclustersRNA_cluster.markers.csv")

DoHeatmap(Zeroandfoursubclusters)

#####Aug112021
#1) sort based on fold change cut off ABS(0.7). 
#New version-Fig.S2F
NewZeroandfoursubclusters.markers <- FindAllMarkers(Zeroandfoursubclusters,
                                                 only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.7)
NewZeroandfoursubclusters.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
write.csv(NewZeroandfoursubclusters.markers, file = "NewZeroandfoursubclusters.markers.csv")

#DoHeatmap(Zeroandfoursubclusters, features = NewZeroandfoursubclusters.markers$gene)

#DoHeatmap(Zeroandfoursubclusters, features = NewZeroandfoursubclusters.markers$gene, slot = "data")

#DoHeatmap(Zeroandfoursubclusters, features = NewZeroandfoursubclusters.markers$gene, 
#          slot = "Zeroandfoursubclusters@assays$RNA")
DoHeatmap(Zeroandfoursubclusters, features = 
            c("Nudt4", "Cd274", "Mafk", "Jun", "Prok2", "Ddit3", "Tra2a", "Gla","Ptp4a1","Cpne2","Sde2","Ciart","Retnlg","Rhov","Ankrd33b",
            "Atf3", "Dedd2","Vps37b","Kpna4","Wfdc17","Suco","S100a9", "Slc31a2","1700017B05Rik", "1810058I24Rik", "Rab9", "Atp6v1c1",
             "Tgif1", "Id2", "F10", "Bnip3l", "Itpr2", "Card19", "Ero1l", "Echdc3", "Eif5","Ifrd1", "Ccl3", "Bri3", "Gadd45b",
             "Ccdc126", "Cstb", "Hilpda", "Fabp5","mt-Co1", "Mfge8", "Lyz2", "H2-Aa", "Ccl7", "H2-Ab1", "Cd74", "Mmp13","Spp1",
             "Cxcl1", "Prdx1", "Ccl2","Rpl32", "Rpsa","Ctss","S100a10","Tnfrsf9","Rpl13","Psap","Lgals1","Mt1", "Ccl9","Msr1",
            "Mmp12","Ctsl","Vim","Tmsb10","Arg1","Il7r","Lcn2","Stfa2l1","Lsp1","Cytip","Cdk2ap2","Retnlg","Stk17b","Selplg",
             "Fgl2", "Tnfaip2","Wfdc21","Mrgpra2a","Trem1","Gmfg","Lst1","G0s2","Gm5150","Ifitm1","Mrgpra2b","Csf3r","S100a8",
            "Cxcr2", "S100a9","Hdc", "S100a11", "Ifitm2","Pglyrp1","Lrg1","Hp","H2-Aa","H2-Ab1","Cd74", "Spp1", "Mt1", "Ccl12",
            "Apoe", "Lyz2","Fn1","Rps20","Rps8","Cst3","Rpl3","Prdx1","Grn","Rpl10a","C1qa", "Ctss","Lgals1","C1qb", "Trem2",
                                               "Ms4a7","Cd81","Ctsc","C1qc","Pmepa1","Ms4a6c","Hexb","Cxcl16","Gatm","Cxcl10","Sp100","Parp14","G0s2", "Ddx60", "Irf7",
                                               "Il1r2", "Gm13822","Plac8","Slfn1","Ifit2","Usp18","Rtp4","Ifit3","Ifit3b","Trim30a","Mxd1","1600014C10Rik","Oasl2",
                                               "Oasl1", "Ifitm3", "Ifi47","Slfn5","Isg20","Gbp2","Cmpk2","Slfn4","Isg15","Rsad2","Ifit1","Ccl4","Pts","Il1rn", "BC100530",
                                               "Ccl6", "Upp1","Stfa2l1","Cxcl3","Cxcl2","F10","Ccl3","Gm5483","Gsn","Il1b","Marcksl1","Bcl2a1b","Ltb","Gpr84","Clec4n",
                                               "Clec4e","Junb","Icam1", "Mmp9", "Il23a"), slot = "data") +
          scale_fill_gradientn(colors = c("white", "red"))#Final version Fig.2SF


####Fig.2E####

freq_table <- prop.table(x = table(Zeroandfoursubclusters@meta.data$seurat_clusters,
                                   Zeroandfoursubclusters@meta.data$HTO_classification_group),margin = 2)

library("scales")

coloridentities <- levels(Zeroandfoursubclusters@meta.data$seurat_clusters)

my_color_palette <- hue_pal()(length(coloridentities))

barplot(height = freq_table, col = my_color_palette)


##vlnplots Fig.2 and Fig.S2:

Idents(Zeroandfoursubclusters) <- "seurat_clusters"

VlnPlot(Zeroandfoursubclusters, features = "Axl", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()



VlnPlot(Zeroandfoursubclusters, features = "Cd72", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Zeroandfoursubclusters, features = "Cx3cr1", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Zeroandfoursubclusters, features = "Ddit3", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

VlnPlot(Zeroandfoursubclusters, features = "Irf8", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Zeroandfoursubclusters, features = "Mmp9", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Zeroandfoursubclusters, features = "Il23a", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Zeroandfoursubclusters, features = "Gpr84", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

VlnPlot(Zeroandfoursubclusters, features = "Atf3", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Zeroandfoursubclusters, features = "Cd274", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

VlnPlot(Zeroandfoursubclusters, features = "Hilpda", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()


VlnPlot(Zeroandfoursubclusters, features = "S100a9", group.by = "seurat_clusters") +
  stat_summary(fun = "mean", geom = "point",shape=23, size=4, fill="red") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means()+
  theme_classic()

save.image("Figure3.SeuratRNArelatedpanels.rda")
save.image("Figure3.CITErelatedpanels.rda")
