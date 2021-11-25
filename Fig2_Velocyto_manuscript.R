##### velocity analysis 

library(dplyr)
library(Seurat)
library(ggplot2)
library(umap)
library(velocyto.R)


load("MCT.singlet_small_analysis_manuscript.rda")

#subclustered: Velocyto_clusters <- subset(MCT.singlet_small, idents = c("0", "1", "3", "4", "7", "8", "10", "12", "14"))
#check for Velocyto_clusters
Velocyto_clusters


MCT_data <- Read10X(data = "filtered_feature_bc_matrix")

MCT <- CreateSeuratObject(counts = MCT_data$"Gene Expression")
MCT$CITE <- CreateAssayObject(counts = 
                                MCT_data$"Antibody Capture")
TotalA.features = read.csv(file = "TotalA_Ref.csv")

HTOs <- c("M-HTO-1","M-HTO-2","M-HTO-3","M-HTO-4","M-HTO-5","M-HTO-6")

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
### run sctransform


MCT.singlet_small <- SCTransform(MCT.singlet_small, vars.to.regress = "percent.mt", 
                                 verbose = TRUE)

##alternatively



MCT.singlet_small <- NormalizeData(MCT.singlet_small, 
                                   assay = "RNA", features = rownames(MCT.singlet_small))




### choose ~1k variable features


MCT.singlet_small <- FindVariableFeatures(MCT.singlet_small,
                                          assay = "RNA", features = rownames(MCT.singlet_small))



### standard scaling (no regression)


MCT.singlet_small <- ScaleData(MCT.singlet_small,
                               assay = "RNA", features = VariableFeatures(MCT.singlet_small))


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


MCT.singlet_small <- NormalizeData(MCT.singlet_small, 
                                   assay = "RNA", features = rownames(MCT.singlet_small))



### choose ~1k variable features


MCT.singlet_small <- FindVariableFeatures(MCT.singlet_small,
                                          assay = "RNA", features = rownames(MCT.singlet_small))



### standard scaling (no regression)


MCT.singlet_small <- ScaleData(MCT.singlet_small,
                               assay = "RNA", features = VariableFeatures(MCT.singlet_small))


# Dimensionality reduction


DefaultAssay(MCT.singlet_small) <- "RNA"

MCT.singlet_small <- RunPCA(MCT.singlet_small, features = VariableFeatures(object = MCT.singlet_small))

MCT.singlet_small <- FindNeighbors(MCT.singlet_small, dims = 1:10)
MCT.singlet_small <- FindClusters(MCT.singlet_small, resolution = 0.8)

MCT.singlet_small <- RunUMAP(MCT.singlet_small, dims = 1:10)

#sub classify myeloid cells

Velocyto_clusters <- subset(MCT.singlet_small, idents = c("0", "1", "3", "4", "7", "8", "10", "12", "14"))



#visualize umap
pdf("DimPlot.UMAP_split.pdf",width=12,height=6,paper='special')
DimPlot(MCT.singlet_small, reduction = "umap", split.by = "HTO_classification_group")
dev.off()

#load loom file

ldat.Velocyto_clusters <- read.loom.matrices("Reanalysis_Bhavana.loom")

#Gather the spliced and unspliced estimates and rename to match Seurat

emat.Velocyto_clusters <- ldat.Velocyto_clusters$spliced
head(colnames(emat.Velocyto_clusters))


# matching names "possorted_genome_bam_3A1T6:AACCCAACATCATCTTx" this is ori name. the goal is to extract only the barcode

colnames(emat.Velocyto_clusters) <- paste(substring(colnames(emat.Velocyto_clusters),20,35),sep="") 
head(colnames(emat.Velocyto_clusters))  # should show: "AACCCAACATCATCTT" ...


nmat.Velocyto_clusters <- ldat.Velocyto_clusters$unspliced
colnames(nmat.Velocyto_clusters) <- paste(substring(colnames(nmat.Velocyto_clusters),20,35),sep="") 

#Spliced expression magnitude distribution across genes:

pdf("ematVelocyto_clusters.pdf",width=6,height=6,paper='special')
hist(log10(colSums(emat.Velocyto_clusters)+1),col='wheat',xlab='colSums(emat.Velocyto_clusters) + 1',main='colSums(emat.Velocyto_clusters)')
dev.off()

pdf("namtVelocyto_clusters.pdf",width=6,height=6,paper='special')
hist(log10(colSums(nmat.Velocyto_clusters)+1),col='wheat',xlab='colSums(nmat.Velocyto_clusters) + 1',main='colSums(nmat.Velocyto_clusters)')
dev.off()

# take cluster label from pmbc object

cluster.labelVelocyto_clusters <- as.factor(as.numeric(Velocyto_clusters@meta.data$seurat_clusters))
names(cluster.labelVelocyto_clusters)  <- row.names(Velocyto_clusters@meta.data)

# filtering 


emat.Velocyto_clusters <- filter.genes.by.cluster.expression(emat.Velocyto_clusters, 
                                                    cluster.labelVelocyto_clusters, min.max.cluster.average = 0.05)
nmat.Velocyto_clusters <- filter.genes.by.cluster.expression(nmat.Velocyto_clusters, cluster.labelVelocyto_clusters, min.max.cluster.average = 0.01)
length(intersect(rownames(emat.Velocyto_clusters),rownames(nmat.Velocyto_clusters)))




# UMAP plot


pdf("DimPlot_Velocyto_clusters.UMAP.pdf",width=12,height=6,paper='special')
DimPlot(Velocyto_clusters, reduction = "umap", split.by = "HTO_classification_group")
dev.off()



# extract embedding 

embVelocyto_clusters <- Velocyto_clusters@reductions$umap@cell.embeddings 



# Estimate the cell-cell distances
cell.distVelocyto_clusters <- as.dist(1-armaCor(t(embVelocyto_clusters)))

# get tSNE colors


gg <- DimPlot(Velocyto_clusters, reduction = "umap", split.by = "HTO_classification_group")
a <- as.data.frame(ggplot_build(gg)$data)
colors <- as.list(a$colour)
names(colors) <- rownames(embVelocyto_clusters)



# Estimate RNA velocity (using gene-relative model with k=50 cell kNN pooling 
#and using top/bottom 2% quantiles for gamma fit): perform gamma 
#fit on a top/bottom quantiles of expression magnitudes


fit.quantile <- 0.02
rvel.cd.Velocyto_clusters <- gene.relative.velocity.estimates(emat.Velocyto_clusters,nmat.Velocyto_clusters,
                                                     deltaT=1, kCells=50, cell.dist=cell.dist,
                                                     fit.quantile=fit.quantile)

# plot velocity embedding
pdf("DimPlot.Velocyto_clusters.UMAP.velocity.pdf",width=6,height=6,paper='special')
p1 <- show.velocity.on.embedding.cor(embVelocyto_clusters,rvel.cd.Velocyto_clusters,n=200,scale="log",
                                     cell.colors=ac(colors,alpha=0.5),
                                     cex=0.8,arrow.scale=6,show.grid.flow=T,
                                     min.grid.cell.mass=1.0,grid.n=60,arrow.lwd=0.7,
                                     do.par=T,cell.border.alpha = 0.1, main="Cell Velocity") 
dev.off()

   


pdf("DimPlot.Velocyto_clusters_UMAP_sqrt.velocity.pdf",width=5.5,height=6,paper='special')
p2 <- show.velocity.on.embedding.cor(embVelocyto_clusters,rvel.cd.Velocyto_clusters, n=200, scale='sqrt', 
                                     cell.colors=ac(colors,alpha=0.8), cex=0.5, arrow.scale=5, 
                                     show.grid.flow=T, min.grid.cell.mass=0.5, grid.n=50, arrow.lwd=0.8, 
                                     do.par=F,  cell.border.alpha = 0.02,
                                     main="Velocyto_clusters.velocity") 

dev.off()   




# pca velocity
pdf("DimPlot.PCA.velocity.pdf",width=12,height=12,paper='special')
pca.velocity.plot(rvel.cd,nPcs=5,plot.cols=2,cell.colors=ac(colors,alpha=0.5), arrow.scale = 0.5, arrow.lwd = 0.5,
                  cex=1,pcount=0.1,pc.multipliers=c(1,-1,-1,-1,-1))

dev.off()             





# trajectory
pdf("DimPlot.Velocyto_clusters_trajectory.pdf",width=12,height=12,paper='special')
x <- show.velocity.on.embedding.eu(embVelocyto_clusters,rvel.cd.Velocyto_clusters,n=100,scale="sqrt",cell.colors=ac(colors,alpha=0.5)
                                   ,nPcs=8,sigma=2.5,show.trajectories=TRUE,diffusion.steps=20,
                                   n.trajectory.clusters=6,ntop.trajectories=1,embedding.knn=T,
                                   control.for.neighborhood.density=TRUE) 

dev.off()   

