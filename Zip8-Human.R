library(Seurat)
library(tidyverse)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(ggpubr)

#load sample
Normal1.Epcam=read.csv("Normal-1.csv")
# Rename the samples on second line with gene symbol
rownames(Normal1.Epcam)=make.names(Normal1.Epcam[, 2], unique = T)
Normal1.Epcam<-Normal1.Epcam[,-c(1,2)]
rownames(Normal1.Epcam)<-gsub("_", "-", rownames(Normal1.Epcam))
#Generate Seurat object for further analysis
Normal1 <- CreateSeuratObject(counts = Normal1.Epcam)
#Setup gateout strategy for data
Normal1 <- subset(Normal1,  subset = nFeature_RNA > 2500 & nFeature_RNA < 9000& percent.mt < 8)
Normal1 <- NormalizeData(Normal1, normalization.method = "LogNormalize", scale.factor = 1000000, verbose = FALSE)
Normal1 <- FindVariableFeatures(Normal1, selection.method = "vst", nfeatures = 4000)


#load sample
Normal2.Epcam=read.csv("Normal-2.csv")
# Rename the samples on second line with gene symbol
rownames(Normal2.Epcam)=make.names(Normal2.Epcam[, 2], unique = T)
Normal2.Epcam<-Normal2.Epcam[,-c(1,2)]
rownames(Normal2.Epcam)<-gsub("_", "-", rownames(Normal2.Epcam))
#Generate Seurat object for further analysis
Normal2 <- CreateSeuratObject(counts = Normal2.Epcam)
#Setup gateout strategy for data
Normal2 <- subset(Normal2,  subset = nFeature_RNA > 4000 & nFeature_RNA < 11000& percent.mt < 7)
Normal2 <- NormalizeData(Normal2, normalization.method = "LogNormalize", scale.factor = 1000000, verbose = FALSE)
Normal2 <- FindVariableFeatures(Normal2, selection.method = "vst", nfeatures = 4000)


#load sample
Normal3.Epcam=read.csv("Normal-3.csv")
# Rename the samples on second line with gene symbol
rownames(Normal3.Epcam)=make.names(Normal3.Epcam[, 2], unique = T)
Normal3.Epcam<-Normal3.Epcam[,-c(1,2)]
rownames(Normal3.Epcam)<-gsub("_", "-", rownames(Normal3.Epcam))
#Generate Seurat object for further analysis
Normal3 <- CreateSeuratObject(counts = Normal3.Epcam)
#Setup gateout strategy for data
Normal3 <- subset(Normal3,  subset = nFeature_RNA > 3000 & nFeature_RNA < 9000& percent.mt < 6)
Normal3 <- NormalizeData(Normal3, normalization.method = "LogNormalize", scale.factor = 1000000, verbose = FALSE)
Normal3 <- FindVariableFeatures(Normal3, selection.method = "vst", nfeatures = 4000)


#load sample
Normal4.Epcam=read.csv("Normal-4.csv")
# Rename the samples on second line with gene symbol
rownames(Normal4.Epcam)=make.names(Normal4.Epcam[, 2], unique = T)
Normal4.Epcam<-Normal4.Epcam[,-c(1,2)]
rownames(Normal4.Epcam)<-gsub("_", "-", rownames(Normal4.Epcam))
#Generate Seurat object for further analysis
Normal4 <- CreateSeuratObject(counts = Normal4.Epcam)
#Setup gateout strategy for data
Normal4 <- subset(Normal4,  subset = nFeature_RNA > 5000 & nFeature_RNA < 8000& percent.mt < 5)
Normal4 <- NormalizeData(Normal4, normalization.method = "LogNormalize", scale.factor = 1000000, verbose = FALSE)
Normal4 <- FindVariableFeatures(Normal4, selection.method = "vst", nfeatures = 4000)


#load sample
IPFA.Epcam=read.csv("IPF-1.csv")
# Rename the samples on second line with gene symbol
rownames(IPFA.Epcam)=make.names(IPFA.Epcam[, 2], unique = T)
IPFA.Epcam<-IPFA.Epcam[,-c(1,2)]
rownames(IPFA.Epcam)<-gsub("_", "-", rownames(IPFA.Epcam))
#Generate Seurat object for further analysis
IPFA <- CreateSeuratObject(counts = IPFA.Epcam)
#Setup gateout strategy for data
IPFA <- subset(IPFA,  subset = nFeature_RNA > 2000 & nFeature_RNA < 7500& percent.mt < 10)
IPFA <- NormalizeData(IPFA, normalization.method = "LogNormalize", scale.factor = 1000000, verbose = FALSE)
IPFA <- FindVariableFeatures(IPFA, selection.method = "vst", nfeatures = 4000)



#load sample
IPFB.Epcam=read.csv("IPF-2.csv")
# Rename the samples on second line with gene symbol
rownames(IPFB.Epcam)=make.names(IPFB.Epcam[, 2], unique = T)
IPFB.Epcam<-IPFB.Epcam[,-c(1,2)]
rownames(IPFB.Epcam)<-gsub("_", "-", rownames(IPFB.Epcam))
#Generate Seurat object for further analysis
IPFB <- CreateSeuratObject(counts = IPFB.Epcam)
#Setup gateout strategy for data
IPFB <- subset(IPFB,  subset = nFeature_RNA > 3000 & nFeature_RNA < 8500& percent.mt < 10)
IPFB <- NormalizeData(IPFB, verbose = FALSE)
IPFB <- FindVariableFeatures(IPFB, selection.method = "vst", nfeatures = 4000)


#load sample
IPFC.Epcam=read.csv("IPF-3.csv")
# Rename the samples on second line with gene symbol
rownames(IPFC.Epcam)=make.names(IPFC.Epcam[, 2], unique = T)
IPFC.Epcam<-IPFC.Epcam[,-c(1,2)]
rownames(IPFC.Epcam)<-gsub("_", "-", rownames(IPFC.Epcam))
#Generate Seurat object for further analysis
IPFC <- CreateSeuratObject(counts = IPFC.Epcam)
#Setup gateout strategy for data
IPFC <- subset(IPFC,  subset = nFeature_RNA > 4000 & nFeature_RNA < 8500& percent.mt < 8)
IPFC <- NormalizeData(IPFC, verbose = FALSE)
IPFC <- FindVariableFeatures(IPFC, selection.method = "vst", nfeatures = 4000)


#load sample
IPFD.Epcam=read.csv("IPF-4.csv")
# Rename the samples on second line with gene symbol
rownames(IPFD.Epcam)=make.names(IPFD.Epcam[, 2], unique = T)
IPFD.Epcam<-IPFD.Epcam[,-c(1,2)]
rownames(IPFD.Epcam)<-gsub("_", "-", rownames(IPFD.Epcam))
#Generate Seurat object for further analysis
IPFD <- CreateSeuratObject(counts = IPFD.Epcam)
#Setup gateout strategy for data
IPFD <- subset(IPFD,  subset = nFeature_RNA > 2500 & nFeature_RNA < 8000& percent.mt <10)
IPFD <- NormalizeData(IPFD, verbose = FALSE)
IPFD <- FindVariableFeatures(IPFD, selection.method = "vst", nfeatures = 4000)



#load sample
IPFE.Epcam=read.csv("IPF-5.csv")
# Rename the samples on second line with gene symbol
rownames(IPFE.Epcam)=make.names(IPFE.Epcam[, 2], unique = T)
IPFE.Epcam<-IPFE.Epcam[,-c(1,2)]
rownames(IPFE.Epcam)<-gsub("_", "-", rownames(IPFE.Epcam))
#Generate Seurat object for further analysis
IPFE <- CreateSeuratObject(counts = IPFE.Epcam)
#Setup gateout strategy for data
IPFE <- subset(IPFE,  subset = nFeature_RNA > 3000 & nFeature_RNA < 7500& percent.mt < 10)
IPFE <- NormalizeData(IPFE, verbose = FALSE)
IPFE <- FindVariableFeatures(IPFE, selection.method = "vst", nfeatures = 4000)



#load sample
IPFF.Epcam=read.csv("IPF-6.csv")
# Rename the samples on second line with gene symbol
rownames(IPFF.Epcam)=make.names(IPFF.Epcam[, 2], unique = T)
IPFF.Epcam<-IPFF.Epcam[,-c(1,2)]
rownames(IPFF.Epcam)<-gsub("_", "-", rownames(IPFF.Epcam))
#Generate Seurat object for further analysis
IPFF <- CreateSeuratObject(counts = IPFF.Epcam)
#Setup gateout strategy for data
IPFF <- subset(IPFF,  subset = nFeature_RNA > 2500 & nFeature_RNA < 7000& percent.mt < 10)
IPFF <- NormalizeData(IPFF, verbose = FALSE)
IPFF <- FindVariableFeatures(IPFF, selection.method = "vst", nfeatures = 4000)

Normal1$group="Health"
Normal2$group="Health"
Normal3$group="Health"
Normal4$group="Health"
Normal5$group="Health"
Normal6$group="Health"
IPFA$group="IPF"
IPFB$group="IPF"
IPFC$group="IPF"
IPFD$group="IPF"
IPFE$group="IPF"
IPFF$group="IPF"




#Find Anchors for data
AAAA.anchors <- FindIntegrationAnchors(object.list = list(Normal1, Normal2, Normal3, Normal4, Normal5, Normal6, IPFA,IPFB, IPFC, IPFD, IPFE), dims = 1:50)

#Combine data
AAAA.combined <- IntegrateData(anchorset = AAAA.anchors, dims = 1:50)



#Perform an integrated analysis
DefaultAssay(object = AAAA.combined) <- "integrated"
#Scale data, all genes need to be scaled for Doheatmap purposes
all.genes <- rownames(AAAA.combined)
# Run the standard workflow for visualization and clustering
AAAA.combined <- ScaleData(object = AAAA.combined, verbose = FALSE)
AAAA.combined <- RunPCA(object = AAAA.combined, npcs = 50, verbose = FALSE)

#After this, set dims to 25, but I think the first 15 clusters are true
#t-SNE and Clustering
AAAA.combined<- FindNeighbors(object = AAAA.combined, reduction = "pca", dims = 1:25)
AAAA.combined <- FindClusters(AAAA.combined, resolution = 1.0)
AAAA.combined <- RunUMAP(object = AAAA.combined, reduction = "pca", dims = 1:25)


#Visualization
p1 <- DimPlot(object = AAAA.combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(object = AAAA.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(object = AAAA.combined, reduction = "umap", split.by = "group")
DimPlot(object = AAAA.combined, reduction = "umap", label = TRUE)
DimPlot(object = AAAA.combined, reduction = "umap")

#Rename Basal and AEC2 CellScatter
AAAA.combined <- RenameIdents(AAAA.combined, `1` = "Basal", `7` = "Basal",`8` = "Basal",`9` = "Basal", `2` = "Club/Goblet", `0` = "AEC2", `3` = "AEC2", `4` = "AEC2", `6` = "Ciliated", `5` = "AEC1", `10` = "Cycling")
AAAA.combined <- RenameIdents(AAAA.combined, `12` = "PNEC")
AAAA.combined <- RenameIdents(AAAA.combined, `11` = "Basal")

#Rename cluster to be list in order
AAAA.combined <- RenameIdents(AAAA.combined, `Basal` = "0-Basal", `AEC2` = "1-AEC2",`AEC1` = "2-AEC1",`Club/Goblet` = "3-Club", `Ciliated` = "4-Ciliated", `PNEC` = "5-PNEC", `Cycling` = "6-Cycling", `CCSP+SPC+` = "7-Unknown")
DimPlot(AAAA.combined)


Normal1$group="Health1"
Normal2$group="Health2"
Normal3$group="Health3"
Normal4$group="Health4"
Normal5$group="Health5"
Normal6$group="Health6"
IPFA$group="IPF1"
IPFB$group="IPF2"
IPFC$group="IPF3"
IPFD$group="IPF4"
IPFE$group="IPF5"
IPFF$group="IPF6"

#Find Anchors for data
ABAA.anchors <- FindIntegrationAnchors(object.list = list(Normal1, Normal2, Normal3, Normal4, Normal5, Normal6, IPFA,IPFB, IPFC, IPFD, IPFE, IPFF), dims = 1:50)

#Combine data
ABAA.combined <- IntegrateData(anchorset = ABAA.anchors, dims = 1:50)

#Perform an integrated analysis
DefaultAssay(object = ABAA.combined) <- "integrated"
#Scale data, all genes need to be scaled for Doheatmap purposes
all.genes <- rownames(ABAA.combined)
# Run the standard workflow for visualization and clustering
ABAA.combined <- ScaleData(object = ABAA.combined, verbose = FALSE)
ABAA.combined <- RunPCA(object = ABAA.combined, npcs = 50, verbose = FALSE)

#After this, set dims to 25, but I think the first 15 clusters are true
#t-SNE and Clustering
ABAA.combined<- FindNeighbors(object = ABAA.combined, reduction = "pca", dims = 1:25)
ABAA.combined <- FindClusters(ABAA.combined, resolution = 1.0)
ABAA.combined <- RunUMAP(object = ABAA.combined, reduction = "pca", dims = 1:25)


#Visualization
p1 <- DimPlot(object = ABAA.combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(object = ABAA.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(object = ABAA.combined, reduction = "umap", split.by = "group")
DimPlot(object = ABAA.combined, reduction = "umap", label = TRUE)
DimPlot(object = ABAA.combined, reduction = "umap")

#Rename Basal and AEC2 CellScatter
ABAA.combined <- RenameIdents(ABAA.combined, `1` = "Basal", `7` = "Basal",`8` = "Basal",`9` = "Basal", `2` = "Club/Goblet", `0` = "AEC2", `3` = "AEC2", `4` = "AEC2", `6` = "Ciliated", `5` = "AEC1", `10` = "Cycling")
ABAA.combined <- RenameIdents(ABAA.combined, `12` = "PNEC")
ABAA.combined <- RenameIdents(ABAA.combined, `11` = "Basal")

#Rename cluster to be list in order
ABAA.combined <- RenameIdents(ABAA.combined, `Basal` = "0-Basal", `AEC2` = "1-AEC2",`AEC1` = "2-AEC1",`Club/Goblet` = "3-Club", `Ciliated` = "4-Ciliated", `PNEC` = "5-PNEC", `Cycling` = "6-Cycling", `CCSP+SPC+` = "7-Unknown")
DimPlot(ABAA.combined)

##Figure 1A
DimPlot(ABAA.combined,, reduction = "umap", group.by = "group")

##Figure 1B
DimPlot(AAAA.combined,, reduction = "umap", split.by = "group")

##Figure 1B
AEC2.cells <- subset(AAAA.combined, idents = "AEC2")
DefaultAssay(AEC2.cells) <- "RNA"

SIRT=read.csv("SIRT.csv", header=T)
colnames(SIRT)="gene"
DefaultAssay(AEC2.cell) <- "RNA"
a=as.data.frame(rownames(AEC2.cell))
colnames(a)='gene'
SIRTlist=intersect(SIRT, a, by="gene")
SIRTs=SIRTlist$gene
SIRTscore <- Matrix::colSums(AEC2.cell [SIRTs, ])
percent.SIRT <- Matrix::colSums(AEC2.cell [SIRTs, ])/Matrix::colSums(AEC2.cell)
AEC2.cell$ SIRTscore = SIRTscore
AEC2.cell$ percent.SIRT = percent.SIRT
FeaturePlot(AEC2.cell, features = "SIRTscore", reduction = "umap", cols= c("grey", "red"), min.cutoff = "q15", max.cutoff = "q80")
VlnPlot(AEC2.cell, features = "SIRTscore", pt.size =0.0 ,cols= c("red", "blue"),group.by = "group", y.max = 40)+ stat_summary(fun = mean, geom='errorbar',aes(ymax = ..y.., ymin = ..y..), width = 0.85, linetype = "dashed", size=1.5, col="black")+ FontSize(x.text = 30, y.text = 30, x.title = 30,y.title = 30, main = 45)+ NoLegend()








