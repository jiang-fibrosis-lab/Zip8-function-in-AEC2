library(Seurat)
library(tidyverse)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(hdf5r)



#load sample
Ctrl1.all=Read10X(data.dir="/Zip8 project/Control1/filtered_feature_bc_matrix/")
#Generate Seurat object for further analy
Ctrl1 <- CreateSeuratObject(counts = Ctrl1.all,project = "zip8", min.cells = 3, min.features = 200)
Ctrl1 <- RenameCells(object = Ctrl1, add.cell.id = "Ctrl1")
Ctrl1[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl1, pattern = "^mt-")
VlnPlot(object = Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = Ctrl1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Ctrl1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(object = Ctrl1, feature1 = "nFeature_RNA", feature2 = "percent.mt")

CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))

#Setup gateout strategy for data
Ctrl1 <- subset(Ctrl1,  subset = nFeature_RNA > 200 & percent.mt < 15& nCount_RNA< 50000)
Ctrl1 <- NormalizeData(Ctrl1, verbose = FALSE)
Ctrl1 <- FindVariableFeatures(Ctrl1, selection.method = "vst", nfeatures = 2000)

#load sample
Ctrl2.all=Read10X(data.dir="/Zip8 project/Control2/filtered_feature_bc_matrix/")
#Generate Seurat object for further analy
Ctrl2 <- CreateSeuratObject(counts = Ctrl2.all,project = "zip8", min.cells = 3, min.features = 200)
Ctrl2 <- RenameCells(object = Ctrl2, add.cell.id = "Ctrl2")
Ctrl2[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl2, pattern = "^mt-")
VlnPlot(object = Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = Ctrl2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Ctrl2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(object = Ctrl2, feature1 = "nFeature_RNA", feature2 = "percent.mt")

CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))

#Setup gateout strategy for data
Ctrl2 <- subset(Ctrl2,  subset = nFeature_RNA > 200 & percent.mt < 15& nCount_RNA< 50000)
Ctrl2 <- NormalizeData(Ctrl2, verbose = FALSE)
Ctrl2 <- FindVariableFeatures(Ctrl2, selection.method = "vst", nfeatures = 2000)


#load sample
Ctrl3.all=Read10X(data.dir="/Zip8 project/Control3/filtered_feature_bc_matrix/")
#Generate Seurat object for further analy
Ctrl3 <- CreateSeuratObject(counts = Ctrl3.all,project = "zip8", min.cells = 3, min.features = 200)
Ctrl3 <- RenameCells(object = Ctrl3, add.cell.id = "Ctrl3")
Ctrl3[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl3, pattern = "^mt-")
VlnPlot(object = Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = Ctrl3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Ctrl3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(object = Ctrl3, feature1 = "nFeature_RNA", feature2 = "percent.mt")

CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))

#Setup gateout strategy for data
Ctrl3 <- subset(Ctrl3,  subset = nFeature_RNA > 200 & percent.mt < 15& nCount_RNA< 50000)
Ctrl3 <- NormalizeData(Ctrl3, verbose = FALSE)
Ctrl3 <- FindVariableFeatures(Ctrl3, selection.method = "vst", nfeatures = 2000)


#load sample
Ctrl4.all=Read10X(data.dir="/Zip8 project/Control4/filtered_feature_bc_matrix/")
#Generate Seurat object for further analy
Ctrl4 <- CreateSeuratObject(counts = Ctrl4.all,project = "zip8", min.cells = 3, min.features = 200)
Ctrl4 <- RenameCells(object = Ctrl4, add.cell.id = "Ctrl4")
Ctrl4[["percent.mt"]] <- PercentageFeatureSet(object = Ctrl4, pattern = "^mt-")
VlnPlot(object = Ctrl4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = Ctrl4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Ctrl4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(object = Ctrl4, feature1 = "nFeature_RNA", feature2 = "percent.mt")

CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))

#Setup gateout strategy for data
Ctrl4 <- subset(Ctrl4,  subset = nFeature_RNA > 200 & percent.mt < 15& nCount_RNA< 50000)
Ctrl4 <- NormalizeData(Ctrl4, verbose = FALSE)
Ctrl4 <- FindVariableFeatures(Ctrl4, selection.method = "vst", nfeatures = 2000)

Zip81.all=Read10X(data.dir="/Zip8 project/Control1/filtered_feature_bc_matrix/")
#Generate Seurat object for further analy
Zip81 <- CreateSeuratObject(counts = Zip81.all,project = "zip8", min.cells = 3, min.features = 200)
Zip81 <- RenameCells(object = Zip81, add.cell.id = "Zip81")
Zip81[["percent.mt"]] <- PercentageFeatureSet(object = Zip81, pattern = "^mt-")
VlnPlot(object = Zip81, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = Zip81, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Zip81, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(object = Zip81, feature1 = "nFeature_RNA", feature2 = "percent.mt")

CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))

#Setup gateout strategy for data
Zip81 <- subset(Zip81,  subset = nFeature_RNA > 200 & percent.mt < 15& nCount_RNA< 50000)
Zip81 <- NormalizeData(Zip81, verbose = FALSE)
Zip81 <- FindVariableFeatures(Zip81, selection.method = "vst", nfeatures = 2000)


Zip82.all=Read10X(data.dir="/Zip8 project/Control1/filtered_feature_bc_matrix/")
#Generate Seurat object for further analy
Zip82 <- CreateSeuratObject(counts = Zip82.all,project = "zip8", min.cells = 3, min.features = 200)
Zip82 <- RenameCells(object = Zip82, add.cell.id = "Zip82")
Zip82[["percent.mt"]] <- PercentageFeatureSet(object = Zip82, pattern = "^mt-")
VlnPlot(object = Zip82, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(object = Zip82, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Zip82, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(object = Zip82, feature1 = "nFeature_RNA", feature2 = "percent.mt")

CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))

#Setup gateout strategy for data
Zip82 <- subset(Zip82,  subset = nFeature_RNA > 200 & percent.mt < 15& nCount_RNA< 50000)
Zip82 <- NormalizeData(Zip82, verbose = FALSE)
Zip82 <- FindVariableFeatures(Zip82, selection.method = "vst", nfeatures = 2000)


Ctrl1$sample="Ctrl1"
Ctrl2$sample="Ctrl2"
Ctrl3$sample="Ctrl3"
Ctrl4$sample="Ctrl4"
Zip81$sample="Zip81"
Zip82$sample="Zip82"

Ctrl1$sample="Ctrl"
Ctrl2$sample="Ctrl"
Ctrl3$sample="Ctrl"
Ctrl4$sample="Ctrl"
Zip81$sample="Zip8"
Zip82$sample="Zip8"

Zip8.anchors <- FindIntegrationAnchors(object.list = list(Ctrl1, Ctrl2, Ctrl3, Ctrl4, Zip81, Zip82), anchor.features = 2000,  dims = 1:30)

save(Zip8.anchors, file="Zip8anchor.rds")
Zip8.combined <- IntegrateData(anchorset = Zip8.anchors, dims = 1:30)

#Perform an integrated analysis
DefaultAssay(object = Zip8.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
Zip8.combined <- ScaleData(object = Zip8.combined, verbose = FALSE)
Zip8.combined <- RunPCA(object = Zip8.combined, npcs = 50, verbose = FALSE)

Zip8.combined <- JackStraw(Zip8.combined, num.replicate = 100)
Zip8.combined <- ScoreJackStraw(Zip8.combined, dims = 1:20)
JackStrawPlot(Zip8.combined, dims = 1:20)
# Method 2, default dims are 20, can set dims with ndims
ElbowPlot(Zip8.combined, ndims = 30)

# Run the standard workflow for visualization and clustering
Zip8.combined<- FindNeighbors(object = Zip8.combined, reduction = "pca", dims = 1:13)
Zip8.combined <- FindClusters(Zip8.combined, resolution = 0.2)
Zip8.combined <- RunUMAP(object = Zip8.combined, reduction = "pca",dims = 1:13)

DimPlot(object = Zip8.combined, reduction = "umap")
DimPlot(object = Zip8.combined, reduction = "umap", label=T)
DimPlot(object = Zip8.combined, reduction = "umap", split.by = "sample", ncol=3, label=T)

DefaultAssay(Zip8.combined) <- "RNA"                              
Zip8.combined.markers <- FindAllMarkers(Zip8.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.20)    

##Figure 5N
AEC2.cell=subset(Zip8.combined, idents="1")
DefaultAssay(AEC2.cell)="RNA"
VlnPlot(AEC2.cell, features = c("Chil1","H2-K1"), pt.size =0.0 ,cols= c("red", "blue"),group.by = "sample", ncol=6)
