library(Seurat)
library(tidyverse)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(hdf5r)

YD01.DATA=read.csv("YD0-1.csv")
# Rename the samples on second line with gene symbol
rownames(YD01.DATA)=make.names(YD01.DATA[, 2], unique = T)
YD01.DATA<-YD01.DATA[,-c(1,2)]
rownames(YD01.DATA)<-gsub("_", "-", rownames(YD01.DATA))
#Generate Seurat object for further analysis
YD01 <- CreateSeuratObject(counts = YD01.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD01[["percent.mt"]] <- PercentageFeatureSet(object = YD01, pattern = "^mt\\.")
VlnPlot(object = YD01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD01, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD01, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD01, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD01 <- subset(YD01,  subset = nFeature_RNA > 1000 & nFeature_RNA < 3000& percent.mt < 10)
YD01 <- NormalizeData(YD01, verbose = FALSE)
YD01 <- FindVariableFeatures(YD01, selection.method = "vst", nfeatures = 2000)

#load sample
YD02.DATA=read.csv("YD0-2.csv")
# Rename the samples on second line with gene symbol
rownames(YD02.DATA)=make.names(YD02.DATA[, 2], unique = T)
YD02.DATA<-YD02.DATA[,-c(1,2)]
rownames(YD02.DATA)<-gsub("_", "-", rownames(YD02.DATA))
#Generate Seurat object for further analysis
YD02 <- CreateSeuratObject(counts = YD02.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD02[["percent.mt"]] <- PercentageFeatureSet(object = YD02, pattern = "^mt\\.")
VlnPlot(object = YD02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD02, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD02, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD02, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD02 <- subset(YD02,  subset = nFeature_RNA > 500 & nFeature_RNA < 2500& percent.mt < 10)
YD02 <- NormalizeData(YD02, verbose = FALSE)
YD02 <- FindVariableFeatures(YD02, selection.method = "vst", nfeatures = 2000)

#load sample
YD03.DATA=read.csv("YD0-3.csv")
# Rename the samples on second line with gene symbol
rownames(YD03.DATA)=make.names(YD03.DATA[, 2], unique = T)
YD03.DATA<-YD03.DATA[,-c(1,2)]
rownames(YD03.DATA)<-gsub("_", "-", rownames(YD03.DATA))
#Generate Seurat object for further analysis
YD03 <- CreateSeuratObject(counts = YD03.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD03[["percent.mt"]] <- PercentageFeatureSet(object = YD03, pattern = "^mt\\.")
VlnPlot(object = YD03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD03, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD03, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD03, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD03 <- subset(YD03,  subset = nFeature_RNA > 1000 & nFeature_RNA < 2500& percent.mt < 10)
YD03 <- NormalizeData(YD03, verbose = FALSE)
YD03 <- FindVariableFeatures(YD03, selection.method = "vst", nfeatures = 2000)

#load sample
YD41.DATA=read.csv("YD4-1.csv")
# Rename the samples on second line with gene symbol
rownames(YD41.DATA)=make.names(YD41.DATA[, 2], unique = T)
YD41.DATA<-YD41.DATA[,-c(1,2)]
rownames(YD41.DATA)<-gsub("_", "-", rownames(YD41.DATA))
#Generate Seurat object for further analysis
YD41 <- CreateSeuratObject(counts = YD41.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD41[["percent.mt"]] <- PercentageFeatureSet(object = YD41, pattern = "^mt\\.")
VlnPlot(object = YD41, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD41, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD41, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD41, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD41 <- subset(YD41,  subset = nFeature_RNA > 1500 & nFeature_RNA < 5000& percent.mt < 10)
YD41 <- NormalizeData(YD41, verbose = FALSE)
YD41 <- FindVariableFeatures(YD41, selection.method = "vst", nfeatures = 2000)


#load sample
YD42.DATA=read.csv("YD4-2.csv")
# Rename the samples on second line with gene symbol
rownames(YD42.DATA)=make.names(YD42.DATA[, 2], unique = T)
YD42.DATA<-YD42.DATA[,-c(1,2)]
rownames(YD42.DATA)<-gsub("_", "-", rownames(YD42.DATA))
#Generate Seurat object for further analysis
YD42 <- CreateSeuratObject(counts = YD42.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD42[["percent.mt"]] <- PercentageFeatureSet(object = YD42, pattern = "^mt\\.")
VlnPlot(object = YD42, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD42, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD42, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD42, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD42 <- subset(YD42,  subset = nFeature_RNA > 1000 & nFeature_RNA < 5000& percent.mt < 10)
YD42 <- NormalizeData(YD42, verbose = FALSE)
YD42 <- FindVariableFeatures(YD42, selection.method = "vst", nfeatures = 2000)

#load sample
YD43.DATA=read.csv("YD4-3.csv")
# Rename the samples on second line with gene symbol
rownames(YD43.DATA)=make.names(YD43.DATA[, 2], unique = T)
YD43.DATA<-YD43.DATA[,-c(1,2)]
rownames(YD43.DATA)<-gsub("_", "-", rownames(YD43.DATA))
#Generate Seurat object for further analysis
YD43 <- CreateSeuratObject(counts = YD43.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD43[["percent.mt"]] <- PercentageFeatureSet(object = YD43, pattern = "^mt\\.")
VlnPlot(object = YD43, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD43, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD43, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD43, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD43 <- subset(YD43,  subset = nFeature_RNA > 1500 & nFeature_RNA < 5000& percent.mt < 8)
YD43 <- NormalizeData(YD43, verbose = FALSE)
YD43 <- FindVariableFeatures(YD43, selection.method = "vst", nfeatures = 2000)

#load sample
YD44.DATA=read.csv("YD4-4.csv")
# Rename the samples on second line with gene symbol
rownames(YD44.DATA)=make.names(YD44.DATA[, 2], unique = T)
YD44.DATA<-YD44.DATA[,-c(1,2)]
rownames(YD44.DATA)<-gsub("_", "-", rownames(YD44.DATA))
#Generate Seurat object for further analysis
YD44 <- CreateSeuratObject(counts = YD44.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD44[["percent.mt"]] <- PercentageFeatureSet(object = YD44, pattern = "^mt\\.")
VlnPlot(object = YD44, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD44, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD44, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD44, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD44 <- subset(YD44,  subset = nFeature_RNA > 2000 & nFeature_RNA < 5000& percent.mt < 8)
YD44 <- NormalizeData(YD44, verbose = FALSE)
YD44 <- FindVariableFeatures(YD44, selection.method = "vst", nfeatures = 2000)

#load sample
OD01.DATA=read.csv("OD0-1.csv")
# Rename the samples on second line with gene symbol
rownames(OD01.DATA)=make.names(OD01.DATA[, 2], unique = T)
OD01.DATA<-OD01.DATA[,-c(1,2)]
rownames(OD01.DATA)<-gsub("_", "-", rownames(OD01.DATA))
#Generate Seurat object for further analysis
OD01 <- CreateSeuratObject(counts = OD01.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD01[["percent.mt"]] <- PercentageFeatureSet(object = OD01, pattern = "^mt\\.")
VlnPlot(object = OD01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD01, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD01, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD01, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
OD01 <- subset(OD01,  subset = nFeature_RNA > 800 & nFeature_RNA < 3000& percent.mt < 10)
OD01 <- NormalizeData(OD01, verbose = FALSE)
OD01 <- FindVariableFeatures(OD01, selection.method = "vst", nfeatures = 2000)

#load sample
OD02.DATA=read.csv("OD0-2.csv")
# Rename the samples on second line with gene symbol
rownames(OD02.DATA)=make.names(OD02.DATA[, 2], unique = T)
OD02.DATA<-OD02.DATA[,-c(1,2)]
rownames(OD02.DATA)<-gsub("_", "-", rownames(OD02.DATA))
#Generate Seurat object for further analysis
OD02 <- CreateSeuratObject(counts = OD02.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD02[["percent.mt"]] <- PercentageFeatureSet(object = OD02, pattern = "^mt\\.")
VlnPlot(object = OD02, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD02, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD02, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD02, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
OD02 <- subset(OD02,  subset = nFeature_RNA > 700 & nFeature_RNA < 3000& percent.mt < 10)
OD02 <- NormalizeData(OD02, verbose = FALSE)
OD02 <- FindVariableFeatures(OD02, selection.method = "vst", nfeatures = 2000)

#load sample
OD03.DATA=read.csv("OD0-3.csv")
# Rename the samples on second line with gene symbol
rownames(OD03.DATA)=make.names(OD03.DATA[, 2], unique = T)
OD03.DATA<-OD03.DATA[,-c(1,2)]
rownames(OD03.DATA)<-gsub("_", "-", rownames(OD03.DATA))
#Generate Seurat object for further analysis
OD03 <- CreateSeuratObject(counts = OD03.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD03[["percent.mt"]] <- PercentageFeatureSet(object = OD03, pattern = "^mt\\.")
VlnPlot(object = OD03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD03, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD03, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD03, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
OD03 <- subset(OD03,  subset = nFeature_RNA > 700 & nFeature_RNA < 3000& percent.mt < 10)
OD03 <- NormalizeData(OD03, verbose = FALSE)
OD03 <- FindVariableFeatures(OD03, selection.method = "vst", nfeatures = 2000)

#load sample
OD41.DATA=read.csv("OD4-1.csv")
# Rename the samples on second line with gene symbol
rownames(OD41.DATA)=make.names(OD41.DATA[, 2], unique = T)
OD41.DATA<-OD41.DATA[,-c(1,2)]
rownames(OD41.DATA)<-gsub("_", "-", rownames(OD41.DATA))
#Generate Seurat object for further analysis
OD41 <- CreateSeuratObject(counts = OD41.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD41[["percent.mt"]] <- PercentageFeatureSet(object = OD41, pattern = "^mt\\.")
VlnPlot(object = OD41, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD41, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD41, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD41, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
OD41 <- subset(OD41,  subset = nFeature_RNA > 400 & nFeature_RNA < 2500& percent.mt < 10)
OD41 <- NormalizeData(OD41, verbose = FALSE)
OD41 <- FindVariableFeatures(OD41, selection.method = "vst", nfeatures = 2000)

#load sample
OD42.DATA=read.csv("OD4-2.csv")
# Rename the samples on second line with gene symbol
rownames(OD42.DATA)=make.names(OD42.DATA[, 2], unique = T)
OD42.DATA<-OD42.DATA[,-c(1,2)]
rownames(OD42.DATA)<-gsub("_", "-", rownames(OD42.DATA))
#Generate Seurat object for further analysis
OD42 <- CreateSeuratObject(counts = OD42.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD42[["percent.mt"]] <- PercentageFeatureSet(object = OD42, pattern = "^mt\\.")
VlnPlot(object = OD42, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD42, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD42, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD42, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
OD42 <- subset(OD42,  subset = nFeature_RNA > 500 & nFeature_RNA < 2500& percent.mt < 8)
OD42 <- NormalizeData(OD42, verbose = FALSE)
OD42 <- FindVariableFeatures(OD42, selection.method = "vst", nfeatures = 2000)

#load sample
OD43.DATA=read.csv("OD4-3.csv")
# Rename the samples on second line with gene symbol
rownames(OD43.DATA)=make.names(OD43.DATA[, 2], unique = T)
OD43.DATA<-OD43.DATA[,-c(1,2)]
rownames(OD43.DATA)<-gsub("_", "-", rownames(OD43.DATA))
#Generate Seurat object for further analysis
OD43 <- CreateSeuratObject(counts = OD43.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD43[["percent.mt"]] <- PercentageFeatureSet(object = OD43, pattern = "^mt\\.")
VlnPlot(object = OD43, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD43, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD43, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD43, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
OD43 <- subset(OD43,  subset = nFeature_RNA > 500 & nFeature_RNA < 2500& percent.mt < 10)
OD43 <- NormalizeData(OD43, verbose = FALSE)
OD43 <- FindVariableFeatures(OD43, selection.method = "vst", nfeatures = 2000)

#load sample
OD44.DATA=read.csv("OD4-4.csv")
# Rename the samples on second line with gene symbol
rownames(OD44.DATA)=make.names(OD44.DATA[, 2], unique = T)
OD44.DATA<-OD44.DATA[,-c(1,2)]
rownames(OD44.DATA)<-gsub("_", "-", rownames(OD44.DATA))
#Generate Seurat object for further analysis
OD44 <- CreateSeuratObject(counts = OD44.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD44[["percent.mt"]] <- PercentageFeatureSet(object = OD44, pattern = "^mt\\.")
VlnPlot(object = OD44, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD44, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD44, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD44, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
OD44 <- subset(OD44,  subset = nFeature_RNA > 700 & nFeature_RNA < 2500& percent.mt < 8)
OD44 <- NormalizeData(OD44, verbose = FALSE)
OD44 <- FindVariableFeatures(OD44, selection.method = "vst", nfeatures = 2000)

#load sample
YD141.DATA=Read10X_h5("YD14-1.h5")
rownames(YD141.DATA)<-gsub("_", "-", rownames(YD141.DATA))

#Generate Seurat object for further analysis
YD141 <- CreateSeuratObject(counts = YD141.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD141[["percent.mt"]] <- PercentageFeatureSet(object = YD141, pattern = "^mt-")
VlnPlot(object = YD141, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD141, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD141, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD141, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data
YD141 <- subset(YD141,  subset = nFeature_RNA > 2000 & nFeature_RNA < 6500& percent.mt < 20)
YD141 <- NormalizeData(YD141, verbose = FALSE)
YD141 <- FindVariableFeatures(YD141, selection.method = "vst", nfeatures = 2000)


#load sample
YD142.DATA=Read10X_h5("YD14-2.h5")
rownames(YD142.DATA)<-gsub("_", "-", rownames(YD142.DATA))
#Generate Seurat object for further analysis
YD142 <- CreateSeuratObject(counts = YD142.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD142[["percent.mt"]] <- PercentageFeatureSet(object = YD142, pattern = "^mt-")
VlnPlot(object = YD142, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD142, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD142, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD142, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data; From QC not a good sample
YD142 <- subset(YD142,  subset = nFeature_RNA > 1800 & nFeature_RNA < 5000& percent.mt < 20)
YD142 <- NormalizeData(YD142, verbose = FALSE)
YD142 <- FindVariableFeatures(YD142, selection.method = "vst", nfeatures = 2000)


#load sample
YD143.DATA=Read10X_h5("YD14-3.h5")
rownames(YD143.DATA)<-gsub("_", "-", rownames(YD143.DATA))
#Generate Seurat object for further analysis
YD143 <- CreateSeuratObject(counts = YD143.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
YD143[["percent.mt"]] <- PercentageFeatureSet(object = YD143, pattern = "^mt-")
VlnPlot(object = YD143, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(YD143, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(YD143, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(YD143, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data; From QC not a good sample
YD143 <- subset(YD143,  subset = nFeature_RNA > 1250 & nFeature_RNA < 7500 &nCount_RNA < 50000& percent.mt < 20)
YD143 <- NormalizeData(YD143, verbose = FALSE)
YD143 <- FindVariableFeatures(YD143, selection.method = "vst", nfeatures = 2000)

#load sample
OD141.DATA=Read10X_h5("OD14-1.h5")
rownames(OD141.DATA)<-gsub("_", "-", rownames(OD141.DATA))
#Generate Seurat object for further analysis
OD141 <- CreateSeuratObject(counts = OD141.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD141[["percent.mt"]] <- PercentageFeatureSet(object = OD141, pattern = "^mt-")
VlnPlot(object = OD141, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD141, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD141, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD141, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data; From QC not a good sample
OD141 <- subset(OD141,  subset = nFeature_RNA > 1800 & nFeature_RNA < 7500 &nCount_RNA < 50000 & percent.mt < 20)
OD141 <- NormalizeData(OD141, verbose = FALSE)
OD141 <- FindVariableFeatures(OD141, selection.method = "vst", nfeatures = 2000)


#load sample
OD142.DATA=Read10X_h5("OD14-2.h5")
rownames(OD142.DATA)<-gsub("_", "-", rownames(OD142.DATA))
#Generate Seurat object for further analysis
OD142 <- CreateSeuratObject(counts = OD142.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD142[["percent.mt"]] <- PercentageFeatureSet(object = OD142, pattern = "^mt-")
VlnPlot(object = OD142, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD142, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD142, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD142, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data; From QC not a good sample
OD142 <- subset(OD142,  subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 &nCount_RNA < 50000 & percent.mt < 20)
OD142 <- NormalizeData(OD142, verbose = FALSE)
OD142 <- FindVariableFeatures(OD142, selection.method = "vst", nfeatures = 2000)


#load sample
OD143.DATA=Read10X_h5("OD14-3.h5")
rownames(OD143.DATA)<-gsub("_", "-", rownames(OD143.DATA))
#Generate Seurat object for further analysis
OD143 <- CreateSeuratObject(counts = OD143.DATA)
# Take a look of mitrochondrial genes and nFeature genes, for gate out purposes;
OD143[["percent.mt"]] <- PercentageFeatureSet(object = OD143, pattern = "^mt-")
VlnPlot(object = OD143, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#QC step
plot1 <- FeatureScatter(OD143, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OD143, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OD143, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1))
CombinePlots(plots = list(plot2))
CombinePlots(plots = list(plot3))
#Setup gateout strategy for data; From QC not a good sample
OD143 <- subset(OD143,  subset = nFeature_RNA > 500 & nFeature_RNA < 7500 &nCount_RNA < 50000 & percent.mt < 20)
OD143 <- NormalizeData(OD143, verbose = FALSE)
OD143 <- FindVariableFeatures(OD143, selection.method = "vst", nfeatures = 2000)


save(YD01,file="YD01.rds")
save(YD02,file="YD02.rds")
save(YD03,file="YD03.rds")
save(YD41,file="YD41.rds")
save(YD42,file="YD42.rds")
save(YD43,file="YD43.rds")
save(YD44,file="YD44.rds")
save(OD01,file="OD01.rds")
save(OD02,file="OD02.rds")
save(OD03,file="OD03.rds")
save(OD41,file="OD41.rds")
save(OD42,file="OD42.rds")
save(OD43,file="OD43.rds")
save(OD44,file="OD44.rds")
save(YD141,file="YD141.rds")
save(YD142,file="YD142.rds")
save(YD143,file="YD143.rds")
save(OD141,file="OD141.rds")
save(OD142,file="OD142.rds")
save(OD143,file="OD143.rds")


load("YD01.rds")
load("YD02.rds")
load("YD03.rds")
load("YD41.rds")
load("YD42.rds")
load("YD43.rds")
load("YD44.rds")
load("YD141.rds")
load("YD142.rds")
load("YD143.rds")
load("OD01.rds")
load("OD02.rds")
load("OD03.rds")
load("OD41.rds")
load("OD42.rds")
load("OD43.rds")
load("OD44.rds")
load("OD141.rds")
load("OD142.rds")
load("OD143.rds")

YD01$age="1YD0"
YD02$group="1YD0"
YD03$group="1YD0"
YD41$group="1YD4"
YD42$group="1YD4"
YD43$group="1YD4"
YD44$group="1YD4"
YD141$group="1YDF4"
YD142$group="1YDF4"
YD143$group="1YDF4"
OD01$group="2OD0"
OD02$group="2OD0"
OD03$group="2OD0"
OD41$group="2OD4"
OD42$group="2OD4"
OD43$group="2OD4"
OD44$group="2OD4"
OD141$group="2ODF4"
OD142$group="2ODF4"
OD143$group="2ODF4"

YD01$group1="D0"
YD02$group1="D0"
YD03$group1="D0"
YD41$group1="D4"
YD42$group1="D4"
YD43$group1="D4"
YD44$group1="D4"
YD141$group1="DF4"
YD142$group1="DF4"
YD143$group1="DF4"
OD01$group1="D0"
OD02$group1="D0"
OD03$group1="D0"
OD41$group1="D4"
OD42$group1="D4"
OD43$group1="D4"
OD44$group1="D4"
OD141$group1="DF4"
OD142$group1="DF4"
OD143$group1="DF4"

YD01$group2="1Young"
YD02$group2="1Young"
YD03$group2="1Young"
YD41$group2="1Young"
YD42$group2="1Young"
YD43$group2="1Young"
YD44$group2="1Young"
YD141$group2="1Young"
YD142$group2="1Young"
YD143$group2="1Young"
OD01$group2="2Old"
OD02$group2="2Old"
OD03$group2="2Old"
OD41$group2="2Old"
OD42$group2="2Old"
OD43$group2="2Old"
OD44$group2="2Old"
OD141$group2="2Old"
OD142$group2="2Old"
OD143$group2="2Old"

YD01$group3="1YD01"
YD02$group3="1YD02"
YD03$group3="1YD03"
YD41$group3="1YD41"
YD42$group3="1YD42"
YD43$group3="1YD43"
YD44$group3="1YD44"
YD141$group3="1YDF41"
YD142$group3="1YDF42"
YD143$group3="1YDF43"
OD01$group3="2OD01"
OD02$group3="2OD02"
OD03$group3="2OD03"
OD41$group3="2OD41"
OD42$group3="2OD42"
OD43$group3="2OD43"
OD44$group3="2OD44"
OD141$group3="2ODF41"
OD142$group3="2ODF42"
OD143$group3="2ODF43"



#Find Anchors for data
AEC.anchors <- FindIntegrationAnchors(object.list = list(YD01, YD02, YD03, YD41, YD42, YD43, YD44,YD141, YD142, YD143, OD01, OD02, OD03, OD41, OD42, OD43, OD44, OD141, OD142, OD143), dims = 1:30)

#Combine data
AEC.combined <- IntegrateData(anchorset = AEC.anchors, dims = 1:30)


#Perform an integrated analysis
DefaultAssay(object = AEC.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
AEC.combined <- ScaleData(object = AEC.combined, features=rownames(AEC.combined), verbose = FALSE)
# Run the standard workflow for visualization and clustering
AEC.combined <- RunPCA(object = AEC.combined, npcs = 50, verbose = FALSE)


#Determine the ??dimensionality?? of the dataset
#Method 1 Jackstrawp
AEC.combined <- JackStraw(AEC.combined, num.replicate = 100)
AEC.combined <- ScoreJackStraw(AEC.combined, dims = 1:20)
JackStrawPlot(AEC.combined, dims = 1:20)
# Method 2
ElbowPlot(AEC.combined, ndims = 30)

#After the above step, set dims to 15
#t-SNE and Clustering
AEC.combined<- FindNeighbors(object = AEC.combined, reduction = "pca", dims = 1:15)
AEC.combined <- FindClusters(AEC.combined, resolution = 0.2)
AEC.combined <- RunUMAP(object = AEC.combined, reduction = "pca", dims = 1:15)

#Visualization of cells;
DimPlot(object = AEC.combined, reduction = "umap")
DimPlot(object = AEC.combined, reduction = "umap", label=T)
DimPlot(object = AEC.combined, reduction = "umap", split.by = "group", ncol=3)
DimPlot(object = AEC.combined, reduction = "umap", split.by = "group1")
DimPlot(object = AEC.combined, reduction = "umap", split.by = "group2")
DimPlot(object = AEC.combined, reduction = "umap", split.by = "group3",ncol=5)
DimPlot(object = AEC.combined, reduction = "umap", group.by = "group")
DimPlot(object = AEC.combined, reduction = "umap", group.by = "group1")
DimPlot(object = AEC.combined, reduction = "umap", group.by = "group2")
DimPlot(object = AEC.combined, reduction = "umap", group.by = "group3")

#Select cells to be sepecific cluster
DimPlot(object = AEC.combined, reduction = "umap", label = TRUE)
plot <- DimPlot(AEC.combined, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select.cells) <- "New1"

plot <- DimPlot(AEC.combined, reduction = "umap")
select2.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select2.cells) <- "New2"

plot <- DimPlot(AEC.combined, reduction = "umap")
select3.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select3.cells) <- "New3"

plot <- DimPlot(AEC.combined, reduction = "umap")
select4.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select4.cells) <- "New4"


plot <- DimPlot(AEC.combined, reduction = "umap")
select5.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select5.cells) <- "New5"


plot <- DimPlot(AEC.combined, reduction = "umap")
select6.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select6.cells) <- "New6"

plot <- DimPlot(AEC.combined, reduction = "umap")
select7.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select7.cells) <- "New7"

plot <- DimPlot(AEC.combined, reduction = "umap")
select8.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select8.cells) <- "New8"

plot <- DimPlot(AEC.combined, reduction = "umap")
select9.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select9.cells) <- "New9"

plot <- DimPlot(AEC.combined, reduction = "umap")
select10.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select10.cells) <- "New10"

plot <- DimPlot(AEC.combined, reduction = "umap")
select11.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select11.cells) <- "New11"

plot <- DimPlot(AEC.combined, reduction = "umap")
select12.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select12.cells) <- "New12"

plot <- DimPlot(AEC.combined, reduction = "umap")
select13.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select13.cells) <- "New13"

plot <- DimPlot(AEC.combined, reduction = "umap")
select14.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select14.cells) <- "New14"


plot <- DimPlot(AEC.combined, reduction = "umap")
select15.cells <- CellSelector(plot = plot)
Idents(AEC.combined, cells = select15.cells) <- "New15"



#Rename clusters
AEC.combined <- RenameIdents(AEC.combined, `0` = "0-AEC2-1",`New12` = "0-AEC2-1", `2` = "0-AEC2-2", `New11` = "0-AEC2-2", `1` = "0-AEC2-3",`New10` = "0-AEC2-3", `New9` = "1-AEC1",`5` = "1-AEC1",`New13` = "1-AEC1", `New8` = "2-Basal", `7` = "2-Basal", `4` = "3-Club-1",`New15` = "3-Club-1", `New14` = "3-Club-2", `New4` = "3-Club-2", `New3` = "4-Goblet", `New5` = "4-Goblet",`3` = "5-Ciliated",`9` = "5-Ciliated", `New2` = "6-CCSP+SPC+", `New7` = "7-PNEC", `8` = "8-Cycling", `New6` = "9-Fibroblasts", `New1` = "91-Immune", `6` = "Unknown")


##Figure 4K
AEC2.cells=subset(AEC.combined, idents=c("0-AEC2-1", "0-AEC2-2", "0-AEC2-3"))
DefaultAssay(AEC2.cells) <- "RNA"
VlnPlot(object = AEC2.cells, features =c("Slc39a8", "Sftpc", "Abca3"),pt.size =0 , cols = c("red","blue"), group.by = "group1", split.by ="group2", split.plot = TRUE, ncol=6)

##Figure 5M
D0.AEC2=subset(AEC2.cells, subset = group1 == "D0",idents=c("0-AEC2-1", "0-AEC2-2", "0-AEC2-3"))
VlnPlot(D0.AEC2, features = c("Chia1","H2-K1"), pt.size =0.0 ,cols= c("red", "blue"),group.by = "age", ncol=6)+ NoLegend()




