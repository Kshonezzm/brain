
# 设置工作目录并读取样本文件夹名称
setwd('/home/zmzhang/schlhproj/brain2')
sample_dirs <- list.files()
samples <- sample_dirs[1:2]  


library(Seurat)
library(stringr)
library(ggplot2)
library(decontX)
library(DoubletFinder)
library(tidyverse)
library(dplyr)
library(patchwork)


seurat_list <- list()


for (i in seq_along(samples)) {
  
  sample_data <- Read10X(data.dir = samples[i])
  
 
  seurat_obj <- CreateSeuratObject(counts = sample_data)
  
  
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = samples[i])
  
 
  seurat_obj$orig.ident <- samples[i]
  
 
  seurat_list[[i]] <- seurat_obj
}


merged_seurat <- Reduce(function(x, y) merge(x, y, merge.data = TRUE), seurat_list)


print(merged_seurat)




#EBVadLA2<-CreateSeuratObject(counts = scRNA22EBVrmLA2)
#scRNA2 = merge(EBVadLA, y = c(EBVadLA2), 
#         add.cell.ids = c("EBV", 'EBVadLA'),
#            project = "HLH",merge.data = TRUE)


scRNA2 <- merged_seurat


#制作分组信息
#orig=as.data.frame((str_split(rownames(scRNA2@meta.data),pattern = '_')))
#orig=as.data.frame(t(orig))
##把分组信息加到meta.data中，变成orig.ident
#scRNA2@meta.data$orig.ident=orig$V1
#table(scRNA2@meta.data$orig.ident) 
#head(scRNA2@meta.data)

#计算线粒体
scRNA2[["percent.mt"]] <- PercentageFeatureSet(scRNA2, pattern = "^mt-")
col.num <- length(levels(as.factor(scRNA2@meta.data$orig.ident)))
minGene=100
maxGene=8500
pctMT=5


scRNA2 <- subset(scRNA2, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(as.factor(scRNA2@meta.data$orig.ident)))


scRNA2 <- NormalizeData(scRNA2)
scRNA2 <- FindVariableFeatures(scRNA2, nfeatures = 3000)
scRNA2 <- ScaleData(scRNA2)
scRNA2 <- RunPCA(scRNA2, npcs = 20) 

scRNA2 <- IntegrateLayers(
  object = scRNA2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
scRNA2[["RNA"]] <- JoinLayers(scRNA2[["RNA"]])

scRNA2 <- FindNeighbors(scRNA2, reduction = "pca", dims = 1:20)
scRNA2 <- FindClusters(scRNA2, resolution = seq(from = 0.1, to = 0.8, by = 0.1))

plot2 <- ElbowPlot(scRNA2, ndims=20, reduction="pca") 

pc.num=1:20

scRNA2 <- RunUMAP(scRNA2, dims = pc.num,
                  reductions = "integrated.cca"
                 #spread = 0.5,
                 #min.dist = 1,
)



save(scRNA2,file = "/home/zmzhang/schlhproj/intergrate/FILTER2.RData")
