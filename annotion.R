#注释
plot1 <- clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
ggsave("/home/zmzhang/schlhproj/QC/clustree.png", plot = plot1, width = 12, height = 8)

scRNA@meta.data$seurat_clusters <- scRNA@meta.data$RNA_snn_res.0.1
Idents(scRNA) <- scRNA@meta.data$seurat_clusters

scRNA <- low_con_scRNA


library(tidyverse)
library(ggplot2)
library(viridis)
library(Seurat)

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#readRDS("/home/zmzhang/schlhproj/script/Umap-scHLH.RData")
DefaultAssay(scRNA) <- "integrated"
genemapData <- read.table(file="/home/zmzhang/schlhproj/script/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/home/zmzhang/schlhproj/script/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(genelenData[,1],genemapData[,1])
lenInfo <- cbind(genemapData[matchIndexes[which(!is.na(matchIndexes))],3],genelenData[which(!is.na(matchIndexes)),2])

rawData <- GetAssayData(object=scRNA,layer="counts")
normData <- GetAssayData(object=scRNA,layer="data") 
cellmatchIndexes <- match(colnames(normData),colnames(rawData))

rawData <- rawData[,cellmatchIndexes]

matchIndexes <- match(rownames(rawData),lenInfo[,1])
countData <- as.matrix(rawData[which(!is.na(matchIndexes)),])
lenVec <- lenInfo[matchIndexes[which(!is.na(matchIndexes))],2]
tpmData <- tpm3(countData,as.numeric(lenVec))

query <- t(tpmData)
model <- readRDS(file="/home/zmzhang/schlhproj/script/major_human_cell_types.rds")
prd <- LoadModel(model)
label <- prd(query)
orgin <- Idents(scRNA)
orgin <- as.numeric(orgin) - 1
orgin <- paste("C",orgin,sep="")

pdf("/home/zmzhang/schlhproj/png/0.3predcelltypes_heatmap_ref.pdf",width=10,height=10)
Confusion_heatmap(orgin,label)
dev.off()

de <- DimPlot(scRNA, reduction = "umap", label=T,group.by = 'celltype')
ggsave("/home/zmzhang/schlhproj/QC/celltype.png", plot = de, width = 12, height = 8)

##手动注释
genes <- c("Cd160","Cd8","Entpd1","Pdcd1","Gzmb","Cx3cr1")
genes <- c("Slc17a7","NRGN","Slc30a3","Slc17a6","Baiap3","Tcf4","Gad1","Gad2","Ctss","C1qb","C1qa","Slc1a3","Gja1","Aqp4","Acsbj1","Aspa","Ermn","Mog","Acsbg1","Cacng4","Flt1","Fldn5")
#画dotplot图，注释利器之一
plot2 <- DotPlot(scRNA, features = genes,assay='RNA') + coord_flip()
#ggsave("/home/zmzhang/schlhproj/QC/DotPlot.png", plot = plot2, width = 12, height = 8)
plot2

FeaturePlot(scRNA, features = c("Cd8b1","Cd8a","Entpd1","Pdcd1","Gzmb","Cx3cr1")
            )


#细胞比例图


Subset_color_panel <- c(
  # Bn
  "c01_Bn_TCL1A" = "#e5d25b",
  "c02_Bn_NR4A2" = "#599014",
  "c03_Bn_IFN-response" = "#e78071",
  # Bm
  "c04_classical-Bm_TXNIP" = "#a82d06",
  "c05_classical-Bm_GPR183" = "#4592bf",
  "c06_Bm_stress-response" = "#d38219",
  "c07_Bm_IFN-response" = "#74a764",
  "c08_ABC_FCRL4" = "#8ca2b4",
  "c09_ABC_FGR" = "#cbb190",
  "c10_Bm_TCL1A" = "#e7ca8d",
  "c11_pre-GC" = "#9d9ec3",
  # Bgc
  "c12_Bgc_LZ-like" = "#593202",
  # ASC
  "c16_PC_IGHG" = "#ebafa4",
  "c17_PC_IGHA" = "#5e8a89",
  "c18_early-PC_MS4A1low" = "#ecd577",
  "c19_early-PC_LTB" = "#7c606c",
  "c20_early-PC_RGS13" = "#5c6489",
  # Bcycling
  "c13_Bgc_DZ-like" = "#ECE4B7",
  "c15_cycling_ASC" = "#D36135",
  "c14_Bm_activated-cycling" = "#467599"
)



annotation <- read_excel("/home/zmzhang/schlhproj/QC/4annotation_bio4.xlsx")
Idents(scRNA) <- scRNA@meta.data$seurat_clusters
new_ids <- annotation$celltype
Idents(scRNA) <- scRNA@meta.data$seurat_clusters
names(new_ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new_ids)
scRNA@meta.data$celltype <- Idents(scRNA)



table(scRNA$orig.ident)#查看各组细胞数
prop.table(table(Idents(scRNA)))
table(Idents(scRNA), scRNA$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(scRNA), scRNA$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

allcolour=c("T cell" = "#e5d25b","NK" = "#a82d06","Monocyte" = "#4592bf","B cell" = "#ECE4B7",
            "DC" = "#9d9ec3",
            "Neutrophil" = "#D36135","Platelet"="#599014")
library(ggplot2)
wc <- ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.3,position = "fill")+
  cowplot::theme_cowplot()+ 
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 5),
    text = element_text(size = 12),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    axis.line = element_line(size = 0.6),
    axis.ticks = element_line(size = 0.6),
    legend.position = "bottom"
  ) 
wc
  
ggsave("/home/zmzhang/schlhproj/QC/proportion.png", plot = wc, width = 4.2, height = 4.5)
save(scRNArmdouble,file = "/home/zmzhang/schlhproj/intergrate/doublefinderfilter2.RData")
save(scRNAnonfilter,file = "/home/zmzhang/schlhproj/intergrate/nonefilter2.RData")
##流向堆叠
library(reshape2)
library(ggplot2)
library(ggalluvial)
library(ggh4x)
Cellratio <- prop.table(table(Idents(subCells), subCells$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colors=c("#FBB463","#80B1D3","#F47F72","#BDBADB","#FBF8B4","#8DD1C6")
p <- ggplot(Cellratio, aes(x = Var2, y = Freq, fill = Var1,
                           stratum = Var1, alluvium = Var1)) +
  geom_col(position = 'stack', width = 0.6) +
  geom_stratum(width = 0.6, color = 'white') +
  geom_alluvium(alpha = 0.4, width = 0.6, color = 'white', linewidth = 1, curve_type = "linear") +
  scale_fill_manual(values = colors) +
  xlab('') + 
  ylab('') +
  scale_y_continuous(expand = c(0, 0))+
  theme_bw(base_size = 12) + 
  theme(
    axis.text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
p
ggsave("stacked.pdf", plot = p, height = 4, width = 5)
