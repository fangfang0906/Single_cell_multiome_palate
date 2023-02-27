library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
set.seed(1234)
seu <- readRDS('data/seu.rds')
time_palette <- RColorBrewer::brewer.pal(9,"BuGn")[c(4,6,7,9)]

# Figure S2A ####
DefaultAssay(seu) <- "SCT"
FeaturePlot(seu,reduction='rna.umap',ncol=4,
            features=c("Prrx1","Krt14","Epcam","Cdh5","Cldn5","Lyz2",
                       "Plp1","Sox10","Myod1","Tubb3","Stmn2","Myf5"))
DefaultAssay(seu) <- "Gene_activity"
Idents(seu) <- seu$celltype_major
VlnPlot(seu,stack=T,features=c("Prrx1","Krt14","Epcam"))
