library(Seurat)
library(ggplot2)
set.seed(1234)
# load data ####
subset <- readRDS('data/subset_w_chondrogenic.rds')
time_palette <- RColorBrewer::brewer.pal(9,"BuGn")[c(4,6,7,9)]

# S4A ####
subset$tmp <- "no"
subset$tmp[subset$stage=="E12.5"] <- "yes"
p1 <- DimPlot(subset,reduction='rna.umap',group.by = 'tmp',
        cols = c("lightgrey",time_palette[1]))+NoLegend()
subset$tmp <- "no"
subset$tmp[subset$stage=="E13.5"] <- "yes"
p2 <- DimPlot(subset,reduction='rna.umap',group.by = 'tmp',
              cols = c("lightgrey",time_palette[2]))+NoLegend()
subset$tmp <- "no"
subset$tmp[subset$stage=="E14.0"] <- "yes"
p3 <- DimPlot(subset,reduction='rna.umap',group.by = 'tmp',
              cols = c("lightgrey",time_palette[3]))+NoLegend()
subset$tmp <- "no"
subset$tmp[subset$stage=="E14.5"] <- "yes"
p4 <- DimPlot(subset,reduction='rna.umap',group.by = 'tmp',
              cols = c("lightgrey",time_palette[4]))+NoLegend()
p1|p2|p3|p4

# S4B ####
DefaultAssay(subset) <- "SCT"
FeaturePlot(subset,reduction="rna.umap",ncol=4,
            c("Shox2","Alx1",'Col2a1','Sp7','Tfap2b',"Tbx22","Wnt16","Tbx15"))