library(Signac)
library(Seurat)
library(ggplot2)
library(viridis)
library(ggpointdensity)
set.seed(1234)
seu <- readRDS('data/seu_no_qc.rds')
time_palette <- RColorBrewer::brewer.pal(9,"BuGn")[c(4,6,7,9)]

# Figure S1AB (QC metrics) ####
DefaultAssay(seu) <- "RNA"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
df <- seu@meta.data
plot_qc <- function(x,yintercept){
        ggplot(df,aes(stage,x))+geom_violin(aes(color=stage),alpha=0.5)+
                geom_point(aes(color=stage),position = position_jitter(seed = 1, width = 0.2))+
                geom_boxplot(aes(color=stage),width=0.2, color="grey", alpha=0.2)+theme_bw()+
                scale_color_manual(values = time_palette)+
                theme(legend.position = "none")+
                geom_hline(yintercept = yintercept,linetype="dashed", color = "gray", size=1)+
                xlab("Developmental stage")
}
p1 <- plot_qc(x=log10(df$nCount_RNA),yintercept = 5)+ylab("Log10(# of reads)")+
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p2 <- plot_qc(x=df$nFeature_RNA,yintercept = 7500)+ylab("# of genes")+
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p3 <- plot_qc(x=df$percent.mt,yintercept = 20)+ylab("% of MT")
p4 <- plot_qc(x=log10(df$nCount_ATAC),yintercept = 5)+ylab("Log10(# of fragments)")+
        geom_hline(yintercept = log10(200),linetype="dashed", color = "gray", size=1)+
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p5 <- plot_qc(x=df$nucleosome_signal,yintercept = 2)+ylab("Nucleosome signal")+
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p6 <- plot_qc(x=df$TSS.enrichment,yintercept = 1)+ylab("TSS enrichment")
(p1/p2/p3)|(p4/p5/p6)

# Figure S1B ####
ggplot(df,aes(log10(nCount_ATAC),TSS.enrichment))+geom_point()+geom_pointdensity()+
        scale_color_viridis()+facet_wrap(~stage)+theme_bw()+
        geom_hline(yintercept = 1,linetype="dashed", color = "darkgray", size=0.5)+
        geom_vline(xintercept = c(log10(200),5),linetype="dashed", color = "darkgray", size=0.5)+
        xlab("Log10(# of fragments)")+ylab("TSS enrichment")

# Figure S1D ####
seu <- readRDS('data/seu.rds')
library(corrplot)
library(dplyr)
rna <- seu@assays$SCT@scale.data
asplit <- split(colnames(seu), seu$library)
aggr <- do.call(cbind, lapply(names(asplit), function(x) {
        rowMeans(rna[,asplit[[x]]])
}))
colnames(aggr) <- names(asplit)
cor <- cor(aggr)
corrplot(cor, tl.col = "black", tl.srt = 45,col.lim = c(min(cor),1),is.corr = F,
         #addCoef.col = T,
         col=colorRampPalette(c("lightgray","#238B45"))(20)) %>% corrRect(c(1,3,6,8,9))
        

atac <- t(seu@reductions$lsi@cell.embeddings)
asplit <- split(colnames(seu), seu$library)
aggr <- do.call(cbind, lapply(names(asplit), function(x) {
        rowMeans(atac[,asplit[[x]]])
}))
colnames(aggr) <- names(asplit)
cor <- cor(aggr)
corrplot(cor, tl.col = "black", tl.srt = 45,col.lim = c(min(cor),1),is.corr = F,
         col=colorRampPalette(c("lightgray","#238B45"))(20)) %>% corrRect(c(1,3,6,8,9))

# Figure S1E ####
seu <- readRDS('data/seu.rds')
p1 <- DimPlot(seu,group.by = 'library',split.by='stage',reduction='rna.umap',label=F)+labs(title="RNA")
p2 <- DimPlot(seu,group.by = 'library',split.by='stage',reduction='atac.umap',label=F)+labs(title='ATAC')
p1/p2

# Figure S1C ####
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
set.seed(1234)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"
counts <- Read10X_h5("data/raw_data/aggr/filtered_feature_bc_matrix.h5")
fragpath <- "data/raw_data/aggr/atac_fragments.tsv.gz"
tmp <- counts$Peaks
tmp <- tmp[sample(rownames(tmp),100),sample(colnames(tmp),100)]
chrom_assay <- CreateChromatinAssay(counts = tmp,sep = c(":", "-"),
                                      fragments = fragpath,annotation = annotation,
                                      genome="mm10")
tmp <- CreateSeuratObject(counts = chrom_assay,assay = "peaks")
tmp <- TSSEnrichment(tmp,fast=FALSE)
TSSPlot(tmp)
FragmentHistogram(tmp,group.by = 'stage')

