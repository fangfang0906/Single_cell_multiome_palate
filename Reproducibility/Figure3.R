library(Seurat)
library(Signac)
library(ggplot2)
set.seed(1234)
time_palette <- RColorBrewer::brewer.pal(9,"BuGn")[c(4,6,7,9)]

# Figure 3A ####
seu_overall <- readRDS('data/seu.rds')
seu_overall$tmp <- "no"
seu_overall$tmp[seu_overall$celltype_major=="CNC-derived mesenchymal"] <- "yes"
DimPlot(seu_overall,reduction='rna.umap',group.by = "tmp",
        cols=c("darkgrey","steelblue"))+NoLegend()

seu_overall$tmp <- "no"
seu_overall$tmp[seu_overall$celltype_major=="Epithelial"] <- "yes"
DimPlot(seu_overall,reduction='rna.umap',group.by = "tmp",
        cols=c("darkgrey","steelblue"))+NoLegend()

seu <- readRDS("data/subset_w_chondrogenic.rds")
seu <- FindClusters(seu,resolution=0.6)
DimPlot(seu, reduction="rna.umap",label = TRUE, repel = TRUE) +NoLegend()
DimPlot(seu, reduction="rna.umap",label = F,group.by = 'stage',cols=time_palette)
DimPlot(seu, reduction="rna.umap",label = F,group.by = 'celltype')+NoLegend()

# Figure 3B ####
FeaturePlot(seu,c("Shox2","Msx1","Meox2","Tbx22"),ncol=2)
FeaturePlot(seu,c("Osr2","Fgf10","Fgf7","Dlx5"),ncol=2)

# Figure 3C-F ####
# load bulk datasets ####
data <- read.csv("star_counts.csv",row.names=1)
group <- factor(unlist(lapply(colnames(data),function(x){substr(x,0,1)})),
                levels=c("P","A"))
levels(group) <- c("Posterior","Anterior")
colData <- data.frame(row.names = colnames(data), group)
countData <- data[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ group)
dds <- DESeq(dds)
res_bulk <- na.omit(results(dds))

# PCA ####
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggrepel)
ggplot(pcaData, aes(PC1, PC2, color=group, shape=group,label=name)) +
        geom_point(size=2) + geom_label_repel(size=3)+
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        ggtitle("PCA")+theme_bw()

# add scRNAseq results ####
seu <- readRDS('seu_RNA_only.rds')
seu$tmp <- NA
seu$tmp[seu$seurat_clusters==2] <- "posterior"
seu$tmp[seu$seurat_clusters %in% c(13,6)] <- "anterior"
seu <- seu[,!is.na(seu$tmp)]
seu <- NormalizeData(seu)
Idents(seu) <- seu$tmp
res_sc <- FindMarkers(seu,min.pct = 0.1,min.diff.pct=0.1,
                      logfc.threshold = 0.1,
                      ident.1="anterior",ident.2 = "posterior")

# fisher test ####
res_bulk <- data.frame(res_bulk)
res_bulk$bulk_sig <- NA
res_bulk$bulk_sig[res_bulk$padj<0.05 & res_bulk$log2FoldChange>0] <- "bulk, up"
res_bulk$bulk_sig[res_bulk$padj<0.05 & res_bulk$log2FoldChange<0] <- "bulk, down"
res_bulk$sc_sig <- NA
res_bulk[rownames(res_bulk) %in% 
                 rownames(res_sc)[res_sc$p_val_adj< 0.05 & res_sc$avg_log2FC>0],'sc_sig'] <- "sc, up"
res_bulk[rownames(res_bulk) %in% 
                 rownames(res_sc)[res_sc$p_val_adj< 0.05 & res_sc$avg_log2FC<0],'sc_sig'] <- "sc, down"

table(res_bulk$bulk_sig,res_bulk$sc_sig)
fisher.test(table(res_bulk$bulk_sig,res_bulk$sc_sig))
# Volcano plot
up <- rownames(res_sc[res_sc$avg_log2FC>0,])
down <- rownames(res_sc[res_sc$avg_log2FC<0,])
res_bulk$label <- rownames(res_bulk)
res_bulk$label[!(rownames(res_bulk) %in% c("Shox2","Meox2"))] <- NA
ggplot()+geom_point(data=subset(res_bulk,select=-sc_sig), 
                    aes(x = log2FoldChange, y = -log10(padj)),color="grey")+
        geom_point(data=res_bulk, aes(x = log2FoldChange, y = -log10(padj),color=sc_sig))+
        geom_point() + theme_bw()+ ggtitle("Anterior vs Posterior")+
        facet_wrap(~sc_sig)

# PCA projection ####
degs <- c(up,down)
expr <- seu@assays$RNA@counts[degs,seu$tmp %in% c("anterior","posterior")]
meta <- seu@meta.data
asplit <- split(rownames(meta), meta$tmp)
pseudobulks <- do.call(cbind, lapply(asplit, function(x) rowMeans(expr[,x])))
pseudobulks <- expr
# match genes
ok <- intersect(rownames(pseudobulks),rownames(data))
pseudobulks_ok <- pseudobulks[ok,]
bulk_ok <- data[ok,]
# voom 
library(limma)
pseudobulks_ok <- voom(pseudobulks_ok)$E
bulk_ok <- voom(bulk_ok)$E
merge <- cbind(pseudobulks_ok,bulk_ok)
library('preprocessCore')
norm_merge <- normalize.quantiles(merge)
colnames(norm_merge) <- colnames(merge)
rownames(norm_merge) <- rownames(merge)
# pca
tmp <- norm_merge[, colnames(pseudobulks)]
pca <- prcomp(t(tmp))
# project bulk into pca sapce
preds <- predict(pca, t(norm_merge[, colnames(data)]))
df <- as.data.frame(rbind(pca$x,preds))
df$group <- "single cell RNA-seq"
df[colnames(data),'group'] <- "bulk RNA-seq"
df$location <- meta[match(rownames(df),rownames(meta)),'tmp']
df$location[rownames(df) %in% c("A1","A2","A3")] <- "anterior"
df$location[rownames(df) %in% c("P1","P2","P3")] <- "posterior"
ggplot()+geom_point(data=subset(df,select=-group),aes(PC1,PC2),color="grey")+
        geom_point(data=df,aes(PC1,PC2,color=location,shape=group))+
        theme_bw()+facet_wrap(~group)

# Boxplot ####
plotGene <- function(gene){
        df <- data.frame(expr=norm_merge[gene,])
        df$group <- "single cell RNA-seq"
        df[colnames(data),'group'] <- "bulk RNA-seq"
        df$location <- meta[match(rownames(df),rownames(meta)),'tmp']
        df$location[rownames(df) %in% c("A1","A2","A3")] <- "anterior"
        df$location[rownames(df) %in% c("P1","P2","P3")] <- "posterior"
        ggplot(df,aes(location,expr,color=location))+
                #geom_jitter(color='grey', size=2, alpha=0.9)+
                geom_boxplot()+facet_wrap(~group)+
                theme_bw()+ggtitle(gene)  
}
plotGene("Shox2")|plotGene("Meox2")
