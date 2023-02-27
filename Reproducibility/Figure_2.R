library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(pheatmap)
library(zoo)
library(ggplot2)

# Load data ####
seu <- readRDS('data/seu.rds')

# downsample ####
sample <- unlist(lapply(unique(seu$celltype_major),function(x){
        cell <- colnames(seu)[seu$celltype_major==x]
        if(length(cell)>2000){sample(cell,2000)} else {cell}
}))
seu <- seu[,sample]

# Putative regulators ####
DefaultAssay(seu) <- 'peaks'
markers_rna <- presto:::wilcoxauc.Seurat(seu, group_by = 'celltype_major', assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(seu, group_by = 'celltype_major', assay = 'data', seurat_assay = 'chromvar')
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(seu, id = markers_motifs$motif.feature)
# additional processing 
library(tidyr)
markers_motifs$gene <- gsub("\\(var.2)|\\(var.3)","",markers_motifs$gene)
markers_motifs <- markers_motifs %>% 
        mutate(gene = strsplit(gene, "::")) %>% 
        unnest(gene)
markers_motifs <- data.frame(markers_motifs)
markers_motifs$gene <- unlist(lapply(markers_motifs$gene,function(x){
        paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}))
topTFs <- function(celltype, padj.cutoff = 0.05) {
        ctmarkers_rna <- dplyr::filter(
                markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
                arrange(-RNA.auc)
        ctmarkers_motif <- dplyr::filter(
                markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
                arrange(-motif.auc)
        top_tfs <- inner_join(
                x = ctmarkers_rna[, c(2, 11, 6, 7)], 
                y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
        )
        top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
        top_tfs <- arrange(top_tfs, -avg_auc)
        return(top_tfs)
}
head(topTFs("Epithelial"), 3)

# Find DEGs bw clusters ####
DefaultAssay(seu) <- "SCT"
Idents(seu) <- seu$celltype_major
de_genes <- FindAllMarkers(seu,only.pos = T,logfc.threshold = 0.1,max.cells.per.ident = 1000)

# Linking peaks to genes ####
DefaultAssay(seu) <- "peaks"
seu <- RegionStats(seu, genome = BSgenome.Mmusculus.UCSC.mm10)
deg <- unique(de_genes[de_genes$p_val_adj<0.05 & de_genes$avg_log2FC>0.1,'gene'])
seu <- LinkPeaks(seu,peak.assay = "peaks",expression.assay = "SCT",genes.use = deg)
link <- Links(seu)
link <- data.frame(link)
link$adj_pval <- p.adjust(link$pvalue,method = 'BH')
link$gene_cluster <- de_genes[match(link$gene,de_genes$gene),'cluster']
lvl <- levels(seu)
link <- arrange(link,factor(gene_cluster,levels=lvl))
link <- link[link$adj_pval<0.05 & link$score>0,]
t <- link %>% group_by(gene) %>% summarise(n=n())
summary(t$n)

# Heatmap ####
link_1 <- link[link$score>0.1,]
link_2 <- link[link$score>0 & link$gene_cluster %in% c("Neuronal","Myocytes"),]
link <- rbind(link_1,link_2)
link <- arrange(link,factor(gene_cluster,levels=lvl))
# gene expression
gene.use <- unique(link$gene)
DefaultAssay(seu) <- "RNA"
seu <- ScaleData(seu,vars.to.regress = c("S.Score", "G2M.Score"),features=gene.use)
rna <- seu@assays$RNA@scale.data[gene.use,]
asplit_cells <- split(colnames(seu), seu$celltype_major)
n <- 10
means <- do.call(cbind, lapply(lvl, function(x){
        df <- rna[,asplit_cells[[x]]]
        t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col <- data.frame(celltype_major)
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
anno_col <- arrange(anno_col,factor(celltype_major,levels=lvl))
anno_col$celltype_major <- factor(anno_col$celltype_major,levels=lvl)
pheatmap(means[link$gene,rownames(anno_col)],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,
         #breaks=seq(-2,2,length=101),
         annotation_col = anno_col,show_colnames = F, show_rownames = F,
         main='Gene expression')
# chromatin accessibility 
peak.use <- unique(link$peak)
DefaultAssay(seu) <- "peaks"
atac <- RunTFIDF(seu@assays$peaks@counts[peak.use,])
asplit_cells <- split(colnames(seu), seu$celltype_major)
means_peak <- do.call(cbind, lapply(lvl, function(x){
        df <- atac[,asplit_cells[[x]]]
        t(apply(df,1,function(x){rollapply(x,n,mean,by=n)}))
}))
celltype_major <- unlist(lapply(lvl, function(x) rep(x, length(asplit_cells[[x]])%/%n)))
anno_col2 <- data.frame(celltype_major)
rownames(anno_col2) <- colnames(means_peak) <- paste0(seq(1:ncol(means_peak)),colnames(means_peak))
anno_col2 <- arrange(anno_col2,factor(celltype_major,levels=lvl))
pheatmap(means_peak[link$peak,rownames(anno_col2)],cluster_rows = F, cluster_cols = F, scale = "row",show_rownames = F,
         breaks = seq(-2, 2, length = 101), color = hcl.colors(100, "BluYl"),
         annotation_col = anno_col,show_colnames = F,
         main='Chromatin accessibility')

# Correlation Plot ####
library(ggpubr)
tmp <- link[link$gene_cluster=="CNC-derived mesenchymal",]
gene <- "Twist1"
tmp <- na.omit(tmp[tmp$gene==gene,])
tmp <- arrange(tmp,-score)
peak1 <- tmp$peak[1]
region <- unlist(lapply(tmp[,'peak'],function(x){as.numeric(strsplit(x,"-",fixed=T)[[1]][c(2,3)])}))
peak2 <- paste(strsplit(tmp[1,'peak'],"-",fixed=T)[[1]][1],
             min(region,na.rm = T),
             max(region,na.rm=T),sep="-")
df <- data.frame(gene_expr=means[gene,],celltype=anno_col$celltype_major)
p1 <- ggplot(df,aes(celltype,gene_expr,col=celltype))+geom_boxplot()+
        theme_bw()+ggtitle(paste0(gene, " expression"))+
        xlab("")+ylab("Relative expression")+
        theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())
df <- data.frame(atac=means_peak[peak1,],celltype=anno_col2$celltype_major)
p2 <- ggplot(df,aes(celltype,atac,col=celltype))+geom_boxplot()+
        theme_bw()+ggtitle(paste0(peak1, " accessibility"))+
        xlab("")+ylab("Chromatin accessibility")+
        theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())
df <- data.frame(gene_expr=means[gene,],atac=means_peak[peak1,])
p3 <- ggplot(df,aes(gene_expr,atac))+geom_point()+geom_smooth(color="grey")+
        theme_bw()+ggtitle(paste0(gene, " & ", peak1))+
        xlab("gene expression")+ylab("chromatin accessibility")+
        stat_cor(method = "spearman", label.x = -0.5, label.y = 3.5)
p1/p2/p3








