library(Signac)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(viridis)
library(zoo)
library(pheatmap)
set.seed(1234)

# Load data ####
time_palette <- RColorBrewer::brewer.pal(9,"BuGn")[c(4,6,7,9)]
seu <- readRDS("data/subset_no_chondrocyte_0617.rds")
DimPlot(seu,group.by='celltype')+NoLegend()

# Load diffusion pseudotime & terminal prob ####
meta <- read.delim("python/CellRank+WOT_CNC_0617_no_chondrocyte/Cellrank+WOT_meta.txt",
                   head=T,sep=",",row.names = 1)
# Load absorption rate ####
prob <- read.delim("python/CellRank+WOT_CNC_0617_no_chondrocyte/cellrank_WOT_terminal_prob.txt",
                   head=T,sep=",",row.names = 1)
seu@meta.data <- cbind(meta[colnames(seu),],prob[colnames(seu),])

# load ancestors and descendants ####
traj <- read.delim("python/CellRank+WOT_CNC_0617_no_chondrocyte/ance_des.txt",
                   head=T,sep=",",row.names = 1)
colnames(traj) <- paste0("transition_",colnames(traj))
seu@meta.data <- cbind(seu@meta.data,traj[colnames(seu),])

seu$tmp <- seu$transition_anterior > summary(seu$transition_anterior)['Mean']
DimPlot(seu,group.by = 'tmp',reduction="rna.umap",cols=c('lightgrey','steelblue')) &
        labs(title="transition prob to anterior > Mean")

tmp <- seu[,seu$tmp]

seu$tmp1 <- seu$dpt_pseudotime
seu$tmp1[!(colnames(seu) %in% colnames(tmp))] <- NA
p1 <- FeaturePlot(seu, 'tmp1',reduction="rna.umap") & scale_color_viridis_c() & labs(title="diffusion pseudotime")
seu$tmp1 <- seu$anterior.palatal.mesenchymal
seu$tmp1[!(colnames(seu) %in% colnames(tmp))] <- NA
p2 <- FeaturePlot(seu, 'tmp1',reduction="rna.umap") & scale_color_viridis_c() & labs(title="absorption rate to anterior")
seu$tmp1 <- seu$terminal_states_probs
seu$tmp1[!(colnames(seu) %in% colnames(tmp))] <- NA
p3 <- FeaturePlot(seu, 'tmp1',reduction="rna.umap") & scale_color_viridis_c() & labs(title="terminal state probabilities")
p1|p2|p3

# load Cellrank+WOT driver genes ####
res <- read.delim("python/CellRank+WOT_CNC_0617_no_chondrocyte/Cellrank+WOT_driver_genes.txt",head=T,sep=",",row.names = 1)
res <- data.frame(res)
res <- na.omit(res[,c(1,2,6,7,11,12,16,17,21,22)])
corr.thres <- 0.05
# anterior
anterior <- na.omit(res[order(res$anterior.palatal.mesenchymal_corr,decreasing = T),])
anterior$adjusted_pval <- p.adjust(anterior$anterior.palatal.mesenchymal_pval,method="BH")
anterior <- anterior[anterior$adjusted_pval < 0.05 & anterior$anterior.palatal.mesenchymal_corr > corr.thres,]

# Heatmap ####
gene.use <- rownames(anterior)
pt <- tmp$dpt_pseudotime
DefaultAssay(tmp) <- "RNA"
tmp <- ScaleData(tmp,vars.to.regress = c("S.Score", "G2M.Score"),features=gene.use)
rna <- tmp@assays$RNA@scale.data[gene.use,names(pt)[order(pt)]]
n <- 100
means <-t(apply(rna,1,function(x){rollapply(x,n,mean,by=n)}))
anno_col <- data.frame(diffusion_pseudotime = rollapply(pt[order(pt)],n,mean,by=n),
                       fate_probability = rollapply(tmp$anterior.palatal.mesenchymal[order(pt)],n,mean,by=n),
                       terminal_likelihood = rollapply(tmp$terminal_states_probs[order(pt)],n,mean,by=n),
                       developmental_stage = rollapply(tmp$stage_new[order(pt)],n,mean,by=n))
rownames(anno_col) <- colnames(means) <- paste0(seq(1:ncol(means)),colnames(means))
peaks <- apply(means, 1, function(x) which(x == max(x)))
rowOrd <- order(peaks)
# cut based on quantiles ####
q <- quantile(anno_col$diffusion_pseudotime,probs = seq(0,1,0.33))
cut1 <- max(which(anno_col$diffusion_pseudotime<=q[2]))
cut1 <- length(peaks[peaks<=cut1])
cut2 <- max(which(anno_col$diffusion_pseudotime<=q[3]))
cut2 <- length(peaks[peaks<=cut2])
pheatmap(means[rowOrd,],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks=seq(-2,2,length=101),
         #gaps_row = c(cut1,cut2),
         annotation_col = anno_col,show_colnames = F, show_rownames = F)

# clusters 
genes <- gene.use[rowOrd]
t1 <-genes[1:cut1]
t2 <- genes[cut1:cut2]
t3 <- genes[cut2:length(genes)]

# single gene expr plot along pseudotime ####
library(reshape2)
plot_multipleGene <- function(gene.use){
        df <- data.frame(expr=t(means[gene.use,]),diffusion_pseudotime=anno_col$diffusion_pseudotime)
        melted <- melt(df,id.vars = 'diffusion_pseudotime')
        ggplot(melted, aes(y = value, x = diffusion_pseudotime,group=variable)) +
                geom_point()+theme_bw()+
                geom_smooth(method = "loess", span = 1,color='darkgrey') +
                ylab("Relative expression") + xlab("Diffusion pseudotime")
}
plot_singleGene <- function(gene){
        df <- data.frame(expr=means[gene,],diffusion_pseudotime=anno_col$diffusion_pseudotime)
        ggplot(df, aes(y = expr, x = diffusion_pseudotime)) +
                geom_point()+theme_bw()+
                geom_smooth(method = "loess", span = 1,color='darkgrey') +
                ylab("Relative expression") + xlab("Diffusion pseudotime")
}
p1 <- plot_multipleGene(t1[1:10])+ggtitle("Start")
p2 <- plot_multipleGene(t2[1:10])+ggtitle("Middle")
p3 <- plot_multipleGene(tail(t3,10))+ggtitle("End")
p1/p2/p3


# Enrichment analysis ####
library("WebGestaltR")
runEnrich <- function(gene,thr,name){
        WebGestaltR(enrichMethod="ORA", organism="mmusculus",
                    enrichDatabase=c("geneontology_Biological_Process_noRedundant"),
                    enrichDatabaseType="genesymbol",
                    interestGene=gene,
                    minNum=5,maxNum = 500,
                    interestGeneType="genesymbol",fdrThr = thr,
                    referenceSet ="genome_protein-coding",
                    referenceGeneType="genesymbol",isOutput = T,
                    #outputDirectory = "GSEA/",
                    projectName=name)  
}
lst <- list(t1,t2,t3)
names(lst) <- c("start","middle","end")
df <- do.call(rbind,lapply(names(lst),function(x){
        write.table(lst[[x]],paste0("GSEA/anterior_",x,".txt"),row.names = F,col.names = F,quote=F,sep="\t")
        a <- runEnrich(gene=lst[[x]],thr=0.05,name=paste0("anterior_",x))
        #if (is.null(a)) {a <- runEnrich(gene=lst[[x]],thr=0.,name=paste0("anterior_",x))}
        a <- a[order(a$enrichmentRatio,decreasing = T),]
        df <- data.frame(a[1:15,c(2,7,9)])
        df$cluster <- x
        df
}
))
df <- na.omit(data.frame(df))
df$cluster <- factor(df$cluster,levels=c("start","middle","end"))
ggplot(df,aes(x=cluster,y = factor(description,levels=rev(unique(description))), 
              color = -log10(FDR), size = enrichmentRatio)) + 
        geom_point() + scale_color_viridis_c(name = '-log10(FDR)') + 
        cowplot::theme_cowplot()+ylab("")+xlab("")

