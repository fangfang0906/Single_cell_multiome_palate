library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
set.seed(1234)
seu <- readRDS('data/seu.rds')
time_palette <- RColorBrewer::brewer.pal(9,"BuGn")[c(4,6,7,9)]

# Figure 1B ####
p1 <- DimPlot(seu,group.by = 'stage',reduction='rna.umap',label=F,cols=time_palette)+
        labs(title="snRNA-seq")
p2 <- DimPlot(seu,group.by = 'stage',reduction='atac.umap',label=F,cols=time_palette)+
        labs(title='snATAC-seq')
p3 <- DimPlot(seu,group.by = 'celltype_major',reduction='rna.umap')+
        labs(title='Cell type')
p4 <- DimPlot(seu,group.by = 'celltype_major',reduction='atac.umap')+
        labs(title='Cell type')
p1/p2
p3|p4

# Figure 1C ####
DefaultAssay(seu) <- "SCT"
Idents(seu) <- seu$celltype_major
m <- FindAllMarkers(seu,only.pos = T,min.pct = 0.2,min.diff.pct=0.2,logfc.threshold = 0.25,
                    max.cells.per.ident = 300)
m <- m[m$pct.1>m$pct.2 & m$p_val_adj<0.001,]
top <- m %>% group_by(cluster) %>% top_n(3, avg_log2FC)
p1 <- DotPlot(seu, features = unique(top$gene)) + RotatedAxis() + 
        theme(axis.title.y = element_blank(),axis.title.x=element_blank(),
              axis.text.x = element_text(size = 10,angle = 30, vjust = 0.8, hjust=0.5),
              axis.text.y = element_text(size = 10))+
        scale_color_viridis_c()+ggtitle("Gene expression (snRNA-seq)")
DefaultAssay(seu) <- "Gene_activity"
p2 <- DotPlot(seu, features = unique(top$gene)) + RotatedAxis() + 
        theme(axis.title.y = element_blank(),axis.title.x=element_blank(),
              axis.text.x = element_text(size = 10,angle = 45, vjust = 0.8, hjust=0.5),
              axis.text.y = element_text(size = 10))+
        scale_color_viridis_c()+ggtitle("Gene activity (snATAC-seq)")
p1|p2

# Figure 1E (cell composition plot) ####
df <- seu@meta.data %>% group_by(stage,library,celltype_major) %>%
        summarise(n = n()) %>% mutate(freq = n / sum(n))
df <- as.data.frame(df)
df$celltype <- factor(df$celltype)
ggplot(df, aes(x=library,y=freq, fill=celltype)) + 
        geom_bar(position="stack", stat="identity")+
        ggtitle("Cell type frequencies")+
        xlab("")+ylab("% of cells")+theme_bw()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(df, aes(x=stage,y=n, fill=celltype)) + 
        geom_bar(position="stack", stat="identity")+
        ggtitle("Cell type frequencies")+
        xlab("")+ylab("# of cells")+theme_bw()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
m <- table(seu$celltype_major,seu$stage)
round(prop.table(m,margin = 2),2)










