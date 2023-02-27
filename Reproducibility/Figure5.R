library(ggplot2)
library(ggrepel)
set.seed(1234)

# load Cellrank+WOT driver genes ####
res <- read.delim("python/CellRank+WOT_CNC_0617_no_chondrocyte/Cellrank+WOT_driver_genes.txt",head=T,sep=",",row.names = 1)
res <- data.frame(res)
res <- na.omit(res[,c(1,2,6,7,11,12,16,17,21,22)])
# res$peak <- apply(res,1,function(x){
#         names(which.max(x[c(1,3,5,7,9)]))
# })
corr.thres <- 0.05
# anterior
anterior <- na.omit(res[order(res$anterior.palatal.mesenchymal_corr,decreasing = T),])
anterior$adjusted_pval <- p.adjust(anterior$anterior.palatal.mesenchymal_pval,method="BH")
anterior <- anterior[anterior$adjusted_pval < 0.05 & anterior$anterior.palatal.mesenchymal_corr > corr.thres,]

# posterior
posterior <- na.omit(res[order(res$posterior.palatal.mesenchymal_corr,decreasing = T),])
posterior$adjusted_pval <- p.adjust(posterior$posterior.palatal.mesenchymal_pval,method="BH")
posterior <- posterior[posterior$adjusted_pval < 0.05 & posterior$posterior.palatal.mesenchymal_corr > corr.thres,]

# osteoblast
osteoblast <- na.omit(res[order(res$Osteoblast_corr,decreasing = T),])
osteoblast$adjusted_pval <- p.adjust(osteoblast$Osteoblast_pval,method="BH")
osteoblast <- osteoblast[osteoblast$adjusted_pval < 0.05 & osteoblast$Osteoblast_corr > corr.thres,]

# dental
dental <- na.omit(res[order(res$dental.mesenchymal_corr,decreasing = T),])
dental$adjusted_pval <- p.adjust(dental$dental.mesenchymal_pval,method="BH")
dental <- dental[dental$adjusted_pval < 0.05 & dental$dental.mesenchymal_corr > corr.thres,]

# perimysial 
perimysial <- na.omit(res[order(res$perimysial_corr,decreasing = T),])
perimysial$adjusted_pval <- p.adjust(perimysial$perimysial_pval,method="BH")
perimysial <- perimysial[perimysial$adjusted_pval < 0.05 & perimysial$perimysial_corr > corr.thres,]

# Scatter plot showing driver genes for anterior vs posterior ####
res$driver_genes <- "non signifiacnt"
res$driver_genes[res$posterior.palatal.mesenchymal_corr>0.05 & res$posterior.palatal.mesenchymal_pval<0.01 & 
                         res$anterior.palatal.mesenchymal_corr < -0.05] <- "driver genes for posterior"
res$driver_genes[res$anterior.palatal.mesenchymal_corr>0.05 & res$anterior.palatal.mesenchymal_pval<0.01 &
                         res$posterior.palatal.mesenchymal_corr < -0.05] <- "driver genes for anterior"
res$label <- rownames(res)
res$label[!(rownames(res) %in% c(rownames(anterior)[1:5],rownames(posterior)[1:5]))] <- NA
ggplot(res,aes(anterior.palatal.mesenchymal_corr,posterior.palatal.mesenchymal_corr,color=driver_genes,label=label))+geom_point()+
        theme_bw()+xlab("fate correlations to anterior")+ylab("fate correlations to posterior")+
        ggtitle("Driver genes for anterior vs posterior")+
        scale_color_manual(values = c("#FB766D","#00BFC4","grey"))+
        geom_vline(xintercept=0,color='lightgrey',linetype='dashed')+
        geom_hline(yintercept=0,color='lightgrey',linetype='dashed')+
        geom_label_repel(size=3,max.overlaps = 10)
# dental vs osteoblast
res$driver_genes <- "non signifiacnt"
res$driver_genes[res$dental.mesenchymal_corr > 0.05 & res$dental.mesenchymal_pval < 0.01] <- "driver genes for dental mesenchymal"
res$driver_genes[res$Osteoblast_corr > 0.05 & res$Osteoblast_pval < 0.01] <- "driver genes for osteoblast"
res$driver_genes[res$dental.mesenchymal_corr > 0.05 & res$Osteoblast_corr>0.05] <- "driver genes for both"
both <- res[res$Osteoblast_corr > 0.1 & res$dental.mesenchymal_corr > 0.2,]
res$label <- rownames(res)
res$label[!(rownames(res) %in% c(rownames(osteoblast)[1:5],rownames(dental)[1:5],
                                 rownames(both)[2:5]))] <- NA
ggplot(res,aes(Osteoblast_corr,dental.mesenchymal_corr,color=driver_genes,label=label))+geom_point()+
        theme_bw()+xlab("fate correlations to osteoblast")+ylab("fate correlations to dental mesenchymal")+
        ggtitle("Driver genes for osteoblast vs dental mesenchymal")+
        scale_color_manual(values = c("#7CAE00","#FB766D","#00BFC4","grey"))+
        geom_vline(xintercept=0,color='lightgrey',linetype='dashed')+
        geom_hline(yintercept=0,color='lightgrey',linetype='dashed')+
        geom_label_repel(size=3,max.overlaps = 10)

# log odds of fate probability ####
library(Seurat)
seu <- readRDS("data/subset_no_chondrocyte_0617.rds")
# Load absorption rate ####
prob <- read.delim("python/CellRank+WOT_CNC_0617_no_chondrocyte/cellrank_WOT_terminal_prob.txt",
                   head=T,sep=",",row.names = 1)
df <- cbind(seu@meta.data,prob[colnames(seu),])
df$log1 <- log(df$anterior.palatal.mesenchymal/df$posterior.palatal.mesenchymal)
df$log2 <- log(df$posterior.palatal.mesenchymal/df$anterior.palatal.mesenchymal)
DefaultAssay(seu) <- "RNA"
seu <- ScaleData(seu,features = c("Shox2",'Meox2'),vars.to.regress = c("S.Score", "G2M.Score"))
df <- cbind(df,t(seu@assays$RNA@scale.data))
p1 <- ggplot(df,aes(stage,log1,color=Shox2>0))+
        geom_point(position = position_jitter(seed = 1, width = 0.1),size=1)+
        theme_bw()+xlab("Developmental Stage")+ylab("log(anterior/posterior)")+
        ggtitle("Shox2>0")+scale_color_manual(values=c("grey","steelblue"))
p2 <- ggplot(df,aes(stage,log1,color=celltype))+
        geom_point(position = position_jitter(seed = 1, width = 0.1),size=1)+
        theme_bw()+xlab("Developmental Stage")+ylab("log(anterior/posterior)")+
        ggtitle("Cell states")
p1|p2

p1 <- ggplot(df,aes(stage,log2,color=Meox2>0))+
        geom_point(position = position_jitter(seed = 1, width = 0.1),size=1)+
        theme_bw()+xlab("Developmental Stage")+ylab("log(posterior/anterior)")+
        ggtitle("Meox2>0")+scale_color_manual(values=c("grey","steelblue"))
p2 <- ggplot(df,aes(stage,log2,color=celltype))+
        geom_point(position = position_jitter(seed = 1, width = 0.1),size=1)+
        theme_bw()+xlab("Developmental Stage")+ylab("log(posterior/anterior)")+
        ggtitle("Cell states")
p1|p2



