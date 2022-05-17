setwd("E:/LabWork/Project/SCANSeq2/batch/6 IGG/")

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(data.table)

source("E:/LabWork/code/MyFunction.R")
ReadGeneInfo("hg38", classify.protein = T)
ReadGeneInfo("mm10")

# 1. summarize data -------------------------------------------------------
# merge sample information
# IGG6
sample.info.IGG6 <- read.table("data/IGG6/summary.txt",header = T, row.names = 1)

sample.info.IGG6$Library <- "IGG6"
sample.info.IGG6$Name <- rownames(sample.info.IGG6)
sample.info.IGG6$Rename <- paste("IGG6", sample.info.IGG6$Name, sep = "_")
sample.info.IGG6$Barcode.5 <- substr(sample.info.IGG6$Name,3,4)
sample.info.IGG6$Barcode.3 <- substr(sample.info.IGG6$Name,8,9)
rownames(sample.info.IGG6) <- sample.info.IGG6$Rename

sample.info.IGG6$Cell_Line <- "none"
sample.info.IGG6[sample.info.IGG6$Barcode.5 %in% as.character(57:68),]$Cell_Line <- "Hela"
sample.info.IGG6[sample.info.IGG6$Barcode.5 %in% as.character(69:80),]$Cell_Line <- "HepG2"
table(sample.info.IGG6$Cell_Line)

sample.info.IGG6$Concentration <- 0
sample.info.IGG6[sample.info.IGG6$Barcode.5 %in% as.character(c(57:62,
                                                                69:74)),]$Concentration <- 10
sample.info.IGG6[sample.info.IGG6$Barcode.5 %in% as.character(c(63:68,
                                                                75:80)),]$Concentration <- 30
table(sample.info.IGG6$Concentration)

sample.info.IGG6$Time <- 6

# IGG24
sample.info.IGG24 <- read.table("data/IGG24/summary.txt",header = T, row.names = 1)

sample.info.IGG24$Library <- "IGG24"
sample.info.IGG24$Name <- rownames(sample.info.IGG24)
sample.info.IGG24$Rename <- paste("IGG24", sample.info.IGG24$Name, sep = "_")
sample.info.IGG24$Barcode.5 <- substr(sample.info.IGG24$Name,3,4)
sample.info.IGG24$Barcode.3 <- substr(sample.info.IGG24$Name,8,9)
rownames(sample.info.IGG24) <- sample.info.IGG24$Rename

sample.info.IGG24$Cell_Line <- "none"
sample.info.IGG24[sample.info.IGG24$Barcode.5 %in% as.character(63:74),]$Cell_Line <- "Hela"
sample.info.IGG24[sample.info.IGG24$Barcode.5 %in% as.character(51:62),]$Cell_Line <- "HepG2"
table(sample.info.IGG24$Cell_Line)

sample.info.IGG24$Concentration <- 0
sample.info.IGG24[sample.info.IGG24$Barcode.5 %in% as.character(c(51:56,
                                                                  63:68)),]$Concentration <- 10
sample.info.IGG24[sample.info.IGG24$Barcode.5 %in% as.character(c(57:62,
                                                                  69:74)),]$Concentration <- 30
table(sample.info.IGG24$Concentration)

sample.info.IGG24$Time <- 24

# IGG48
sample.info.IGG48 <- read.table("data/IGG48/summary.txt",header = T, row.names = 1)

sample.info.IGG48$Library <- "IGG48"
sample.info.IGG48$Name <- rownames(sample.info.IGG48)
sample.info.IGG48$Rename <- paste("IGG48", sample.info.IGG48$Name, sep = "_")
sample.info.IGG48$Barcode.5 <- substr(sample.info.IGG48$Name,3,4)
sample.info.IGG48$Barcode.3 <- substr(sample.info.IGG48$Name,8,9)
rownames(sample.info.IGG48) <- sample.info.IGG48$Rename

sample.info.IGG48$Cell_Line <- "none"
sample.info.IGG48[sample.info.IGG48$Barcode.5 %in% as.character(c(51:53,78:80,
                                                                  66:74)),]$Cell_Line <- "Hela"
sample.info.IGG48[sample.info.IGG48$Barcode.5 %in% as.character(c(54:56,75:77,
                                                                  57:65)),]$Cell_Line <- "HepG2"
table(sample.info.IGG48$Cell_Line)

sample.info.IGG48$Concentration <- 0
sample.info.IGG48[sample.info.IGG48$Barcode.5 %in% as.character(c(60:62,
                                                                  69:71)),]$Concentration <- 10
sample.info.IGG48[sample.info.IGG48$Barcode.5 %in% as.character(c(63:65,
                                                                  72:74)),]$Concentration <- 30
table(sample.info.IGG48$Concentration)

sample.info.IGG48$Time <- 48
sample.info.IGG48[sample.info.IGG48$Barcode.5 %in% as.character(c(51:56)),]$Time <- 6
sample.info.IGG48[sample.info.IGG48$Barcode.5 %in% as.character(c(75:80)),]$Time <- 24
table(sample.info.IGG48$Time)

# merge
sample.info <- rbind(sample.info.IGG6,
                     sample.info.IGG24,
                     sample.info.IGG48)

sample.info$Cell_Line <- factor(sample.info$Cell_Line,
                                levels = c("HepG2","Hela"))
sample.info$Concentration <- factor(sample.info$Concentration, 
                                    levels = c(0,10,30))
sample.info$Time <- factor(sample.info$Time, 
                           levels = c(6,24,48))

# merge umi_count matrix
# IGG6
rc_gene_IGG6 <- as.data.frame(fread("data/IGG6/Gene.Count.txt",header = T))
rc_trans_IGG6 <- as.data.frame(fread("data/IGG6/Trans.Count.txt",header = T))
rownames(rc_gene_IGG6) <- rc_gene_IGG6$Name
rownames(rc_trans_IGG6) <- rc_trans_IGG6$Name
rc_gene_IGG6 <- rc_gene_IGG6[,-1]
rc_trans_IGG6 <- rc_trans_IGG6[,-1]

rc_gene_IGG6 <- rc_gene_IGG6[,sample.info.IGG6$Name]
rc_trans_IGG6 <- rc_trans_IGG6[,sample.info.IGG6$Name]

colnames(rc_gene_IGG6) <- sample.info.IGG6$Rename
colnames(rc_trans_IGG6) <- sample.info.IGG6$Rename

# IGG24
rc_gene_IGG24 <- as.data.frame(fread("data/IGG24/Gene.Count.txt",header = T))
rc_trans_IGG24 <- as.data.frame(fread("data/IGG24/Trans.Count.txt",header = T))
rownames(rc_gene_IGG24) <- rc_gene_IGG24$Name
rownames(rc_trans_IGG24) <- rc_trans_IGG24$Name
rc_gene_IGG24 <- rc_gene_IGG24[,-1]
rc_trans_IGG24 <- rc_trans_IGG24[,-1]

rc_gene_IGG24 <- rc_gene_IGG24[,sample.info.IGG24$Name]
rc_trans_IGG24 <- rc_trans_IGG24[,sample.info.IGG24$Name]

colnames(rc_gene_IGG24) <- sample.info.IGG24$Rename
colnames(rc_trans_IGG24) <- sample.info.IGG24$Rename

# IGG48
rc_gene_IGG48 <- as.data.frame(fread("data/IGG48/Gene.Count.txt",header = T))
rc_trans_IGG48 <- as.data.frame(fread("data/IGG48/Trans.Count.txt",header = T))
rownames(rc_gene_IGG48) <- rc_gene_IGG48$Name
rownames(rc_trans_IGG48) <- rc_trans_IGG48$Name
rc_gene_IGG48 <- rc_gene_IGG48[,-1]
rc_trans_IGG48 <- rc_trans_IGG48[,-1]

rc_gene_IGG48 <- rc_gene_IGG48[,sample.info.IGG48$Name]
rc_trans_IGG48 <- rc_trans_IGG48[,sample.info.IGG48$Name]

colnames(rc_gene_IGG48) <- sample.info.IGG48$Rename
colnames(rc_trans_IGG48) <- sample.info.IGG48$Rename

# merge
rc_gene <- cbind(rc_gene_IGG6, rc_gene_IGG24, rc_gene_IGG48)
rc_trans <- cbind(rc_trans_IGG6, rc_trans_IGG24, rc_trans_IGG48)
identical(rownames(sample.info), colnames(rc_gene))
identical(rownames(sample.info), colnames(rc_trans))
rm(sample.info.IGG6, sample.info.IGG24, sample.info.IGG48,
   rc_gene_IGG6, rc_gene_IGG24, rc_gene_IGG48,
   rc_trans_IGG6, rc_trans_IGG24, rc_trans_IGG48)
gc()

# remove ERCC rows
ERCC <- grep("ERCC|RGC",rownames(rc_gene),value = T)
rc_ERCC <- rc_gene[ERCC,]
rc_gene <- rc_gene[!rownames(rc_gene)%in%ERCC,]
rc_trans <- rc_trans[!rownames(rc_trans)%in%ERCC,]

human.gene <- grep("ENSG",rownames(rc_gene),value = T)
mouse.gene <- grep("ENSMU",rownames(rc_gene),value = T)
human.trans <- grep("ENST",rownames(rc_trans),value = T)
mouse.trans <- grep("ENSMU",rownames(rc_trans),value = T)

rc_gene <- rc_gene[human.gene,]
rc_trans <- rc_trans[human.trans,]

sample.info$Gene_Detected <- colSums(rc_gene[, rownames(sample.info)]>=1)
sample.info$Trans_Detected <- colSums(rc_trans[, rownames(sample.info)]>=1)

sample.info$Mapped_percent <- sample.info$Mapped_Reads / sample.info$Raw_Reads
sample.info$QC_Reads <- NULL

color.cell <- c("HepG2" = "#ffed6f",
                "Hela" = "#BEBADA")
color.time <- brewer.pal(3,"Set2")
color.concentration <- brewer.pal(3,"Set1")


# 3. QC of Cell -----------------------------------------------------------
# basic statistics
pdf("QC.stat.pdf",5,3)
ggplot(sample.info)+
  geom_violin(aes(x=Time,y=log10(Raw_Reads),fill=Cell_Line,alpha=Time),size=0.4,show.legend = F)+
  geom_boxplot(aes(x=Time,y=log10(Raw_Reads),color=Cell_Line),size=0.4,alpha=0,show.legend = F)+
  facet_grid(vars(Cell_Line), vars(Concentration))+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_alpha_manual(values = c(0.9, 0.6, 0.3))+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Time,y=log10(Mapped_Reads),fill=Cell_Line,alpha=Time),size=0.4,show.legend = F)+
  geom_boxplot(aes(x=Time,y=log10(Mapped_Reads),color=Cell_Line),size=0.4,alpha=0,show.legend = F)+
  facet_grid(vars(Cell_Line), vars(Concentration))+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_alpha_manual(values = c(0.9, 0.6, 0.3))+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Time,y=Mapped_percent,fill=Cell_Line,alpha=Time),size=0.4,show.legend = F)+
  geom_boxplot(aes(x=Time,y=Mapped_percent,color=Cell_Line),size=0.4,alpha=0,show.legend = F)+
  facet_grid(vars(Cell_Line), vars(Concentration))+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_alpha_manual(values = c(0.9, 0.6, 0.3))+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Time,y=Gene_Detected,fill=Cell_Line,alpha=Time),size=0.4,show.legend = F)+
  geom_boxplot(aes(x=Time,y=Gene_Detected,color=Cell_Line),size=0.4,alpha=0,show.legend = F)+
  facet_grid(vars(Cell_Line), vars(Concentration))+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_alpha_manual(values = c(0.9, 0.6, 0.3))+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Time,y=Trans_Detected,fill=Cell_Line,alpha=Time),size=0.4,show.legend = F)+
  geom_boxplot(aes(x=Time,y=Trans_Detected,color=Cell_Line),size=0.4,alpha=0,show.legend = F)+
  facet_grid(vars(Cell_Line), vars(Concentration))+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_alpha_manual(values = c(0.9, 0.6, 0.3))+
  theme_lzy()
dev.off()

pdf("Dot.Gene_Trans.pdf",5,4)
ggplot(sample.info)+
  geom_point(aes(x=Mapped_Reads, y=Gene_Detected, color=Time, shape=Cell_Line))+
  scale_color_manual(values = color.time)+
  theme_bw()+theme(panel.grid = element_blank())
dev.off()

aggregate(sample.info$Raw_Reads,list("Time" = sample.info$Time,
                                     "Cell" = sample.info$Cell_Line),FUN=median)
aggregate(sample.info$Mapped_percent,list("Time" = sample.info$Time,
                                          "Cell" = sample.info$Cell_Line),FUN=median)
aggregate(sample.info$Gene_Detected,list("Time" = sample.info$Time,
                                         "Cell" = sample.info$Cell_Line),FUN=median)
aggregate(sample.info$Trans_Detected,list("Time" = sample.info$Time,
                                          "Cell" = sample.info$Cell_Line),FUN=median)
aggregate(sample.info$Gene_Detected,list("Time" = sample.info$Time,
                                         "Cell" = sample.info$Cell_Line),FUN=mean)
aggregate(sample.info$Trans_Detected,list("Time" = sample.info$Time,
                                          "Cell" = sample.info$Cell_Line),FUN=mean)

# rename genes and trans
#human gene
human.gene <- substr(human.gene,1,15)
rownames(rc_gene) <- human.gene

rc_gene <- rc_gene[rownames(rc_gene)%in%gene.info.hg38$`Gene stable ID`,]
rc_gene <- rc_gene[!duplicated(gene.info.hg38[rownames(rc_gene),]$`Gene name`),]
rownames(rc_gene) <- gene.info.hg38[rownames(rc_gene),]$`Gene name`

#human trans
temp <- substr(rownames(rc_trans),1,15)
temp <- duplicated(temp)
rc_trans <- rc_trans[!temp,]
rownames(rc_trans) <- substr(rownames(rc_trans),1,15)
table(rownames(rc_trans) %in% gene2trans$trans.id)

rc_trans <- rc_trans[rownames(rc_trans) %in% gene2trans$trans.id,]
rownames(rc_trans) <- gene2trans[rownames(rc_trans),]$trans.symbol

saveRDS(rc_gene,"data/rc_gene.RDS")
saveRDS(rc_trans,"data/rc_trans.RDS")
saveRDS(rc_ERCC,"data/rc_ERCC.RDS")
gc()

# Create Object
mix.obj <- CreateSeuratObject(rc_gene,meta.data = sample.info,
                              project = "IGG",assay = "Gene")
mix.obj
Idents(mix.obj) <- mix.obj$Cell_Line

mix.obj$percent.mt <- PercentageFeatureSet(mix.obj, pattern = "^MT-|^mt-")
sample.info$percent.mt <- mix.obj$percent.mt

sample.info$QC_Density <- get_density(sample.info$Gene_Detected,
                                      sample.info$percent.mt,
                                      n = 100)

pdf("Scatter.QC.pdf",8,3)
p1 <- ggplot(sample.info, aes(x=Gene_Detected,y=percent.mt))+
  geom_point(aes(color = QC_Density),
             show.legend = T, size=0.01)+
  viridis::scale_color_viridis()+
  geom_vline(xintercept = 3500, color="blue",lty=2, size = 0.2)+
  geom_hline(yintercept = 15, color="blue",lty=2, size = 0.2)+
  theme_lzy()+coord_equal(210)
p2 <- ggplot(sample.info, aes(x=Gene_Detected,y=percent.mt))+
  geom_point(aes(color = Cell_Line),
             show.legend = T, size=0.01)+
  scale_color_manual(values = color.cell)+
  geom_vline(xintercept = 3500, color="blue",lty=2, size = 0.2)+
  geom_hline(yintercept = 15, color="blue",lty=2, size = 0.2)+
  theme_lzy()+coord_equal(210)
multiplot(p1,p2,cols = 2)
dev.off()

table(mix.obj$Gene_Detected> 3500 & mix.obj$percent.mt < 15)

mix.obj <- subset(mix.obj, subset = Gene_Detected > 3500 & percent.mt < 15)

sample.info$Pass_QC <- 0
sample.info[sample.info$Rename%in%colnames(mix.obj),]$Pass_QC <- 1

table(mix.obj$Time, mix.obj$Concentration,mix.obj$Cell_Line)

pdf("QC.stat.QC.pdf",5,3)
ggplot(sample.info[sample.info$Pass_QC==1,])+
  geom_violin(aes(x=Time,y=Gene_Detected,fill=Cell_Line,alpha=Time),size=0.4,show.legend = F)+
  geom_boxplot(aes(x=Time,y=Gene_Detected,color=Cell_Line),size=0.4,alpha=0,show.legend = F)+
  facet_grid(vars(Cell_Line), vars(Concentration))+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_alpha_manual(values = c(0.9, 0.6, 0.3))+
  theme_lzy()
ggplot(sample.info[sample.info$Pass_QC==1,])+
  geom_violin(aes(x=Time,y=Trans_Detected,fill=Cell_Line,alpha=Time),size=0.4,show.legend = F)+
  geom_boxplot(aes(x=Time,y=Trans_Detected,color=Cell_Line),size=0.4,alpha=0,show.legend = F)+
  facet_grid(vars(Cell_Line), vars(Concentration))+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_alpha_manual(values = c(0.9, 0.6, 0.3))+
  theme_lzy()
dev.off()

median(sample.info[sample.info$Pass_QC==1&
                     sample.info$Cell_Line=="HepG2",]$Gene_Detected) # 6126
median(sample.info[sample.info$Pass_QC==1&
                     sample.info$Cell_Line=="Hela",]$Trans_Detected) # 8369.5

# 4. gene-level analysis -----------------------------------------------------
mix.obj <- NormalizeData(mix.obj, normalization.method = "LogNormalize",
                             scale.factor = 10000)

mix.obj <- FindVariableFeatures(mix.obj, selection.method = "vst", nfeatures = 3000)

#VariableFeaturePlot(mix.obj)

mix.obj <- ScaleData(mix.obj)

mix.obj <- RunPCA(mix.obj,npcs = 30,
                      features = VariableFeatures(object = mix.obj))

pdf("PCA.CellType.pdf",4.8,3.5)
DimPlot(mix.obj, reduction = "pca",label = T,cols = color.cell,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.3)
DimPlot(mix.obj, reduction = "pca",label = F,cols = color.cell,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.3)
dev.off()

#DimHeatmap(mix.obj, dims = 1:10, balanced = TRUE)
ElbowPlot(mix.obj) #10~15

pc.use=1:4

# UMAP
mix.obj <- RunUMAP(mix.obj, dims = pc.use,
                   n.neighbors = 50, min.dist = 0.1)

pdf("UMAP.CellType.pdf",4.8,3.5)
DimPlot(mix.obj, reduction = "umap",group.by = "Cell_Line",
        label = T,cols = color.cell,pt.size = 0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
DimPlot(mix.obj, reduction = "umap",group.by = "Concentration",
        label = F,cols = color.concentration,pt.size = 0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
DimPlot(mix.obj, reduction = "umap",group.by = "Time",
        label = T,cols = color.time,pt.size = 0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
dev.off()

pdf("UMAP.Quality.pdf",4.5,3.5)
FeaturePlot(mix.obj,"nFeature_Gene", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
FeaturePlot(mix.obj,"nCount_Gene", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
FeaturePlot(mix.obj,"Mapped_Reads", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
FeaturePlot(mix.obj,"Mapped_percent", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
FeaturePlot(mix.obj,"percent.mt", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1.4)
dev.off()

mix.obj <- FindNeighbors(mix.obj,reduction = "pca",dims = pc.use,k.param = 20)
mix.obj <- FindClusters(mix.obj,resolution = 0.2)

pdf("UMAP.Cluster.pdf",4.4,3.5)
DimPlot(mix.obj, reduction = "umap",label = T,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values = brewer.pal(10,"Set3"))+
  coord_equal(1.4)
DimPlot(mix.obj, reduction = "umap",label = F,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values = brewer.pal(10,"Set3"))+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  coord_equal(1.4)
dev.off()

# 7. ERCC ---------------------------------------------------------------
# read ERCC info
ERCC.info <- as.data.frame(fread("E:/LabWork/genome/Thermo_ERCC.txt"))
ERCC.info <- ERCC.info[,c(1:5)]
colnames(ERCC.info) <- c("order","id","group","concentration_mix_1",
                         "concentration_mix_2")
rc_ERCC <- rc_ERCC[ERCC.info$id,]

ggplot()+
  geom_point(aes(x=log10(1+ERCC.info$concentration_mix_1),
                 y=log10(1+rc_ERCC$IGG6_Bc57_Bc01),color="IGG6_Bc57_Bc01"))+
  geom_point(aes(x=log10(1+ERCC.info$concentration_mix_1),
                 y=log10(1+rc_ERCC$IGG6_Bc61_Bc15),color="IGG6_Bc61_Bc15"))+
  geom_point(aes(x=log10(1+ERCC.info$concentration_mix_1),
                 y=log10(1+rc_ERCC$IGG6_Bc57_Bc24),color="IGG6_Bc57_Bc24"))+
  geom_point(aes(x=log10(1+ERCC.info$concentration_mix_1),
                 y=log10(1+rc_ERCC$IGG6_Bc59_Bc31),color="IGG6_Bc59_Bc31"))+
  xlab("log10(ERCC concentration (attomoles/ul))")+
  ylab("log(UMI Count + 1)")+
  coord_fixed()+
  theme_lzy()

sample.info$cor_ERCC_rc <- 0
for (cell in rownames(sample.info)) {
  sample.info[cell,]$cor_ERCC_rc <- cor(rc_ERCC[,cell],
                                        ERCC.info$concentration_mix_1)
}

pdf("Density.cor.ERCC.pdf",4.5,3)
ggplot(sample.info[sample.info$Pass_QC==1,])+
  geom_density(aes(x = cor_ERCC_rc, color = Cell_Line), alpha = 0.5)+
  scale_color_manual(values = color.cell)+
  theme_lzy()
dev.off()

# write sample info
write.table(sample.info,"sample.info.txt",quote = F,sep = "\t")
saveRDS(mix.obj,"obj.IGG.rds")
