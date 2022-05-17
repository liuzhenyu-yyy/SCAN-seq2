setwd("E:/LabWork/Project/SCANSeq2/batch/1 9CL/")

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
sample.info <- read.table("data/summary.txt",header = T, row.names = 1)
sample.info$Mapped_percent <- sample.info$Mapped_Reads / sample.info$Raw_Reads
sample.info$QC_Reads <- NULL

sample.info$Library <- "9CL"
sample.info$Name <- rownames(sample.info)
sample.info$Rename <- paste("9CL", sample.info$Name, sep = "_")
sample.info$Barcode.5 <- substr(sample.info$Name,3,4)
sample.info$Barcode.3 <- substr(sample.info$Name,8,9)
rownames(sample.info) <- sample.info$Rename

sample.info$Cell_Line <- "none"
sample.info[sample.info$Barcode.5 %in% as.character(53:55),]$Cell_Line <- "AtT20"
sample.info[sample.info$Barcode.5 %in% as.character(56:58),]$Cell_Line <- "MEF"
sample.info[sample.info$Barcode.5 %in% as.character(59:61),]$Cell_Line <- "3T3"
sample.info[sample.info$Barcode.5 %in% as.character(62:64),]$Cell_Line <- "Hela"
sample.info[sample.info$Barcode.5 %in% as.character(65:67),]$Cell_Line <- "HepG2"
sample.info[sample.info$Barcode.5 %in% as.character(68:70),]$Cell_Line <- "293T"
sample.info[sample.info$Barcode.5 %in% as.character(71:73),]$Cell_Line <- "H9"
sample.info[sample.info$Barcode.5 %in% as.character(74:79),]$Cell_Line <- "K562"
sample.info[sample.info$Barcode.5 %in% as.character(80:82),]$Cell_Line <- "GM12878"
table(sample.info$Cell_Line)

sample.info$Cell_Line <- factor(sample.info$Cell_Line,
                                levels = c("GM12878","HepG2","Hela","H9","K562","293T",
                                           "AtT20","MEF","3T3"))

sample.info$Organism <- ifelse(sample.info$Cell_Line %in% c("GM12878","HepG2","Hela","H9","K562","293T"),
                               "Human","Mouse")
table(sample.info$Organism)

# merge umi_count matrix
rc_gene <- as.data.frame(fread("data/Gene.Count.txt",header = T))
rc_trans <- as.data.frame(fread("data/Trans.Count.txt",header = T))
rownames(rc_gene) <- rc_gene$Name
rownames(rc_trans) <- rc_trans$Name
rc_gene <- rc_gene[,-1]
rc_trans <- rc_trans[,-1]

rc_gene <- rc_gene[,sample.info$Name]
rc_trans <- rc_trans[,sample.info$Name]

colnames(rc_gene) <- sample.info$Rename
colnames(rc_trans) <- sample.info$Rename

gc()

# remove ERCC rows
ERCC <- grep("ERCC|RGC",rownames(rc_gene),value = T)
rc_ERCC <- rc_gene[ERCC,]
rc_gene <- rc_gene[!rownames(rc_gene)%in%ERCC,]
rc_trans <- rc_trans[!rownames(rc_trans)%in%ERCC,]

human.cell <- sample.info[sample.info$Organism=="Human",]$Rename
mouse.cell <- sample.info[sample.info$Organism=="Mouse",]$Rename
human.gene <- grep("ENSG",rownames(rc_gene),value = T)
mouse.gene <- grep("ENSMU",rownames(rc_gene),value = T)
human.trans <- grep("ENST",rownames(rc_trans),value = T)
mouse.trans <- grep("ENSMU",rownames(rc_trans),value = T)

rc_gene_human <- rc_gene[human.gene,]
rc_trans_human <- rc_trans[human.trans,]
rc_gene_mouse <- rc_gene[mouse.gene,]
rc_trans_mouse <- rc_trans[mouse.trans,]
rm(rc_gene, rc_trans)

sample.info.human <- sample.info[human.cell,]
sample.info.mouse <- sample.info[mouse.cell,]

sample.info.human$Gene_Detected <- colSums(rc_gene_human[, human.cell]>=1)
sample.info.human$Trans_Detected <- colSums(rc_trans_human[, human.cell]>=1)
sample.info.mouse$Gene_Detected <- colSums(rc_gene_mouse[, mouse.cell]>=1)
sample.info.mouse$Trans_Detected <- colSums(rc_trans_mouse[, mouse.cell]>=1)

sample.info <- rbind(sample.info.human, sample.info.mouse)

mycolor <- brewer.pal(11,"Set3")[c(1:8,11)]
names(mycolor) <- c("GM12878","HepG2","Hela","H9","K562","293T",
                    "AtT20","MEF","3T3")
# mycolor["3T3"] <- "gray"
mycolor["HepG2"] <- "#ffed6f"

# 2. demultiplex report -------------------------------------------------------
TotalReads <- 84801426

plot.data <- aggregate(sample.info$Raw_Reads,
                       list(sample.info$Cell_Line),
                       FUN = sum)
colnames(plot.data) <- c("Cell_Line", "Reads")
temp <- t(as.data.frame(c("unclassified",
          TotalReads-sum(plot.data$Reads))))
colnames(temp) <- c("Cell_Line", "Reads")
plot.data <- rbind(plot.data, temp)
rownames(plot.data) <- plot.data$Cell_Line
plot.data$Freq <- as.numeric(plot.data$Reads) / TotalReads
plot.data$label <- paste(plot.data$Cell_Line, "\n",
                         round(plot.data$Freq,3)*100, "%", sep = "")

pdf("Pie.demultiplex.pdf",6,6)
ggpubr::ggdonutchart(plot.data, "Freq",
                     label = "label",                               
                     fill = "Cell_Line",                            
                     color = "white",
                     lab.pos = "in",
                     palette = c(mycolor,
                                 "unclassified" = "gray75"))
dev.off()  

mean(sample.info$Raw_Reads) # 67025.06
quantile(sample.info$Raw_Reads) # 62712.00   

sum(sample.info$Mapped_Reads) / sum(sample.info$Raw_Reads)
quantile(sample.info$Mapped_Reads / sample.info$Raw_Reads)

# 3. QC of Cell -----------------------------------------------------------
# basic statistics
pdf("Vln.QC.stat.pdf",6,2)
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=log10(Raw_Reads),fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=log10(Raw_Reads),color=Cell_Line),alpha=0,show.legend = F)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  theme_bw()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=log10(Mapped_Reads),fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=log10(Mapped_Reads),color=Cell_Line),alpha=0,show.legend = F)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  theme_bw()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Mapped_percent,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Mapped_percent,color=Cell_Line),alpha=0,show.legend = F)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  theme_bw()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Gene_Detected,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Gene_Detected,color=Cell_Line),alpha=0,show.legend = F)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  theme_bw()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Trans_Detected,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Trans_Detected,color=Cell_Line),alpha=0,show.legend = F)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  theme_bw() 
dev.off()

pdf("Dot.Gene_Trans.pdf",5,3.5)
ggplot(sample.info)+
  geom_point(aes(x=Mapped_Reads, y=Gene_Detected, color=Cell_Line, shape=Organism))+
  scale_color_manual(values = mycolor)+
  theme_bw()+theme(panel.grid = element_blank())
dev.off()

quantile(sample.info$Mapped_Reads) # 52374.50         
median(sample.info[sample.info$Organism=="Human",]$Gene_Detected) # 5483.5
median(sample.info[sample.info$Organism=="Mouse",]$Trans_Detected) # 5251

mean(sample.info[sample.info$Organism=="Human",]$Gene_Detected) # 5053.818
mean(sample.info[sample.info$Organism=="Mouse",]$Trans_Detected) # 5076.059

aggregate(sample.info$Raw_Reads,list(sample.info$Cell_Line), FUN=mean)
aggregate(sample.info$Mapped_percent,list(sample.info$Cell_Line), FUN=mean)
aggregate(sample.info$Gene_Detected,list(sample.info$Cell_Line), FUN=mean)
aggregate(sample.info$Trans_Detected,list(sample.info$Cell_Line), FUN=mean)

# rename genes and trans
#human gene
human.gene <- substr(human.gene,1,15)
rownames(rc_gene_human) <- human.gene

rc_gene_human <- rc_gene_human[rownames(rc_gene_human)%in%gene.info.hg38$`Gene stable ID`,]
rc_gene_human <- rc_gene_human[!duplicated(gene.info.hg38[rownames(rc_gene_human),]$`Gene name`),]
rownames(rc_gene_human) <- gene.info.hg38[rownames(rc_gene_human),]$`Gene name`

#human trans
temp <- substr(rownames(rc_trans_human),1,15)
temp <- duplicated(temp)
rc_trans_human <- rc_trans_human[!temp,]
rownames(rc_trans_human) <- substr(rownames(rc_trans_human),1,15)
table(rownames(rc_trans_human) %in% gene2trans$trans.id)

rc_trans_human <- rc_trans_human[rownames(rc_trans_human) %in% gene2trans$trans.id,]
rownames(rc_trans_human) <- gene2trans[rownames(rc_trans_human),]$trans.symbol

#mouse genes
mouse.gene <- substr(mouse.gene,1,18)
rownames(rc_gene_mouse) <- mouse.gene

rc_gene_mouse <- rc_gene_mouse[rownames(rc_gene_mouse)%in%gene.info.mm10$`Gene stable ID`,]
rc_gene_mouse <- rc_gene_mouse[!duplicated(gene.info.mm10[rownames(rc_gene_mouse),]$`Gene name`),]
rownames(rc_gene_mouse) <- gene.info.mm10[rownames(rc_gene_mouse),]$`Gene name`

#merge
identical(colnames(rc_gene_mouse),colnames(rc_gene_human))
identical(colnames(rc_trans_mouse),colnames(rc_trans_human))
rc_gene <- rbind(rc_gene_human, rc_gene_mouse)
rc_trans <- rbind(rc_trans_human, rc_trans_mouse)
rm(rc_gene_mouse, rc_gene_human,
   rc_trans_human, rc_trans_mouse)
saveRDS(rc_gene,"data/rc_gene.RDS")
saveRDS(rc_trans,"data/rc_trans.RDS")
saveRDS(rc_ERCC,"data/rc_ERCC.RDS")
saveRDS(sample.info,"data/sample.info.RDS")
gc()

# Create Object
mix.obj <- CreateSeuratObject(rc_gene,meta.data = sample.info,
                              project = "9CL",assay = "Gene")
mix.obj
Idents(mix.obj) <- mix.obj$Cell_Line

mix.obj$percent.mt <- PercentageFeatureSet(mix.obj, pattern = "^MT-|^mt-")

pdf("Scatter.QC.pdf",5,4)
FeatureScatter(mix.obj, feature1 = "Gene_Detected", feature2 = "percent.mt",cols = mycolor)+
  geom_hline(yintercept = 15, lty=2, color = "red")+
  geom_vline(xintercept = 2000, lty=2, color = "red")
dev.off()

table(mix.obj$Gene_Detected> 2000 & mix.obj$percent.mt < 15)

mix.obj <- subset(mix.obj, subset = Gene_Detected > 2000 & percent.mt < 15)

QC.g2k <- table(mix.obj$Cell_Line)
table(sample.info$Cell_Line)
QC.g2k / table(sample.info$Cell_Line)

sample.info$Pass_QC <- 0
sample.info[sample.info$Rename%in%colnames(mix.obj),]$Pass_QC <- 1

pdf("Vln.gene_trans.QC.pdf",6,2)
ggplot(sample.info[sample.info$Pass_QC==1,])+
  geom_violin(aes(x=Cell_Line,y=(Gene_Detected),fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=(Gene_Detected),color=Cell_Line),alpha=0,show.legend = F)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  theme_bw()
ggplot(sample.info[sample.info$Pass_QC==1,])+
  geom_violin(aes(x=Cell_Line,y=(Trans_Detected),fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=(Trans_Detected),color=Cell_Line),alpha=0,show.legend = F)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  theme_bw() 
dev.off()

median(sample.info[sample.info$Pass_QC==1&
                   sample.info$Organism=="Human",]$Gene_Detected) # 5674
median(sample.info[sample.info$Pass_QC==1&
                   sample.info$Organism=="Human",]$Trans_Detected) # 6695

median(sample.info[sample.info$Pass_QC==1&
                     sample.info$Organism=="Mouse",]$Gene_Detected) # 4901
median(sample.info[sample.info$Pass_QC==1&
                     sample.info$Organism=="Mouse",]$Trans_Detected) # 5418

# 3. Full length ----------------------------------------------------------
read.length <- fread("data/temp.100w.fa.len",header = F,sep = "\t")
hist(read.length$V2)
mean(read.length$V2)
pdf("Hist.length.pdf",4,3)
ggplot(read.length)+
  geom_histogram(aes(x=V2),fill="#92c5de",alpha = 1,color = "black",size=0.5,
                 binwidth = 100)+
  xlim(c(0, 4000))+
  theme_lzy()
dev.off()

# 4. gene-level analysis -----------------------------------------------------
mix.obj <- NormalizeData(mix.obj, normalization.method = "LogNormalize",
                             scale.factor = 10000)

mix.obj <- FindVariableFeatures(mix.obj, selection.method = "vst", nfeatures = 3000)

#VariableFeaturePlot(mix.obj)

mix.obj <- ScaleData(mix.obj)

mix.obj <- RunPCA(mix.obj,npcs = 30,
                      features = VariableFeatures(object = mix.obj))

pdf("PCA.CellType.pdf",4.8,3.5)
DimPlot(mix.obj, reduction = "pca",label = T,cols = mycolor,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(0.9)
DimPlot(mix.obj, reduction = "pca",label = F,cols = mycolor,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(0.9)
dev.off()

#DimHeatmap(mix.obj, dims = 1:10, balanced = TRUE)
ElbowPlot(mix.obj) #10~15

pc.use=1:8

# UMAP
mix.obj <- RunUMAP(mix.obj, dims = pc.use,
                   n.neighbors = 50, min.dist = 0.1)

pdf("UMAP.CellType.pdf",4.8,3.5)
DimPlot(mix.obj, reduction = "umap",group.by = "Cell_Line",
        label = T,cols = mycolor,pt.size = 0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1)
DimPlot(mix.obj, reduction = "umap",group.by = "Cell_Line",
        label = F,cols = mycolor,pt.size = 0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  coord_equal(1)+NoLegend()
dev.off()

pdf("UMAP.Quality.pdf",4.5,3.5)
FeaturePlot(mix.obj,"nFeature_Gene", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1)
FeaturePlot(mix.obj,"nCount_Gene", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1)
FeaturePlot(mix.obj,"Mapped_Reads", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1)
FeaturePlot(mix.obj,"Mapped_percent", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1)
FeaturePlot(mix.obj,"percent.mt", reduction = "umap", pt.size = 0.6,
            cols = c("gray","red"))+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1)
dev.off()

mix.obj <- FindNeighbors(mix.obj,reduction = "pca",dims = pc.use,k.param = 20)
mix.obj <- FindClusters(mix.obj,resolution = 0.2)

pdf("UMAP.Cluster.pdf",4.4,3.5)
DimPlot(mix.obj, reduction = "umap",label = T,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values = brewer.pal(10,"Set3"))+
  coord_equal(1)
DimPlot(mix.obj, reduction = "umap",label = F,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values = brewer.pal(10,"Set3"))+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  coord_equal(1)
dev.off()

pdf("UMAP.Cell_Type_Cluster.pdf",8,4)
p1 <- DimPlot(mix.obj, reduction = "umap",group.by = "Cell_Line",
        label = T,cols = mycolor,pt.size = 0.3)+
  theme_bw()+theme(panel.grid = element_blank())+
  coord_equal(1)+NoLegend()
p2 <- DimPlot(mix.obj, reduction = "umap",group.by = "Gene_snn_res.0.2",
              label = T,pt.size = 0.6)+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values = brewer.pal(10,"Set3"))+
  coord_equal(1)+NoLegend()
multiplot(p1,p2,cols = 2)
dev.off()

# marker heatmap
Idents(mix.obj) <- mix.obj$Cell_Line
marker.all <- FindAllMarkers(mix.obj, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)
table(marker.all$cluster)

top10 <- marker.all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
table(top10$cluster)

pdf("Heatmap.Marker.all.pdf",6,5)
mix.obj <- ScaleData(mix.obj,features =  top10$gene)
DoHeatmap(mix.obj, features = top10$gene,
          draw.lines = F,
          group.colors = mycolor)+
  scale_fill_gradientn(colors = rev(brewer.pal(7,"RdYlBu")))
mix.obj <- ScaleData(mix.obj)
dev.off()

top5 <- marker.all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
marker.list <- c(top5$gene)

pdf("UMAP.marker.pdf",4,3)
for (one in marker.list) {
  p <- FeaturePlot(mix.obj,one,order=F, pt.size = 1)+
    theme_lzy()+
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())+
    theme(plot.title = element_text(face = "bold.italic",size = 20,
                                    hjust = 0.5))+
    coord_equal(1)
  plot(p)
}
one <- "Scg5"
FeaturePlot(mix.obj,one,order=F, pt.size = 1)+
  theme_lzy()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  theme(plot.title = element_text(face = "bold.italic",size = 20,
                                  hjust = 0.5))+
  coord_equal(1)
dev.off()

pdf("Vln.marker.cluster.pdf",10,8)
VlnPlot(mix.obj,features = c("IGLC2","ALB","PRSS56",
                             "CD3G","HBA1","UCHL1",
                             "Cep85","Col1a1","Eno1b"),cols = mycolor,ncol = 3)
dev.off()

# 5. Barnyard plot --------------------------------------------------------
mouse.trans <- rownames(rc_trans)[grep("ENSMUST",rownames(rc_trans))]
human.trans <- rownames(rc_trans)[-grep("ENSMUST",rownames(rc_trans))]
nrow(rc_trans) - length(mouse.trans) - length(human.trans)

Barnyard.data <- data.frame(row.names = colnames(rc_trans),
                            "cell"=colnames(rc_trans),
                            "human_transcripts"=colSums(rc_trans[human.trans,]),
                            "mouse_transcripts"=colSums(rc_trans[mouse.trans,]))

Barnyard.data <- Barnyard.data[rownames(sample.info),]

rm(mouse.trans,human.trans)

Barnyard.data$Cell_Line <- sample.info$Cell_Line

Barnyard.data$ratio_human <- Barnyard.data$human_transcripts/(Barnyard.data$human_transcripts+
                                                                Barnyard.data$mouse_transcripts)
Barnyard.data$ratio_mouse <- Barnyard.data$mouse_transcripts/(Barnyard.data$human_transcripts+
                                                               Barnyard.data$mouse_transcripts)

Barnyard.data$ratio_human %>% hist()

Barnyard.data$Organism <- "Mixed"
Barnyard.data[Barnyard.data$ratio_human > 0.9,]$Organism <- "Human"
Barnyard.data[Barnyard.data$ratio_mouse > 0.9,]$Organism <- "Mouse"
table(Barnyard.data$Organism)
table(Barnyard.data$Organism) / nrow(Barnyard.data)

Barnyard.data.sub <- Barnyard.data[colnames(mix.obj),]
table(Barnyard.data.sub$Organism) 
table(Barnyard.data.sub$Organism, Barnyard.data.sub$Cell_Line) 
table(Barnyard.data.sub$Organism) / nrow(Barnyard.data.sub)

pdf("Hist_Ratio_Human.cDNA.pdf",8,6)
ggplot(Barnyard.data.sub)+
  geom_histogram(aes(x =ratio_human,fill=Cell_Line),binwidth = 0.01)+
  geom_vline(xintercept = 0.1,lty=2)+
  geom_vline(xintercept = 0.9,lty=2)+
  facet_wrap(~Cell_Line)+
  theme_bw()+xlab("Percent of Human UMIs")+
  scale_fill_manual(values = mycolor)
dev.off()


pdf("Barnyard.cDNA.pdf",4,3)
ggplot(Barnyard.data.sub)+
  geom_point(aes(x=human_transcripts/1000,y=mouse_transcripts/1000,color=Organism),
             size=1)+
  theme_bw()+
  xlab("Human UMI counts (k)")+
  ylab("Mouse UMI counts (k)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_equal(2)+
  scale_color_manual(values = c("#d8756e","#808183","#6dbbc1"))+
  geom_abline(slope = 1/9,lty=2)+geom_abline(slope = 9,lty=2)
ggplot(Barnyard.data)+
  geom_point(aes(x=human_transcripts/1000,y=mouse_transcripts/1000,color=Cell_Line),
             size=0.3)+
  theme_bw()+
  scale_color_manual(values = mycolor)+
  xlab("Human Transcripts (10^3)")+
  ylab("Mouse Transcripts (10^3)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_equal(2)+
  geom_abline(slope = 1/9,lty=2)+geom_abline(slope = 9,lty=2)
dev.off()

mix.obj$ratio_human <- Barnyard.data[colnames(mix.obj),]$ratio_human

sum(Barnyard.data[human.cell,]$human_transcripts) /
  sum(Barnyard.data[human.cell, c("human_transcripts","mouse_transcripts")]) # 0.9918084
sum(Barnyard.data[mouse.cell,]$mouse_transcripts) /
  sum(Barnyard.data[mouse.cell, c("human_transcripts","mouse_transcripts")]) # 0.9725694

# 6. Correlation Among Cells ----------------------------------------------
library(pheatmap)

anno <- data.frame(row.names = rownames(sample.info),
                   "cell_line"=sample.info$Cell_Line,
                   "temp" = 1)
anno <- anno[order(anno$cell_line),]
anno <- anno[rownames(anno) %in% colnames(mix.obj),]
anno$temp <- NULL

rc_gene_var <- rc_gene[,
                       rownames(anno)]

identical(rownames(anno),colnames(rc_gene_var))

cor.gene <- matrix(nrow = nrow(anno),
                   ncol = nrow(anno))
rownames(cor.gene) <- rownames(anno)
colnames(cor.gene) <- rownames(anno)

for(i in 1:nrow(anno)){
  for (j in 1:nrow(anno)){
    cor.gene[i,j]=cor(rc_gene_var[,i],rc_gene_var[,j])
  }
}

pdf("Heatmap.cor.all.pdf",7,6)
pheatmap(cor.gene, main = "Correlation at Gene Level",
         cluster_rows = F,cluster_cols = F,
         annotation_row = anno, annotation_col = anno,
         annotation_names_row = F, annotation_names_col = F,
         annotation_colors = list(cell_line = mycolor),
         show_rownames = F,show_colnames = F)
dev.off()

pdf("Heatmap.cor.human.pdf",7,6)
pheatmap(cor.gene[rownames(cor.gene) %in% human.cell,
                  colnames(cor.gene) %in% human.cell], main = "Correlation at Gene Level",
         cluster_rows = F,cluster_cols = F,
         annotation_row = anno, annotation_col = anno,
         annotation_names_row = F, annotation_names_col = F,
         annotation_colors = list(cell_line = mycolor),
         show_rownames = F,show_colnames = F)
dev.off()
pdf("Heatmap.cor.mouse.pdf",7,6)
pheatmap(cor.gene[rownames(cor.gene) %in% mouse.cell,
                  colnames(cor.gene) %in% mouse.cell], main = "Correlation at Gene Level",
         cluster_rows = F,cluster_cols = F,
         annotation_row = anno, annotation_col = anno,
         annotation_names_row = F, annotation_names_col = F,
         annotation_colors = list(cell_line = mycolor),
         show_rownames = F,show_colnames = F)
dev.off()

rm(rc_gene_hg38, rc_gene_mm10)
rm(rc_gene, rc_trans, rc_gene_var)

sample.info.qc <- sample.info[sample.info$Pass_QC==1,]

plot.data <- data.frame()
for(one in unique(sample.info.qc$Cell_Line)){
  temp <- data.frame("cell_line" = one,
                     "Cor" = as.numeric(cor.gene[sample.info.qc[sample.info.qc$Cell_Line==one,]$Rename,
                                             sample.info.qc[sample.info.qc$Cell_Line==one,]$Rename]))
  plot.data <- rbind(plot.data, temp)
}
table(plot.data$cell_line)
plot.data <- plot.data[plot.data$Cor < 1,]

aggregate(plot.data$Cor, list(plot.data$cell_line), median)

pdf("Density.Cor.pdf",4,3)
ggplot(plot.data)+
  geom_density(aes(x = Cor, color = cell_line),)+
  scale_color_manual(values = mycolor)+
  theme_lzy()
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
                 y=log10(1+rc_ERCC$`9CL_Bc53_Bc01`),color="Bc75_Bc01"))+
  geom_point(aes(x=log10(1+ERCC.info$concentration_mix_1),
                 y=log10(1+rc_ERCC$`9CL_Bc53_Bc02`),color="Bc75_Bc02"))+
  geom_point(aes(x=log10(1+ERCC.info$concentration_mix_1),
                 y=log10(1+rc_ERCC$`9CL_Bc82_Bc31`),color="Bc75_Bc03"))+
  geom_point(aes(x=log10(1+ERCC.info$concentration_mix_1),
                 y=log10(1+rc_ERCC$`9CL_Bc82_Bc28`),color="Bc75_Bc04"))+
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
  scale_color_manual(values = mycolor)+
  theme_lzy()
dev.off()

# write sample info
write.table(sample.info,"sample.info.txt",quote = F,sep = "\t")
saveRDS(mix.obj,"obj.9CL.rds")

for (one in unique(sample.info.qc$Cell_Line)) {
  write.table(sample.info.qc[sample.info.qc$Cell_Line==one,]$Name,
              paste0("sample.list.9CL.",one),quote = F,sep = "\t",row.names = F,col.names = F)
}
