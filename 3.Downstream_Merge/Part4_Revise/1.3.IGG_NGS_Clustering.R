setwd("E:/LabWork/Project/SCANSeq2/revise/1.3")

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

source("E:/LabWork/code/MyFunction.R")
ReadGeneInfo("hg38", classify.protein = T)

# color for cell lines
color.cell <- brewer.pal(11, "Set3")[c(1:8, 10)]
names(color.cell) <- c("GM12878",
                       "HepG2",
                       "Hela",
                       "H9",
                       "K562",
                       "293T",
                       "AtT20",
                       "MEF",
                       "3T3")
color.cell["HepG2"] <- "#ffed6f"
color.cell <- color.cell[c("Hela", "HepG2")]
color.time <- brewer.pal(6, "Set2")[1:3]
color.concentration <- brewer.pal(6, "Set2")[4:6]
color.Phase <- brewer.pal(3, "Set1")
color.cluster <- brewer.pal(9, "Set3")[c(1, 3:8)]

# 1. load NGS and TGS data ------------------------------------------------

## 1.1. SCAN-seq2 TGS data ----
sample.info.TGS <-
  read.table("data/TGS/summary.txt",
             header = T,
             row.names = 1)

sample.info.TGS$Library <- "SCAN2-IGG"
sample.info.TGS$Name <- rownames(sample.info.TGS)
sample.info.TGS$Rename <-
  paste("SCAN2-IGG", sample.info.TGS$Name, sep = "_")
sample.info.TGS$Barcode.5 <- substr(sample.info.TGS$Name, 3, 4)
sample.info.TGS$Barcode.3 <- substr(sample.info.TGS$Name, 8, 9)
sample.info.TGS$Method <- "SCAN-seq2"
rownames(sample.info.TGS) <- sample.info.TGS$Rename

sample.info.TGS$Cell_Line <- "Hela"

sample.info.TGS$Concentration <- 0
sample.info.TGS[sample.info.TGS$Barcode.5 %in%
                  as.character(c(54:56, 63:65, 72:74)), ]$Concentration <-
  10
sample.info.TGS[sample.info.TGS$Barcode.5 %in%
                  as.character(c(57:59, 66:68, 75)), ]$Concentration <-
  30
table(sample.info.TGS$Concentration)

sample.info.TGS$Time <- 0
sample.info.TGS[sample.info.TGS$Barcode.5 %in%
                  as.character(51:59), ]$Time <- 6
sample.info.TGS[sample.info.TGS$Barcode.5 %in%
                  as.character(60:68), ]$Time <- 24
sample.info.TGS[sample.info.TGS$Barcode.5 %in%
                  as.character(69:75), ]$Time <- 48
table(sample.info.TGS$Time)
table(sample.info.TGS$Time, sample.info.TGS$Concentration)

sample.info.TGS$nGenes <- colSums(rc.TGS >= 1)

# rc
rc.TGS <- as.data.frame(fread("data/TGS/Gene.Count.txt", header = T))
rownames(rc.TGS) <- rc.TGS$Name
rc.TGS <- rc.TGS[, -1]
rc.TGS <- rc.TGS[, sample.info.TGS$Name]
colnames(rc.TGS) <- sample.info.TGS$Rename

# remove ERCC rows
ERCC <- grep("ERCC|RGC", rownames(rc.TGS), value = T)
rc.TGS <- rc.TGS[!rownames(rc.TGS) %in% ERCC, ]

human.gene <- grep("ENSG", rownames(rc.TGS), value = T)
rc.TGS <- rc.TGS[human.gene, ]

# rename genes
human.gene <- substr(human.gene, 1, 15)
rownames(rc.TGS) <- human.gene
rc.TGS <-
  rc.TGS[rownames(rc.TGS) %in% gene.info.hg38$`Gene stable ID`, ]
rownames(gene.info.hg38) <- gene.info.hg38$`Gene stable ID`
rc.TGS <-
  rc.TGS[!duplicated(gene.info.hg38[rownames(rc.TGS), ]$`Gene name`), ]
rownames(rc.TGS) <- gene.info.hg38[rownames(rc.TGS), ]$`Gene name`

## 1.2. STRT NGS data ----
samples <- dir("data/NGS/") %>% gsub("\\..+?$", "", .) %>% unique
sample.info.NGS <- data.frame()
rc.NGS <- data.frame()

for (one in samples) {
  message(paste0(Sys.time(), ". Loading data for library: ", one))
  sample.info <-
    read.table(paste0("data/NGS/", one, ".gene_tran_num.xls"), header = T)
  colnames(sample.info) <- c("Rename", "nGenes", "UMI_count")
  sample.info$Rename <- gsub("-", ".", sample.info$Rename)
  sample.info$Library <- one
  rownames(sample.info) <- sample.info$Rename
  sample.info.NGS <- rbind(sample.info.NGS, sample.info)

  rc <-
    read.table(
      paste0("data/NGS/", one, ".umi_counts.xls"),
      header = T,
      row.names = 1
    )
  rc <- as.data.frame(t(rc))
  rc.NGS <- rbind(rc.NGS, rc)
}
rm(sample.info, rc)
rc.NGS <- as.data.frame(t(rc.NGS))

sample.info.NGS$Method <- "STRT-seq"
sample.info.NGS$Cell_Line <- "Hela"

sample.info.NGS$Concentration <- 0
sample.info.NGS[grep("2|5|8", sample.info.NGS$Library), ]$Concentration <-
  10
sample.info.NGS[grep("3|6|9", sample.info.NGS$Library), ]$Concentration <-
  30

sample.info.NGS$Time <- 0
sample.info.NGS[grep("1|2|3", sample.info.NGS$Library), ]$Time <- 6
sample.info.NGS[grep("4|5|6", sample.info.NGS$Library), ]$Time <- 24
sample.info.NGS[grep("7|8|9", sample.info.NGS$Library), ]$Time <- 48
table(sample.info.NGS$Concentration, sample.info.NGS$Time)


# 2. Primary Seurat pipeline ---------------------------------------------

## 2.1. Quality Control  ----
TGS.obj <- CreateSeuratObject(rc.TGS, meta.data = sample.info.TGS)
NGS.obj <- CreateSeuratObject(rc.NGS, meta.data = sample.info.NGS)

TGS.obj$percent.mt <-
  PercentageFeatureSet(TGS.obj, pattern = "^MT-|^mt-")
NGS.obj$percent.mt <-
  PercentageFeatureSet(NGS.obj, pattern = "^MT-|^mt-")

pdf("Scatter.QC.pdf", 10, 4)
FeatureScatter(TGS.obj,
               feature1 = "nGenes",
               feature2 = "percent.mt",
               group.by = "Concentration") + ggtitle("SCAN-seq2") +
  FeatureScatter(NGS.obj,
                 feature1 = "nGenes",
                 feature2 = "percent.mt",
                 group.by = "Concentration") + ggtitle("STRT-seq")
FeatureScatter(TGS.obj,
               feature1 = "nGenes",
               feature2 = "percent.mt",
               group.by = "Time") + ggtitle("SCAN-seq2") +
  FeatureScatter(NGS.obj,
                 feature1 = "nGenes",
                 feature2 = "percent.mt",
                 group.by = "Time") + ggtitle("STRT-seq")
dev.off()

TGS.obj <-
  subset(TGS.obj, subset = nGenes > 2500 &
           percent.mt < 15) # 579 / 793
NGS.obj <-
  subset(NGS.obj, subset = nGenes > 2500 &
           percent.mt < 15) # 581 / 816

TGS.obj <- TGS.obj %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(nfeatures = 3000, ) %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors()
ElbowPlot(TGS.obj)

NGS.obj <- NGS.obj %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(nfeatures = 3000, ) %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors()
ElbowPlot(NGS.obj)

TGS.obj <- RunUMAP(TGS.obj,
                   dims = 1:5,
                   n.neighbors = 50,
                   min.dist = 0.1)
NGS.obj <- RunUMAP(NGS.obj,
                   dims = 1:5,
                   n.neighbors = 50,
                   min.dist = 0.1)

TGS.obj <- FindClusters(TGS.obj, resolution = 0.3)
NGS.obj <- FindClusters(NGS.obj, resolution = 0.3)

## 2.2. Cell cycle regression
TGS.obj <- CellCycleScoring(
  TGS.obj,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = T
)
NGS.obj <- CellCycleScoring(
  NGS.obj,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = T
)

pdf("UMAP.befor.regression.pdf", 8, 6)
DimPlot(TGS.obj, group.by = "seurat_clusters", cols = color.cluster) +
  DimPlot(TGS.obj, group.by = "Concentration", cols = color.concentration) +
  DimPlot(TGS.obj, group.by = "Time", cols = color.time) +
  DimPlot(TGS.obj, group.by = "Phase", cols = color.Phase)

DimPlot(NGS.obj, group.by = "seurat_clusters", cols = color.cluster) +
  DimPlot(NGS.obj, group.by = "Concentration", cols = color.concentration) +
  DimPlot(NGS.obj, group.by = "Time", cols = color.time) +
  DimPlot(NGS.obj, group.by = "Phase", cols = color.Phase)
dev.off()


TGS.obj <- ScaleData(
  TGS.obj,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(TGS.obj)
) %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:5,
          n.neighbors = 50,
          min.dist = 0.1)
NGS.obj <- ScaleData(
  NGS.obj,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(NGS.obj)
) %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:5,
          n.neighbors = 50,
          min.dist = 0.1)

TGS.obj <- FindClusters(TGS.obj, resolution = 0.3)
NGS.obj <- FindClusters(NGS.obj, resolution = 0.3)

pdf("UMAP.after.regression.pdf", 8, 6)
DimPlot(TGS.obj, group.by = "seurat_clusters", cols = color.cluster) +
  DimPlot(TGS.obj, group.by = "Concentration", cols = color.concentration) +
  DimPlot(TGS.obj, group.by = "Time", cols = color.time) +
  DimPlot(TGS.obj, group.by = "Phase", cols = color.Phase)

DimPlot(NGS.obj, group.by = "seurat_clusters", cols = color.cluster) +
  DimPlot(NGS.obj, group.by = "Concentration", cols = color.concentration) +
  DimPlot(NGS.obj, group.by = "Time", cols = color.time) +
  DimPlot(NGS.obj, group.by = "Phase", cols = color.Phase)
dev.off()
