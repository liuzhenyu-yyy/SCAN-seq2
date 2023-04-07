setwd("E:/LabWork/Project/SCANSeq2/batch/2 9CL_Mix/")

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
sample.info <-
  read.table("data/summary.txt", header = T, row.names = 1)
sample.info$Mapped_percent <-
  sample.info$Mapped_Reads / sample.info$Raw_Reads
sample.info$QC_Reads <- NULL

sample.info$Library <- "9CL_Mix"
sample.info$Name <- rownames(sample.info)
sample.info$Rename <- paste("9CL_Mix", sample.info$Name, sep = "_")
sample.info$Barcode.5 <- substr(sample.info$Name, 3, 4)
sample.info$Barcode.3 <- substr(sample.info$Name, 8, 9)
rownames(sample.info) <- sample.info$Rename

sample.info$Cell_Line <- "none"
sample.info$Organism <- "none"

# merge umi_count matrix
rc_gene <- as.data.frame(fread("data/Gene.Count.txt", header = T))
rc_trans <- as.data.frame(fread("data/Trans.Count.txt", header = T))
rownames(rc_gene) <- rc_gene$Name
rownames(rc_trans) <- rc_trans$Name
rc_gene <- rc_gene[, -1]
rc_trans <- rc_trans[, -1]

rc_gene <- rc_gene[, sample.info$Name]
rc_trans <- rc_trans[, sample.info$Name]

colnames(rc_gene) <- sample.info$Rename
colnames(rc_trans) <- sample.info$Rename

gc()

# remove ERCC rows
ERCC <- grep("ERCC|RGC", rownames(rc_gene), value = T)
rc_ERCC <- rc_gene[ERCC, ]
rc_gene <- rc_gene[!rownames(rc_gene) %in% ERCC, ]
rc_trans <- rc_trans[!rownames(rc_trans) %in% ERCC, ]

human.gene <- grep("ENSG", rownames(rc_gene), value = T)
mouse.gene <- grep("ENSMU", rownames(rc_gene), value = T)
human.trans <- grep("ENST", rownames(rc_trans), value = T)
mouse.trans <- grep("ENSMU", rownames(rc_trans), value = T)

rc_gene_human <- rc_gene[human.gene, ]
rc_trans_human <- rc_trans[human.trans, ]
rc_gene_mouse <- rc_gene[mouse.gene, ]
rc_trans_mouse <- rc_trans[mouse.trans, ]
rm(rc_gene, rc_trans)

n1 <- colSums(rc_gene_human >= 1)
n2 <- colSums(rc_gene_mouse >= 1)
sample.info$Gene_Detected <- pmax(n1, n2)
n1 <- colSums(rc_trans_human >= 1)
n2 <- colSums(rc_trans_mouse >= 1)
sample.info$Trans_Detected <- pmax(n1, n2)
rm(n1, n2)
quantile(sample.info$Gene_Detected)
quantile(sample.info$Trans_Detected)

mycolor <- brewer.pal(11, "Set3")[c(1:8, 11)]
names(mycolor) <- c("GM12878",
                    "HepG2",
                    "Hela",
                    "H9",
                    "K562",
                    "293T",
                    "AtT20",
                    "MEF",
                    "3T3")
# mycolor["3T3"] <- "gray"
mycolor["HepG2"] <- "#ffed6f"

# 2. QC of Cell -----------------------------------------------------------
# basic statistics
pdf("Vln.QC.stat.pdf", 2, 2)
ggplot(sample.info) +
  geom_violin(
    aes(
      x = Cell_Line,
      y = log10(Raw_Reads),
      fill = Cell_Line
    ),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(
      x = Cell_Line,
      y = log10(Raw_Reads),
      color = Cell_Line
    ),
    alpha = 0,
    show.legend = F
  ) +
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  theme_bw()
ggplot(sample.info) +
  geom_violin(
    aes(
      x = Cell_Line,
      y = log10(Mapped_Reads),
      fill = Cell_Line
    ),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(
      x = Cell_Line,
      y = log10(Mapped_Reads),
      color = Cell_Line
    ),
    alpha = 0,
    show.legend = F
  ) +
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  theme_bw()
ggplot(sample.info) +
  geom_violin(
    aes(x = Cell_Line, y = Mapped_percent, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Mapped_percent, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  theme_bw()
ggplot(sample.info) +
  geom_violin(
    aes(x = Cell_Line, y = Gene_Detected, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Gene_Detected, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  theme_bw()
ggplot(sample.info) +
  geom_violin(
    aes(x = Cell_Line, y = Trans_Detected, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Trans_Detected, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  theme_bw()
dev.off()

pdf("Dot.Gene_Trans.pdf", 4, 3.5)
ggplot(sample.info) +
  geom_point(aes(x = Mapped_Reads, y = Gene_Detected)) +
  scale_color_manual(values = mycolor) +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

quantile(sample.info$Mapped_Reads) # 26064.50

# rename genes and trans
#human gene
human.gene <- substr(human.gene, 1, 15)
rownames(rc_gene_human) <- human.gene

rc_gene_human <-
  rc_gene_human[rownames(rc_gene_human) %in% gene.info.hg38$`Gene stable ID`, ]
rc_gene_human <-
  rc_gene_human[!duplicated(gene.info.hg38[rownames(rc_gene_human), ]$`Gene name`), ]
rownames(rc_gene_human) <-
  gene.info.hg38[rownames(rc_gene_human), ]$`Gene name`

#human trans
temp <- substr(rownames(rc_trans_human), 1, 15)
temp <- duplicated(temp)
rc_trans_human <- rc_trans_human[!temp, ]
rownames(rc_trans_human) <- substr(rownames(rc_trans_human), 1, 15)
table(rownames(rc_trans_human) %in% gene2trans$trans.id)

rc_trans_human <-
  rc_trans_human[rownames(rc_trans_human) %in% gene2trans$trans.id, ]
rownames(rc_trans_human) <-
  gene2trans[rownames(rc_trans_human), ]$trans.symbol

#mouse genes
mouse.gene <- substr(mouse.gene, 1, 18)
rownames(rc_gene_mouse) <- mouse.gene

rc_gene_mouse <-
  rc_gene_mouse[rownames(rc_gene_mouse) %in% gene.info.mm10$`Gene stable ID`, ]
rc_gene_mouse <-
  rc_gene_mouse[!duplicated(gene.info.mm10[rownames(rc_gene_mouse), ]$`Gene name`), ]
rownames(rc_gene_mouse) <-
  gene.info.mm10[rownames(rc_gene_mouse), ]$`Gene name`

#merge
identical(colnames(rc_gene_mouse), colnames(rc_gene_human))
identical(colnames(rc_trans_mouse), colnames(rc_trans_human))
rc_gene <- rbind(rc_gene_human, rc_gene_mouse)
rc_trans <- rbind(rc_trans_human, rc_trans_mouse)
rm(rc_gene_mouse, rc_gene_human,
   rc_trans_human, rc_trans_mouse)
saveRDS(rc_gene, "data/rc_gene.RDS")
saveRDS(rc_trans, "data/rc_trans.RDS")
saveRDS(rc_ERCC, "data/rc_ERCC.RDS")
gc()

# Create Object
mix.obj <- CreateSeuratObject(rc_gene,
                              meta.data = sample.info,
                              project = "9CL_Mix",
                              assay = "Gene")
mix.obj

mix.obj$percent.mt <-
  PercentageFeatureSet(mix.obj, pattern = "^MT-|^mt-")

pdf("Scatter.QC.pdf", 5, 4)
FeatureScatter(mix.obj, feature1 = "Gene_Detected", feature2 = "percent.mt") +
  geom_hline(yintercept = 15,
             lty = 2,
             color = "red") +
  geom_vline(xintercept = 2000,
             lty = 2,
             color = "red")
dev.off()

table(mix.obj$Gene_Detected > 2000 & mix.obj$percent.mt < 15)

mix.obj <-
  subset(mix.obj, subset = Gene_Detected > 2000 & percent.mt < 15)

QC.g2k <- table(mix.obj$Cell_Line)
table(sample.info$Cell_Line)
QC.g2k / table(sample.info$Cell_Line)

sample.info$Pass_QC <- 0
sample.info[sample.info$Rename %in% colnames(mix.obj), ]$Pass_QC <-
  1

pdf("Vln.gene_trans.QC.pdf", 2, 2)
ggplot(sample.info[sample.info$Pass_QC == 1, ]) +
  geom_violin(aes(
    x = Cell_Line,
    y = (Gene_Detected),
    fill = Cell_Line
  ),
  alpha = 0.7,
  show.legend = F) +
  geom_boxplot(aes(
    x = Cell_Line,
    y = (Gene_Detected),
    color = Cell_Line
  ),
  alpha = 0,
  show.legend = F) +
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  theme_bw()
ggplot(sample.info[sample.info$Pass_QC == 1, ]) +
  geom_violin(aes(
    x = Cell_Line,
    y = (Trans_Detected),
    fill = Cell_Line
  ),
  alpha = 0.7,
  show.legend = F) +
  geom_boxplot(aes(
    x = Cell_Line,
    y = (Trans_Detected),
    color = Cell_Line
  ),
  alpha = 0,
  show.legend = F) +
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) +
  theme_bw()
dev.off()

# 3. gene-level analysis -----------------------------------------------------
mix.obj <-
  NormalizeData(mix.obj,
                normalization.method = "LogNormalize",
                scale.factor = 10000)

mix.obj <-
  FindVariableFeatures(mix.obj, selection.method = "vst", nfeatures = 3000)

#VariableFeaturePlot(mix.obj)

mix.obj <- ScaleData(mix.obj)

mix.obj <- RunPCA(mix.obj,
                  npcs = 30,
                  features = VariableFeatures(object = mix.obj))

#DimHeatmap(mix.obj, dims = 1:10, balanced = TRUE)
ElbowPlot(mix.obj) #10~15

pc.use = 1:8

# UMAP
mix.obj <- RunUMAP(mix.obj,
                   dims = pc.use,
                   n.neighbors = 20,
                   min.dist = 0.1)
mix.obj <-
  FindNeighbors(mix.obj,
                reduction = "pca",
                dims = pc.use,
                k.param = 20)
mix.obj <- FindClusters(mix.obj, resolution = 0.3)

pdf("UMAP.Cluster.pdf", 4.4, 3.5)
DimPlot(mix.obj,
        reduction = "umap",
        label = T,
        pt.size = 0.6) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_color_manual(values = brewer.pal(10, "Set3")) +
  coord_equal(0.7)
DimPlot(mix.obj,
        reduction = "umap",
        label = F,
        pt.size = 0.6) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_color_manual(values = brewer.pal(10, "Set3")) +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) +
  coord_equal(0.7)
dev.off()

marker.list <- c("IGHM",
                 "ALB",
                 "LY6K",
                 # GM12878, HepG2, Hela
                 "CD3D",
                 "HBA1",
                 "BEX3",
                 # H9, K562, 293T
                 "Apobec3",
                 "Col1a1",
                 "Lgals3" # mES, 3T3)

                 pdf("UMAP.marker.pdf", 4.5, 4)
                 for (one in marker.list) {
                   p <- FeaturePlot(mix.obj, one, order = F) +
                     coord_equal(0.8)
                   plot(p)
                 }
                 dev.off()
                 FeaturePlot(mix.obj, "UCHL1", order = F) +
                   coord_equal(0.8)

                 pdf("Vln.marker.pdf", 8, 6)
                 VlnPlot(mix.obj, features = marker.list, cols = brewer.pal(10, "Set3"))
                 dev.off()

                 # cell type
                 p <- DimPlot(mix.obj,
                              reduction = "umap",
                              label = T,
                              pt.size = 0.6) +
                   theme_bw() + theme(panel.grid = element_blank())
                 GM12878.cell <- CellSelector(p)
                 HepG2.cell <- CellSelector(p)
                 Hela.cell <- CellSelector(p)
                 H9.cell <- CellSelector(p)
                 K562.cell <- CellSelector(p)
                 cell.293T <- CellSelector(p)
                 AtT20.cell <- CellSelector(p)
                 MEF.cell <- CellSelector(p)
                 cell.3T3 <- CellSelector(p)

                 mix.obj@meta.data[GM12878.cell, ]$Cell_Line <- "GM12878"
                 mix.obj@meta.data[HepG2.cell, ]$Cell_Line <- "HepG2"
                 mix.obj@meta.data[Hela.cell , ]$Cell_Line <- "Hela"
                 mix.obj@meta.data[H9.cell, ]$Cell_Line <- "H9"
                 mix.obj@meta.data[K562.cell, ]$Cell_Line <- "K562"
                 mix.obj@meta.data[cell.293T , ]$Cell_Line <- "293T"
                 mix.obj@meta.data[AtT20.cell, ]$Cell_Line <- "AtT20"
                 mix.obj@meta.data[MEF.cell, ]$Cell_Line <- "MEF"
                 mix.obj@meta.data[cell.3T3, ]$Cell_Line <- "3T3"
                 table(mix.obj$Cell_Line)

                 Idents(mix.obj) <- mix.obj$Cell_Line
                 mix.obj@meta.data$Cell_Line <- factor(
                   mix.obj@meta.data$Cell_Line,
                   levels = c(
                     "GM12878",
                     "HepG2",
                     "Hela",
                     "H9",
                     "K562",
                     "293T",
                     "AtT20",
                     "MEF",
                     "3T3"
                   )
                 )
                 pdf("UMAP.CellType.pdf", 4.8, 3.5)
                 DimPlot(
                   mix.obj,
                   reduction = "umap",
                   group.by = "Cell_Line",
                   label = T,
                   cols = mycolor,
                   pt.size = 0.3
                 ) +
                   theme_bw() + theme(panel.grid = element_blank()) +
                   coord_equal(0.8)
                 DimPlot(
                   mix.obj,
                   reduction = "umap",
                   group.by = "Cell_Line",
                   label = F,
                   cols = mycolor,
                   pt.size = 0.3
                 ) +
                   theme_bw() + theme(panel.grid = element_blank()) +
                   theme(
                     axis.ticks = element_blank(),
                     axis.text = element_blank(),
                     axis.title = element_blank()
                   ) +
                   coord_equal(0.8) + NoLegend()
                 dev.off()

                 pdf("UMAP.Quality.pdf", 4.5, 3.5)
                 FeaturePlot(
                   mix.obj,
                   "nFeature_Gene",
                   reduction = "umap",
                   pt.size = 0.6,
                   cols = c("gray", "red")
                 ) +
                   theme_bw() + theme(panel.grid = element_blank()) +
                   coord_equal(0.8)
                 FeaturePlot(
                   mix.obj,
                   "nCount_Gene",
                   reduction = "umap",
                   pt.size = 0.6,
                   cols = c("gray", "red")
                 ) +
                   theme_bw() + theme(panel.grid = element_blank()) +
                   coord_equal(0.8)
                 FeaturePlot(
                   mix.obj,
                   "Mapped_Reads",
                   reduction = "umap",
                   pt.size = 0.6,
                   cols = c("gray", "red")
                 ) +
                   theme_bw() + theme(panel.grid = element_blank()) +
                   coord_equal(0.8)
                 FeaturePlot(
                   mix.obj,
                   "Mapped_percent",
                   reduction = "umap",
                   pt.size = 0.6,
                   cols = c("gray", "red")
                 ) +
                   theme_bw() + theme(panel.grid = element_blank()) +
                   coord_equal(0.8)
                 FeaturePlot(
                   mix.obj,
                   "percent.mt",
                   reduction = "umap",
                   pt.size = 0.6,
                   cols = c("gray", "red")
                 ) +
                   theme_bw() + theme(panel.grid = element_blank()) +
                   coord_equal(0.8)
                 dev.off()


                 # marker heatmap
                 Idents(mix.obj) <- mix.obj$Cell_Line
                 marker.all <- FindAllMarkers(
                   mix.obj,
                   only.pos = TRUE,
                   min.pct = 0.25,
                   logfc.threshold = 0.25
                 )
                 table(marker.all$cluster)

                 top10 <-
                   marker.all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                 table(top10$cluster)

                 pdf("Heatmap.Marker.all.pdf", 6, 5)
                 mix.obj <- ScaleData(mix.obj, features =  top10$gene)
                 DoHeatmap(
                   mix.obj,
                   features = top10$gene,
                   draw.lines = F,
                   group.colors = mycolor
                 ) +
                   scale_fill_gradientn(colors = rev(brewer.pal(7, "RdYlBu")))
                 mix.obj <- ScaleData(mix.obj)
                 dev.off()

                 # 5. Barnyard plot --------------------------------------------------------
                 sample.info[colnames(mix.obj), ]$Cell_Line <-
                   as.character(mix.obj$Cell_Line)
                 table(sample.info$Cell_Line)
                 sample.info$Organism <-
                   ifelse(
                     sample.info$Cell_Line %in% c("GM12878", "HepG2", "Hela", "H9", "K562", "293T"),
                     "Human",
                     "Mouse"
                   )
                 sample.info[sample.info$Pass_QC == 0, ]$Organism <- "none"

                 mouse.trans <-
                   rownames(rc_trans)[grep("ENSMUST", rownames(rc_trans))]
                 human.trans <-
                   rownames(rc_trans)[-grep("ENSMUST", rownames(rc_trans))]
                 nrow(rc_trans) - length(mouse.trans) - length(human.trans)

                 Barnyard.data <- data.frame(
                   row.names = colnames(rc_trans),
                   "cell" = colnames(rc_trans),
                   "human_transcripts" = colSums(rc_trans[human.trans, ]),
                   "mouse_transcripts" = colSums(rc_trans[mouse.trans, ])
                 )

                 Barnyard.data <- Barnyard.data[rownames(sample.info), ]

                 rm(mouse.trans, human.trans)

                 Barnyard.data$Cell_Line <- sample.info$Cell_Line

                 Barnyard.data$ratio_human <-
                   Barnyard.data$human_transcripts / (Barnyard.data$human_transcripts +
                                                        Barnyard.data$mouse_transcripts)
                 Barnyard.data$ratio_mouse <-
                   Barnyard.data$mouse_transcripts / (Barnyard.data$human_transcripts +
                                                        Barnyard.data$mouse_transcripts)

                 Barnyard.data$ratio_human %>% hist()

                 Barnyard.data$Organism <- "Mixed"
                 Barnyard.data[Barnyard.data$ratio_human > 0.9, ]$Organism <- "Human"
                 Barnyard.data[Barnyard.data$ratio_mouse > 0.9, ]$Organism <- "Mouse"
                 table(Barnyard.data$Organism)
                 table(Barnyard.data$Organism) / nrow(Barnyard.data)

                 Barnyard.data.sub <- Barnyard.data[colnames(mix.obj), ]
                 table(Barnyard.data.sub$Organism)
                 table(Barnyard.data.sub$Organism, Barnyard.data.sub$Cell_Line)
                 table(Barnyard.data.sub$Organism) / nrow(Barnyard.data.sub)

                 pdf("Hist_Ratio_Human.cDNA.pdf", 8, 6)
                 ggplot(Barnyard.data.sub) +
                   geom_histogram(aes(x = ratio_human, fill = Cell_Line), binwidth = 0.01) +
                   geom_vline(xintercept = 0.1, lty = 2) +
                   geom_vline(xintercept = 0.9, lty = 2) +
                   facet_wrap( ~ Cell_Line) +
                   theme_bw() + xlab("Percent of Human UMIs") +
                   scale_fill_manual(values = mycolor)
                 dev.off()


                 pdf("Barnyard.cDNA.pdf", 4, 3)
                 ggplot(Barnyard.data.sub) +
                   geom_point(aes(
                     x = human_transcripts / 1000,
                     y = mouse_transcripts / 1000,
                     color = Organism
                   ),
                   size = 1) +
                   theme_bw() +
                   xlab("Human UMI counts (k)") +
                   ylab("Mouse UMI counts (k)") +
                   theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()) +
                   coord_equal(2) +
                   scale_color_manual(values = c("#d8756e", "#808183", "#6dbbc1")) +
                   geom_abline(slope = 1 / 9, lty = 2) + geom_abline(slope = 9, lty = 2)
                 ggplot(Barnyard.data) +
                   geom_point(aes(
                     x = human_transcripts / 1000,
                     y = mouse_transcripts / 1000,
                     color = Cell_Line
                   ),
                   size = 0.3) +
                   theme_bw() +
                   scale_color_manual(values = mycolor) +
                   xlab("Human Transcripts (10^3)") +
                   ylab("Mouse Transcripts (10^3)") +
                   theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()) +
                   coord_equal(2) +
                   geom_abline(slope = 1 / 9, lty = 2) + geom_abline(slope = 9, lty = 2)
                 dev.off()

                 mix.obj$ratio_human <- Barnyard.data[colnames(mix.obj), ]$ratio_human

                 # 7. ERCC ---------------------------------------------------------------
                 # read ERCC info
                 ERCC.info <-
                   as.data.frame(fread("E:/LabWork/genome/Thermo_ERCC.txt"))
                 ERCC.info <- ERCC.info[, c(1:5)]
                 colnames(ERCC.info) <- c("order",
                                          "id",
                                          "group",
                                          "concentration_mix_1",
                                          "concentration_mix_2")
                 rc_ERCC <- rc_ERCC[ERCC.info$id, ]

                 ggplot() +
                   geom_point(aes(
                     x = log10(1 + ERCC.info$concentration_mix_1),
                     y = log10(1 + rc_ERCC$`9CL_Mix_Bc33_Bc01`),
                     color = "Cell1"
                   )) +
                   geom_point(aes(
                     x = log10(1 + ERCC.info$concentration_mix_1),
                     y = log10(1 + rc_ERCC$`9CL_Mix_Bc33_Bc09`),
                     color = "Cell2"
                   )) +
                   geom_point(aes(
                     x = log10(1 + ERCC.info$concentration_mix_1),
                     y = log10(1 + rc_ERCC$`9CL_Mix_Bc33_Bc32`),
                     color = "Cell3"
                   )) +
                   geom_point(aes(
                     x = log10(1 + ERCC.info$concentration_mix_1),
                     y = log10(1 + rc_ERCC$`9CL_Mix_Bc34_Bc01`),
                     color = "Cell4"
                   )) +
                   xlab("log10(ERCC concentration (attomoles/ul))") +
                   ylab("log(UMI Count + 1)") +
                   coord_fixed() +
                   theme_lzy()

                 sample.info$cor_ERCC_rc <- 0
                 for (cell in rownames(sample.info)) {
                   sample.info[cell, ]$cor_ERCC_rc <- cor(rc_ERCC[, cell],
                                                          ERCC.info$concentration_mix_1)
                 }

                 pdf("Density.cor.ERCC.pdf", 4.5, 3)
                 ggplot(sample.info[sample.info$Pass_QC == 1, ]) +
                   geom_density(aes(x = cor_ERCC_rc, color = Cell_Line), alpha = 0.5) +
                   scale_color_manual(values = mycolor) +
                   theme_lzy()
                 dev.off()

                 # write sample info
                 write.table(sample.info,
                             "sample.info.txt",
                             quote = F,
                             sep = "\t")
                 saveRDS(mix.obj, "obj.9CL_Mix.rds")

                 table(sample.info$Cell_Line, sample.info$Organism)
                 table(sample.info$Organism)
                 write.table(
                   sample.info[sample.info$Organism == "Human", ]$Name,
                   "sample.list.hg38",
                   quote = F,
                   sep = "\t",
                   row.names = F,
                   col.names = F
                 )
                 write.table(
                   sample.info[sample.info$Organism == "Mouse", ]$Name,
                   "sample.list.mm10",
                   quote = F,
                   sep = "\t",
                   row.names = F,
                   col.names = F
                 )
                 sample.info.qc <- sample.info[sample.info$Pass_QC == 1, ]
                 for (one in unique(sample.info.qc$Cell_Line)) {
                   write.table(
                     sample.info.qc[sample.info.qc$Cell_Line == one, ]$Name,
                     paste0("sample.list.9Mix.", one),
                     quote = F,
                     sep = "\t",
                     row.names = F,
                     col.names = F
                   )
                 }
