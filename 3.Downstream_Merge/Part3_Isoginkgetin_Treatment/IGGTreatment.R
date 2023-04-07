setwd("E:/LabWork/Project/SCANSeq2/merge/Part3/")

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DTUrtle)
library(Vennerable)

source("E:/LabWork/code/MyFunction.R")
ReadGeneInfo("hg38", classify.protein = T)
ReadGeneInfo("mm10")
gtf <- import_gtf(
  "E:/LabWork/genome/Homo_sapiens.GRCh38.101.gtf",
  feature_type = NULL,
  out_df = FALSE
)

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

# 1. Dimensional reduction and clustering ---------------------------------
dir.create("DR_Clustering")

obj.IGG <- readRDS("../../batch/6 IGG/obj.IGG.rds")
rc_trans <-
  readRDS("E:/LabWork/Project/SCANSeq2/batch/6 IGG/data/rc_trans.RDS")

sample.info <- read.table("../../batch/6 IGG/sample.info.txt")
sample.info <- sample.info[sample.info$Pass_QC == 1,]

pdf("DR_Clustering/UMAP_Cell_Line.pdf", 5, 4)
DimPlot(obj.IGG, group.by = "Cell_Line", cols = color.cell) +
  coord_equal(1.3)
dev.off()

pdf("DR_Clustering/UMAP_Treatment.pdf", 5, 4)
DimPlot(obj.IGG, group.by = "Concentration", cols = color.concentration) +
  coord_equal(1.3)
DimPlot(obj.IGG, group.by = "Time", cols = color.time) +
  coord_equal(1.3)
dev.off()

obj.IGG <- FindClusters(obj.IGG, resolution = 0.3)
pdf("DR_Clustering/UMAP_Cluster.pdf", 5, 4)
DimPlot(obj.IGG, group.by = "Gene_snn_res.0.2",
        cols = color.cluster) +
  coord_equal(1.3)
dev.off()

# add cell cycle score
obj.IGG <- CellCycleScoring(
  obj.IGG,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = T
)
table(obj.IGG$Phase)

pdf("DR_Clustering/UMAP_Cell_Cycle_Phase.pdf", 5, 4)
DimPlot(obj.IGG, group.by = "Phase", cols = color.Phase) +
  coord_equal(1.3)
dev.off()

pdf("temp.pdf", 6.6, 3)
multiplot(p1, p2, cols = 2)
dev.off()



# 2. subset cell lines ----------------------------------------------------
obj.Hela <- subset(obj.IGG, subset = Cell_Line == "Hela")
obj.HepG2 <- subset(obj.IGG, subset = Cell_Line == "HepG2")

isoform_IGG6 <-
  read.table("../../batch/6 IGG/data/IGG6/Isoform.category.tsv", header = T)
isoform_IGG24 <-
  read.table("../../batch/6 IGG/data/IGG24/Isoform.category.tsv", header = T)
isoform_IGG48 <-
  read.table("../../batch/6 IGG/data/IGG48/Isoform.category.tsv", header = T)
isoform_IGG6$Library <- "IGG6"
isoform_IGG24$Library <- "IGG24"
isoform_IGG48$Library <- "IGG48"
isoform <- rbind(isoform_IGG6, isoform_IGG24, isoform_IGG48)
isoform$Rename <- paste(isoform$Library, isoform$cell, sep = "_")
rownames(isoform) <- isoform$Rename

isoform$percent.FSM <- isoform$FSM / isoform$Total_Trans
isoform$percent.ISM <- isoform$ISM / isoform$Total_Trans
isoform$percent.NIC <- isoform$NIC / isoform$Total_Trans
isoform$percent.NNC <- isoform$NNC / isoform$Total_Trans
isoform$percent.CJ <- isoform$NIC_CJ / isoform$Total_Trans
isoform$percent.CS <- isoform$NIC_CS / isoform$Total_Trans
isoform$percent.IR <- isoform$NIC_IR / isoform$Total_Trans
isoform$percent.NIC.CJ <- isoform$NIC_CJ / isoform$NIC
isoform$percent.NIC.CS <- isoform$NIC_CS / isoform$NIC
isoform$percent.NIC.IR <- isoform$NIC_IR / isoform$NIC

rm(isoform_IGG6, isoform_IGG24, isoform_IGG48)

## 2.1. Hela ---------------------------------------
dir.create("Hela")
obj.Hela <- obj.Hela %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>% FindNeighbors()
ElbowPlot(obj.Hela)

obj.Hela <- RunUMAP(
  obj.Hela,
  dims = 1:5,
  n.neighbors = 30,
  min.dist = 0.1
)
obj.Hela <- FindClusters(obj.Hela, resolution = 0.2)

pdf("Hela/UMAP.befor.regression.pdf", 8, 6)
multiplot(
  DimPlot(obj.Hela, group.by = "Gene_snn_res.0.2",
          cols = color.cluster) + coord_equal(1.5),
  DimPlot(obj.Hela, group.by = "Phase", cols = color.Phase) + coord_equal(1.5),
  DimPlot(obj.Hela, group.by = "Concentration",
          cols = color.concentration) + coord_equal(1.5),
  DimPlot(obj.Hela, group.by = "Time", cols = color.time) + coord_equal(1.5),
  cols = 2
)
dev.off()

# regress cell cycle
obj.Hela <- ScaleData(
  obj.Hela,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(obj.Hela)
)

obj.Hela <- obj.Hela %>%
  RunPCA() %>% FindNeighbors()
ElbowPlot(obj.Hela)

obj.Hela <- RunUMAP(
  obj.Hela,
  dims = 1:5,
  n.neighbors = 20,
  min.dist = 0.1
)
obj.Hela <- FindClusters(obj.Hela, resolution = 0.3)
obj.Hela$seurat_clusters <- as.numeric(obj.Hela$seurat_clusters)

pdf("Hela/UMAP.after.regression.pdf", 8, 6)
multiplot(
  DimPlot(obj.Hela, group.by = "seurat_clusters",
          cols = color.cluster) + coord_equal(1.3),
  DimPlot(obj.Hela, group.by = "Phase", cols = color.Phase) + coord_equal(1.3),
  DimPlot(obj.Hela, group.by = "Concentration",
          cols = color.concentration) + coord_equal(1.3),
  DimPlot(obj.Hela, group.by = "Time", cols = color.time) + coord_equal(1.3),
  cols = 2
)
dev.off()

# marker for each cluster
Idents(obj.Hela)
marker.Hela <- FindAllMarkers(
  obj.Hela,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
table(marker.Hela$cluster)
table(marker.Hela$p_val_adj < 0.05)
marker.Hela <- marker.Hela[marker.Hela$p_val_adj < 0.05, ]
pdf("Hela/Heatmap.Marker.cluster.pdf", 6, 5)
DoHeatmap(
  obj.Hela,
  features = marker.Hela$gene,
  label = T,
  draw.lines = F,
  group.colors = color.cluster,
  disp.min = -2.5,
  disp.max = 2.5,
  angle = 0
) +
  scale_fill_gradientn(colors = rev(brewer.pal(7, "RdYlBu")))
dev.off()

# GO for cluster markers
GO.Hela.C0 <- enrichGO(
  marker.Hela[marker.Hela$cluster == 0, ]$gene,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
GO.Hela.C0 <- simplify(GO.Hela.C0,
                       cutoff = 0.7,
                       by = "p.adjust",
                       select_fun = min)
GO.Hela.C1 <- enrichGO(
  marker.Hela[marker.Hela$cluster == 1, ]$gene,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
GO.Hela.C1 <- simplify(GO.Hela.C1,
                       cutoff = 0.7,
                       by = "p.adjust",
                       select_fun = min)
GO.Hela.C2 <- enrichGO(
  marker.Hela[marker.Hela$cluster == 2, ]$gene,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
GO.Hela.C2 <- simplify(GO.Hela.C2,
                       cutoff = 0.7,
                       by = "p.adjust",
                       select_fun = min)
write.csv(GO.Hela.C2@result,
          "GO.Hela.C2.csv")
pdf("Hela/Bar_GO_Clusters.pdf", 9, 4)
dotplot(GO.Hela.C0, showCategory = 20) + ggtitle("GO of C0 Makers")
dotplot(GO.Hela.C1, showCategory = 20) + ggtitle("GO of C1 Makers")
dotplot(GO.Hela.C2, showCategory = 20) + ggtitle("GO of C2 Makers")
dev.off()

pdf("Hela/Tree_GO_Clusters.pdf", 9, 5)
edox2 <- pairwise_termsim(GO.Hela.C0)
treeplot(edox2, showCategory = 30, xlim = c(0, 40)) + ggtitle("GO of C0 Makers")
edox2 <- pairwise_termsim(GO.Hela.C1)
treeplot(edox2, showCategory = 30, xlim = c(0, 40)) + ggtitle("GO of C1 Makers")
edox2 <- pairwise_termsim(GO.Hela.C2)
treeplot(edox2, showCategory = 30, xlim = c(0, 40)) + ggtitle("GO of C2 Makers")
dev.off()

View(GO.Hela.C1@result)
plot.data <-
  rbind(GO.Hela.C0@result[c(
    "GO:0006839",
    "GO:2000573",
    "GO:0043462",
    "GO:0050821",
    "GO:1904951",
    "GO:0007004",
    "GO:0002446"
  ), ],
  GO.Hela.C1@result[c(
    "GO:0016126",
    "GO:0006066",
    "GO:0009165",
    "GO:0006790",
    "GO:0034976",
    "GO:0009615",
    "GO:0006767"
  ), ],
  GO.Hela.C2@result[c(
    "GO:0000184",
    "GO:0006413",
    "GO:0042255",
    "GO:0006364",
    "GO:1903311",
    "GO:0006457",
    "GO:2001244"
  ), ])
plot.data$cluster <- rep(c(0, 1, 2), times = c(7, 7, 7))

frac2number <- function(x) {
  temp <- strsplit(x , "/")
  temp <-
    lapply(temp, function(x) {
      return(as.numeric(x[1]) / as.numeric(x[2]))
    })
  return(unlist(temp))
}
plot.data$Enrichment <-
  frac2number(plot.data$GeneRatio) / frac2number(plot.data$BgRatio)

plot.data <- plot.data[order(plot.data$Enrichment), ]
plot.data <- plot.data[order(plot.data$cluster, decreasing = T), ]
plot.data$Description <- factor(plot.data$Description,
                                levels = plot.data$Description)
pdf("Hela/Bar_GO_All.pdf", 10, 3.5)
ggplot(plot.data) +
  geom_bar(aes(
    x = Description,
    y = log10(Enrichment),
    fill = as.character(cluster)
  ),
  stat = "identity") +
  scale_fill_manual(values = color.cluster) +
  coord_flip() +
  theme_bw()
dev.off()

# isoform composition
features <- grep("percent", colnames(isoform), value = T)
for (one in features) {
  obj.Hela@meta.data[, one] <-
    isoform[rownames(obj.Hela@meta.data), one]
}

pdf("Hela/Vln.Isoform.classification.pdf", 10, 4)
VlnPlot(obj.Hela,
        features = features,
        cols = color.cluster,
        ncol = 5)
dev.off()

###  DTU analysis----------------------------
marker.Hela.C2 <- FindMarkers(
  object = obj.Hela,
  ident.1 = 2,
  ident.2 = c(0, 1),
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = F
)
marker.Hela.C2 <- marker.Hela.C2[marker.Hela.C2$p_val_adj < 0.05, ]

gene2trans <- gene2trans[, c(4, 3, 1, 2)]
rownames(gene2trans) <- gene2trans$trans.symbol

rc_trans.Hela <- rc_trans[, colnames(obj.Hela)]
table(rownames(rc_trans.Hela) %in% gene2trans$trans.symbol)
pb <- obj.Hela@meta.data
table(pb$seurat_clusters)
pb$seurat_clusters <- gsub("0|1", "C0_C1", pb$seurat_clusters)
pb$seurat_clusters <- gsub("2", "C2", pb$seurat_clusters)
pb$seurat_clusters_2 <- pb$seurat_clusters

DTU.Hela <- run_drimseq(
  counts = as.matrix(rc_trans.Hela),
  tx2gene = gene2trans,
  pd = pb,
  id_col = "Rename",
  cond_col = "seurat_clusters",
  cond_levels = c("C2", "C0_C1"),
  filtering_strategy = "own",
  filter_only = F,
  BPPARAM = BiocParallel::SnowParam(6),
  min_samps_gene_expr = 25,
  min_gene_expr = 5,
  min_samps_feature_expr = 25,
  min_feature_expr = 5
)
DTU.Hela <-
  posthoc_and_stager(dturtle = DTU.Hela,
                     ofdr = 0.05,
                     posthoc = 0.1)
DTU.Hela$meta_table_gene[1:5, ]
DTU.Hela$meta_table_tx[1:5, ]
length(DTU.Hela$sig_gene)
length(DTU.Hela$sig_tx)

pdf("Hela/DGE_DTU.pdf", 4, 3)
plot(Venn(list(
  "DGE" = rownames(marker.Hela.C2),
  "DTU" = DTU.Hela$sig_gene
)),
doWeights = TRUE,
show = list(SetLabels = T, Faces = FALSE))
dev.off()

DTU.Hela <- create_dtu_table(dturtle = DTU.Hela,
                             add_tx_metadata = list(
                               tx_expr_in_max = c("exp_in", max),
                               gene.id = c("gene.id", function(x) {
                                 return(names(which.max(table(x))))
                               })
                             ))
head(DTU.Hela$dtu_table, n = 20)

temp <-
  plot_proportion_barplot(dturtle = DTU.Hela, genes = "TFAP2A")
temp$TFAP2A
temp <-
  plot_proportion_pheatmap(
    dturtle = DTU.Hela,
    genes = "TFAP2A",
    sample_meta_table_columns = c("sample_id", "seurat_clusters_2"),
    include_expression = TRUE,
    treeheight_col = 20,
    BPPARAM = BiocParallel::SnowParam(1)
  )
temp$TFAP2A
plot_transcripts_view(
  dturtle = DTU.Hela,
  genes = "TFAP2A",
  gtf = gtf,
  genome = NULL,
  reduce_introns_min_size = 100,
  arrow_colors = c("#f46d43", "#74add1"),
  one_to_one = TRUE,
  savepath = "./Hela",
  filename_ext  = "_transcripts.pdf",
  width = 5,
  height = 3
)

# functional annotation of DTU
head(DTU.Hela$dtu_table)
head(DTU.Hela$meta_table_gene)
head(DTU.Hela$meta_table_tx)

DTU.Hela.sig <-
  DTU.Hela$meta_table_tx[DTU.Hela$meta_table_tx$tx %in%
                           DTU.Hela$sig_tx, ]
DTU.Hela.sig <- DTU.Hela.sig[!DTU.Hela.sig$gene %in%
                               rownames(marker.Hela.C2), ] # 390
length(unique(DTU.Hela.sig$gene))
table(table(DTU.Hela.sig$gene))
gene.selected <- table(DTU.Hela.sig$gene)
gene.selected <- gene.selected[gene.selected >= 2]
gene.selected <- names(gene.selected)

DTU.Hela.sig.2 <- data.frame()
for (one in gene.selected) {
  one <- DTU.Hela.sig[DTU.Hela.sig$gene == one, ]
  trans_C2 <- one$tx[which.max(one$exp_in_C2)]
  trans_C01 <- one$tx[which.max(one$exp_in_C0_C1)]
  temp <- data.frame(
    row.names = one$gene[1],
    "Gene" = one$gene[1],
    "Isoform_C2" = trans_C2,
    "Isoform_C0_C1" = trans_C01,
    "Percent_C2" = one[trans_C2, ]$exp_in_C2,
    "Percent_C0_C1" = one[trans_C01, ]$exp_in_C0_C1,
    "Type_C2" = trans.info.hg38[trans_C2, ]$`Transcript type`,
    "Type_C0_C1" = trans.info.hg38[trans_C01, ]$`Transcript type`,
    "Protein_C2" = protien.info.hg38[trans_C2, ]$`Protein stable ID`,
    "Protein_C0_C1" = protien.info.hg38[trans_C01, ]$`Protein stable ID`,
    "Location_C2" = protien.info.hg38[trans_C2, ]$location,
    "Location_C0_C1" = protien.info.hg38[trans_C01, ]$location
  )
  DTU.Hela.sig.2 <- rbind(DTU.Hela.sig.2, temp)
}

temp <-
  DTU.Hela.sig.2[DTU.Hela.sig.2$Type_C2 != DTU.Hela.sig.2$Type_C0_C1, ]
write.csv(temp, "Hela/DTU.biotype_change.csv")

temp <- DTU.Hela.sig.2[DTU.Hela.sig.2$Type_C2 == "protein_coding" &
                         DTU.Hela.sig.2$Type_C0_C1 == "protein_coding", ]
temp <- temp[temp$Protein_C0_C1 != temp$Protein_C2, ]
write.csv(temp, "Hela/DTU.protein_change.csv")

temp <- DTU.Hela.sig.2[DTU.Hela.sig.2$Type_C2 == "protein_coding" &
                         DTU.Hela.sig.2$Type_C0_C1 == "protein_coding", ]
temp <- temp[temp$Location_C0_C1 != temp$Location_C2, ]
write.csv(temp, "Hela/DTU.subcelluar_location_change.csv")

## 2.2. HepG2 ---------------------------------------
dir.create("HepG2")
obj.HepG2 <- obj.HepG2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>% FindNeighbors()
ElbowPlot(obj.HepG2)

obj.HepG2 <- RunUMAP(
  obj.HepG2,
  dims = 1:5,
  n.neighbors = 40,
  min.dist = 0.05
)
obj.HepG2 <- FindClusters(obj.HepG2, resolution = 0.2)

pdf("HepG2/UMAP.befor.regression.pdf", 8, 6)
multiplot(
  DimPlot(obj.HepG2, group.by = "Gene_snn_res.0.2",
          cols = color.cluster) + coord_equal(1.3),
  DimPlot(obj.HepG2, group.by = "Phase", cols = color.Phase) + coord_equal(1.3),
  DimPlot(obj.HepG2, group.by = "Concentration",
          cols = color.concentration) + coord_equal(1.3),
  DimPlot(obj.HepG2, group.by = "Time", cols = color.time) + coord_equal(1.3),
  cols = 2
)
dev.off()

# regress cell cycle
obj.HepG2 <- ScaleData(
  obj.HepG2,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(obj.HepG2)
)

obj.HepG2 <- obj.HepG2 %>%
  RunPCA() %>% FindNeighbors()
ElbowPlot(obj.HepG2)

obj.HepG2 <- RunUMAP(
  obj.HepG2,
  dims = 1:5,
  n.neighbors = 30,
  min.dist = 0.1
)
obj.HepG2 <- FindClusters(obj.HepG2, resolution = 0.2)

# swap cluster
temp <- as.character(obj.HepG2$Gene_snn_res.0.2)
temp[temp == "0"] <- "4"
temp[temp == "1"] <- "0"
temp[temp == "4"] <- "1"
table(temp)
table(obj.HepG2$Gene_snn_res.0.2)
obj.HepG2$Gene_snn_res.0.2 <- as.numeric(temp)
Idents(obj.HepG2) <- obj.HepG2$Gene_snn_res.0.2
Idents(obj.HepG2) <- factor(Idents(obj.HepG2), levels = c(0, 1, 2))

pdf("HepG2/UMAP.after.regression.pdf", 8, 6)
multiplot(
  DimPlot(obj.HepG2, group.by = "Gene_snn_res.0.2",
          cols = color.cluster) + coord_equal(1.1),
  DimPlot(obj.HepG2, group.by = "Phase", cols = color.Phase) + coord_equal(1.1),
  DimPlot(obj.HepG2, group.by = "Concentration",
          cols = color.concentration) + coord_equal(1.1),
  DimPlot(obj.HepG2, group.by = "Time", cols = color.time) + coord_equal(1.1),
  cols = 2
)
dev.off()

Idents(obj.HepG2)
marker.HepG2 <- FindAllMarkers(
  obj.HepG2,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
table(marker.HepG2$cluster)
table(marker.HepG2$p_val_adj < 0.05)
marker.HepG2 <- marker.HepG2[marker.HepG2$p_val_adj < 0.05, ]
pdf("HepG2/Heatmap.Marker.cluster.pdf", 6, 5)
DoHeatmap(
  obj.HepG2,
  features = marker.HepG2$gene,
  label = T,
  draw.lines = F,
  group.colors = color.cluster,
  disp.min = -2.5,
  disp.max = 2.5,
  angle = 0
) +
  scale_fill_gradientn(colors = rev(brewer.pal(7, "RdYlBu")))
dev.off()

GO.HepG2.C0 <- enrichGO(
  marker.HepG2[marker.HepG2$cluster == 0, ]$gene,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
GO.HepG2.C0 <- simplify(GO.HepG2.C0,
                        cutoff = 0.7,
                        by = "p.adjust",
                        select_fun = min)
GO.HepG2.C1 <- enrichGO(
  marker.HepG2[marker.HepG2$cluster == 1, ]$gene,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
GO.HepG2.C1 <- simplify(GO.HepG2.C1,
                        cutoff = 0.7,
                        by = "p.adjust",
                        select_fun = min)
GO.HepG2.C2 <- enrichGO(
  marker.HepG2[marker.HepG2$cluster == 2, ]$gene,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
GO.HepG2.C2 <- simplify(GO.HepG2.C2,
                        cutoff = 0.7,
                        by = "p.adjust",
                        select_fun = min)
write.csv(GO.HepG2.C2@result,
          "GO.HepG2.C2.csv")
pdf("HepG2/Bar_GO_Clusters.pdf", 9, 4)
dotplot(GO.HepG2.C0, showCategory = 20) + ggtitle("GO of C0 Makers")
dotplot(GO.HepG2.C1, showCategory = 20) + ggtitle("GO of C1 Makers")
dotplot(GO.HepG2.C2, showCategory = 20) + ggtitle("GO of C2 Makers")
dev.off()

pdf("HepG2/Tree_GO_Clusters.C3.pdf", 9, 5)
edox2 <- pairwise_termsim(GO.HepG2.C0)
treeplot(edox2, showCategory = 30, xlim = c(0, 40)) + ggtitle("GO of C0 Makers")
edox2 <- pairwise_termsim(GO.HepG2.C1)
treeplot(edox2, showCategory = 30, xlim = c(0, 40)) + ggtitle("GO of C1 Makers")
edox2 <- pairwise_termsim(GO.HepG2.C2)
treeplot(edox2, showCategory = 30, xlim = c(0, 40)) + ggtitle("GO of C2 Makers")
dev.off()

View(GO.HepG2.C1@result)
plot.data <-
  rbind(GO.HepG2.C0@result[c(
    "GO:0022613",
    "GO:0002446",
    "GO:0140053",
    "GO:0072655",
    "GO:0050821",
    "GO:1902903",
    "GO:0070125"
  ), ],
  GO.HepG2.C1@result[c(
    "GO:0008202",
    "GO:0006066",
    "GO:0046395",
    "GO:0044282",
    "GO:0006790",
    "GO:1901605",
    "GO:0098754"
  ), ],
  GO.HepG2.C2@result[c(
    "GO:0000184",
    "GO:0006413",
    "GO:0042255",
    "GO:0006364",
    "GO:0022618",
    "GO:0030968",
    "GO:0006986"
  ), ])
plot.data$cluster <- rep(c(0, 1, 2), times = c(7, 7, 7))

plot.data$Enrichment <-
  frac2number(plot.data$GeneRatio) / frac2number(plot.data$BgRatio)

plot.data <- plot.data[order(plot.data$Enrichment), ]
plot.data <- plot.data[order(plot.data$cluster, decreasing = T), ]
plot.data$Description <- factor(plot.data$Description,
                                levels = plot.data$Description)
pdf("HepG2/Bar_GO_All.pdf", 10, 3.5)
ggplot(plot.data) +
  geom_bar(aes(
    x = Description,
    y = log10(Enrichment),
    fill = as.character(cluster)
  ),
  stat = "identity") +
  scale_fill_manual(values = color.cluster) +
  coord_flip() +
  theme_bw()
dev.off()

# isoform classification
for (one in features) {
  obj.HepG2@meta.data[, one] <-
    isoform[rownames(obj.HepG2@meta.data), one]
}

pdf("HepG2/Vln.Isoform.classification.pdf", 10, 4)
VlnPlot(obj.HepG2,
        features = features,
        cols = color.cluster,
        ncol = 5)
dev.off()

# DTU analysis
marker.HepG2.C2 <- FindMarkers(
  object = obj.HepG2,
  ident.1 = 2,
  ident.2 = c(0, 1),
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = F
)

rc_trans.HepG2 <- rc_trans[, colnames(obj.HepG2)]
table(rownames(rc_trans.HepG2) %in% gene2trans$trans.symbol)
pb <- obj.HepG2@meta.data
table(pb$seurat_clusters)
pb$seurat_clusters <- gsub("0|1", "C0_C1", pb$seurat_clusters)
pb$seurat_clusters <- gsub("2", "C2", pb$seurat_clusters)
pb$seurat_clusters_2 <- pb$seurat_clusters

DTU.HepG2 <- run_drimseq(
  counts = as.matrix(rc_trans.HepG2),
  tx2gene = gene2trans,
  pd = pb,
  id_col = "Rename",
  cond_col = "seurat_clusters",
  cond_levels = c("C2", "C0_C1"),
  filtering_strategy = "own",
  filter_only = F,
  BPPARAM = BiocParallel::SnowParam(6),
  min_samps_gene_expr = 25,
  min_gene_expr = 5,
  min_samps_feature_expr = 25,
  min_feature_expr = 5
)

DTU.HepG2 <-
  posthoc_and_stager(dturtle = DTU.HepG2,
                     ofdr = 0.05,
                     posthoc = 0.1)
DTU.HepG2$meta_table_gene[1:5, ]
DTU.HepG2$meta_table_tx[1:5, ]
length(DTU.HepG2$sig_gene)
length(DTU.HepG2$sig_tx)

pdf("HepG2/DGE_DTU.pdf", 4, 3)
plot(Venn(list(
  "DGE" = rownames(marker.HepG2.C2),
  "DTU" = DTU.HepG2$sig_gene
)),
doWeights = TRUE,
show = list(SetLabels = T, Faces = FALSE))
dev.off()

DTU.HepG2 <- create_dtu_table(dturtle = DTU.HepG2,
                              add_tx_metadata = list(
                                tx_expr_in_max = c("exp_in", max),
                                gene.id = c("gene.id", function(x) {
                                  return(names(which.max(table(x))))
                                })
                              ))
head(DTU.HepG2$dtu_table, n = 20)

temp <-
  plot_proportion_barplot(dturtle = DTU.HepG2, genes = "EIF4G1")
temp$EIF4G1
temp <-
  plot_proportion_pheatmap(
    dturtle = DTU.HepG2,
    genes = "EIF4G1",
    sample_meta_table_columns = c("sample_id", "seurat_clusters_2"),
    include_expression = TRUE,
    treeheight_col = 20,
    BPPARAM = BiocParallel::SnowParam(1)
  )
temp$EIF4G1
plot_transcripts_view(
  dturtle = DTU.HepG2,
  genes = "EIF4G1",
  gtf = gtf,
  genome = NULL,
  reduce_introns_min_size = 100,
  arrow_colors = c("#f46d43", "#74add1"),
  one_to_one = TRUE,
  savepath = "./HepG2",
  filename_ext  = "_transcripts.pdf",
  width = 5,
  height = 3
)

# 3. compare different cell lines -----------------------------------------
dir.create("Compare_Cell_Line")
library(Vennerable)

marker.Hela.C2

## DEG and DTU genes --------------------------
plot(Venn(
  list("Hela C0" = marker.Hela[marker.Hela$cluster == 0, ]$gene,
       "HepG2 C0" = marker.HepG2[marker.HepG2$cluster == 0, ]$gene)
))
plot(Venn(
  list("Hela C1" = marker.Hela[marker.Hela$cluster == 1, ]$gene,
       "HepG2 C1" = marker.HepG2[marker.HepG2$cluster == 1, ]$gene)
))
plot(Venn(
  list("Hela C2" = marker.Hela[marker.Hela$cluster == 2, ]$gene,
       "HepG2 C2" = marker.HepG2[marker.HepG2$cluster == 2, ]$gene)
))

plot(Venn(list(
  "Hela C2" = rownames(marker.Hela.C2)[marker.Hela.C2$avg_log2FC > 0],
  "HepG2 C2" = rownames(marker.HepG2.C2)[marker.HepG2.C2$avg_log2FC >
                                           0]
)))

pdf("Compare_Cell_Line/Venn.Marker.C2.pdf", 4, 3)
# plot(Venn(list("Hela" = marker.Hela[marker.Hela$cluster==2,]$gene,
#                "HepG2" = marker.HepG2[marker.HepG2$cluster==2,]$gene)),
#      doWeights = TRUE,
#      show = list( SetLabels = T, Faces = FALSE))
plot(Venn(list(
  "Hela" = rownames(marker.Hela.C2)[marker.Hela.C2$avg_log2FC > 0],
  "HepG2" = rownames(marker.HepG2.C2)[marker.HepG2.C2$avg_log2FC >
                                        0]
)),
doWeights = TRUE,
show = list(SetLabels = T, Faces = FALSE))
plot(Venn(list(
  "Hela" = rownames(marker.Hela.C2)[marker.Hela.C2$avg_log2FC < 0],
  "HepG2" = rownames(marker.HepG2.C2)[marker.HepG2.C2$avg_log2FC <
                                        0]
)),
doWeights = TRUE,
show = list(SetLabels = T, Faces = FALSE))
dev.off()

pdf("Compare_Cell_Line/Venn.DTU.pdf", 4, 3)
DTU.up.Hela <-
  DTU.Hela$sig_tx[DTU.Hela$meta_table_tx[DTU.Hela$sig_tx, ]$exp_in_C2
                  > DTU.Hela$meta_table_tx[DTU.Hela$sig_tx, ]$exp_in_C0_C1]
DTU.down.Hela <-
  DTU.Hela$sig_tx[DTU.Hela$meta_table_tx[DTU.Hela$sig_tx, ]$exp_in_C2
                  < DTU.Hela$meta_table_tx[DTU.Hela$sig_tx, ]$exp_in_C0_C1]
DTU.up.HepG2 <-
  DTU.HepG2$sig_tx[DTU.HepG2$meta_table_tx[DTU.HepG2$sig_tx, ]$exp_in_C2
                   > DTU.HepG2$meta_table_tx[DTU.HepG2$sig_tx, ]$exp_in_C0_C1]
DTU.down.HepG2 <-
  DTU.HepG2$sig_tx[DTU.HepG2$meta_table_tx[DTU.HepG2$sig_tx, ]$exp_in_C2
                   < DTU.HepG2$meta_table_tx[DTU.HepG2$sig_tx, ]$exp_in_C0_C1]
plot(Venn(list("Hela" = DTU.up.Hela,
               "HepG2" = DTU.up.HepG2)),
     doWeights = TRUE,
     show = list(SetLabels = T, Faces = FALSE))
plot(Venn(list("Hela" = DTU.down.Hela,
               "HepG2" = DTU.down.HepG2)),
     doWeights = TRUE,
     show = list(SetLabels = T, Faces = FALSE))
dev.off()

plot(Venn(list(
  "DEG" = rownames(marker.Hela.C2),
  "DTU" = DTU.Hela$sig_gene
)),
doWeights = TRUE,
show = list(SetLabels = T, Faces = FALSE))
plot(Venn(list(
  "DEG" = rownames(marker.HepG2.C2),
  "DTU" = DTU.HepG2$sig_gene
)),
doWeights = TRUE,
show = list(SetLabels = T, Faces = FALSE))

## go terms ---------------------
common.term <- intersect(GO.Hela.C2@result$ID,
                         GO.HepG2.C2@result$ID)
View(GO.Hela.C2@result)

marker.Hela$Cell_Line <- "Hela"
marker.HepG2$Cell_Line <- "HepG2"

plot.data.count <- data.frame(
  row.names = common.term,
  "Decription" = GO.Hela.C2@result[common.term, ]$Description,
  "Hela_0" = 1000,
  "Hela_1" = 1000,
  "Hela_2" = 1000,
  "HepG2_0" = 1000,
  "HepG2_1" = 1000,
  "HepG2_2" = 1000
)
plot.data.OR <- data.frame(
  row.names = common.term,
  "Decription" = GO.Hela.C2@result[common.term, ]$Description,
  "Hela_0" = 1000,
  "Hela_1" = 1000,
  "Hela_2" = 1000,
  "HepG2_0" = 1000,
  "HepG2_1" = 1000,
  "HepG2_2" = 1000
)

for (one in common.term) {
  gene.one <- unique(
    select(
      org.Hs.eg.db,
      keys = org.Hs.egGO2ALLEGS[[one]],
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "ENTREZID"
    )$SYMBOL
  )
  ratio.bg <- frac2number(GO.Hela.C2@result[one, ]$BgRatio)

  # hela c0
  markers <- marker.Hela[marker.Hela$cluster == 0, ]$gene
  hits <- length(intersect(markers, gene.one))
  plot.data.count[one, ]$Hela_0 <- hits
  plot.data.OR[one, ]$Hela_0 <- (hits / length(markers)) / ratio.bg
  # hela c1
  markers <- marker.Hela[marker.Hela$cluster == 1, ]$gene
  hits <- length(intersect(markers, gene.one))
  plot.data.count[one, ]$Hela_1 <- hits
  plot.data.OR[one, ]$Hela_1 <- (hits / length(markers)) / ratio.bg
  # hela c2
  markers <- marker.Hela[marker.Hela$cluster == 2, ]$gene
  hits <- length(intersect(markers, gene.one))
  plot.data.count[one, ]$Hela_2 <- hits
  plot.data.OR[one, ]$Hela_2 <- (hits / length(markers)) / ratio.bg

  # HepG2 c0
  markers <- marker.HepG2[marker.HepG2$cluster == 0, ]$gene
  hits <- length(intersect(markers, gene.one))
  plot.data.count[one, ]$HepG2_0 <- hits
  plot.data.OR[one, ]$HepG2_0 <- (hits / length(markers)) / ratio.bg
  # HepG2 c1
  markers <- marker.HepG2[marker.HepG2$cluster == 1, ]$gene
  hits <- length(intersect(markers, gene.one))
  plot.data.count[one, ]$HepG2_1 <- hits
  plot.data.OR[one, ]$HepG2_1 <- (hits / length(markers)) / ratio.bg
  # HepG2 c2
  markers <- marker.HepG2[marker.HepG2$cluster == 2, ]$gene
  hits <- length(intersect(markers, gene.one))
  plot.data.count[one, ]$HepG2_2 <- hits
  plot.data.OR[one, ]$HepG2_2 <- (hits / length(markers)) / ratio.bg
}

common.term <-
  common.term[order(plot.data.OR$Hela_2, decreasing = T)]
plot.data.count <- melt(plot.data.count)
plot.data.OR <- melt(plot.data.OR)

plot.data <- plot.data.OR
colnames(plot.data) <- c("Description", "Cluster", "Enrichment")
plot.data$Count <- plot.data.count$value
rm(plot.data.count, plot.data.OR)
plot.data$Cell_Line <- gsub("_.+?$", "", plot.data$Cluster)
plot.data$Cluster_1 <- gsub("^.+?_", "", plot.data$Cluster)
plot.data <- plot.data[order(plot.data$Enrichment), ]
plot.data$Description <- factor(plot.data$Description,
                                levels = GO.Hela.C2@result[common.term, ]$Description)

plot.data.GO <- plot.data

pdf("Compare_Cell_Line/Dot.Compare_GO.pdf", 10, 4)
ggplot(plot.data) +
  geom_point(aes(
    x = paste0("C", Cluster_1),
    y = Description,
    size = Count,
    color = log10(Enrichment + 0.1)
  )) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  facet_wrap( ~ Cell_Line) +
  theme_lzy()
dev.off()

rm(p1, p2, p3)

## percent for NIC merge plot -----------------------
plot.data <- rbind(obj.Hela@meta.data, obj.HepG2@meta.data)
plot.data <-
  plot.data[, c(
    "Cell_Line",
    "seurat_clusters",
    "percent.NIC.CJ",
    "percent.NIC.CS",
    "percent.NIC.IR"
  )]
colnames(plot.data)[3:5] <- c("CJ", "CS", "IR")
plot.data <- melt(plot.data)
plot.data$Seurat_Clusters <-
  as.numeric(unfactor(plot.data$seurat_clusters)) + 1
plot.data$Seurat_Clusters <-
  paste("C", plot.data$Seurat_Clusters, sep = "")

pdf("Compare_Cell_Line/Vln.Percent_NIC.pdf", 4, 3)
ggplot(plot.data, aes(x = Seurat_Clusters, y = value * 100)) +
  geom_violin(aes(fill = Seurat_Clusters), show.legend = F) +
  geom_jitter(width = 0.3,
              height = 0,
              size = 0.00001) +
  scale_fill_manual(values = color.cluster) +
  facet_grid(rows = vars(variable),
             cols = vars(Cell_Line)) +
  ggsignif::geom_signif(
    comparisons = list(c("C3", "C1"),
                       c("C3", "C2")),
    map_signif_level = T,
    textsize = 1.5,
    test = wilcox.test,
    step_increase = 0.2,
    margin_top = 0.15
  ) +
  theme_lzy() + ylab("Percent of NIC") + ylim(c(0, 105))
dev.off()

c2 <- plot.data[plot.data$Cell_Line == "HepG2" &
                  plot.data$Seurat_Clusters == "C2" &
                  plot.data$variable == "CJ", ]$value
c3 <- plot.data[plot.data$Cell_Line == "HepG2" &
                  plot.data$Seurat_Clusters == "C3" &
                  plot.data$variable == "CJ", ]$value
c1 <- plot.data[plot.data$Cell_Line == "HepG2" &
                  plot.data$Seurat_Clusters == "C1" &
                  plot.data$variable == "CJ", ]$value
a <- wilcox.test(c2, c1)
a$p.value

## examples for DTU --------------
DTU.common <- intersect(DTU.Hela$sig_gene, DTU.HepG2$sig_gene)

CC.gene <- read.table("CellCycleGene.txt", sep = "\t")
CC.gene <- unlist(strsplit(CC.gene$V1, ","))
CC.gene <- CC.gene[CC.gene != ""]
CC.gene <- intersect(CC.gene, DTU.common)

plot_transcripts_view(
  dturtle = DTU.Hela,
  genes = CC.gene,
  gtf = gtf,
  genome = NULL,
  reduce_introns_min_size = 100,
  arrow_colors = c("#f46d43", "#74add1"),
  one_to_one = TRUE,
  savepath = "./Compare_Cell_Line/CellCycle",
  filename_ext  = "_Hela.transcripts.pdf",
  width = 5,
  height = 3
)
plot_transcripts_view(
  dturtle = DTU.HepG2,
  genes = CC.gene,
  gtf = gtf,
  genome = NULL,
  reduce_introns_min_size = 100,
  arrow_colors = c("#f46d43", "#74add1"),
  one_to_one = TRUE,
  savepath = "./Compare_Cell_Line/CellCycle",
  filename_ext  = "_HepG2.transcripts.pdf",
  width = 5,
  height = 3
)

Liver_cancer_gene <-
  read.table("frequently-mutated-genes.2022-04-14.tsv",
             sep = "\t",
             header = T)
liver_gene_selected <- intersect(DTU.HepG2$sig_gene,
                                 Liver_cancer_gene[1:100, ]$Symbol)
plot_transcripts_view(
  dturtle = DTU.HepG2,
  genes = liver_gene_selected,
  gtf = gtf,
  genome = NULL,
  reduce_introns_min_size = 100,
  arrow_colors = c("#f46d43", "#74add1"),
  one_to_one = TRUE,
  savepath = "./Compare_Cell_Line",
  filename_ext  = "_HepG2.transcripts.pdf",
  width = 5,
  height = 3
)


# 4. supplementary tables -------------------------------------------------
sample.info.Hela <- obj.Hela@meta.data
sample.info.HepG2 <- obj.HepG2@meta.data

sample.info.Hela$seurat_clusters <-
  paste("C", as.numeric(unfactor(sample.info.Hela$seurat_clusters)) + 1,
        sep = "")
sample.info.HepG2$seurat_clusters <-
  paste("C", as.numeric(unfactor(sample.info.HepG2$seurat_clusters)) + 1,
        sep = "")
write.table(sample.info.Hela,
            "sample.info.Hela.tsv", sep = "\t")
write.table(sample.info.HepG2,
            "sample.info.HepG2.tsv", sep = "\t")


marker.Hela$cluster <-
  paste("C", as.numeric(unfactor(marker.Hela$cluster)) + 1,
        sep = "")
marker.HepG2$cluster <-
  paste("C", as.numeric(unfactor(marker.HepG2$cluster)) + 1,
        sep = "")
write.table(marker.Hela,
            "marker.Hela.tsv", sep = "\t")
write.table(marker.HepG2,
            "marker.HepG2.tsv", sep = "\t")

FDR <- DTU.Hela$FDR_table
rownames(FDR) <- FDR$txID
out <- DTU.Hela$meta_table_tx[DTU.Hela$meta_table_tx$tx %in%
                                DTU.Hela$sig_tx, ]
out$gene_FDR_value <- FDR[out$tx, ]$gene
out$tx_FDR_value <- FDR[out$tx, ]$transcript
out$Transcript_Type <- trans.info.hg38[out$tx, ]$`Transcript type`
write.table(out,
            "DTU.Hela.tsv", sep = "\t")

FDR <- DTU.HepG2$FDR_table
rownames(FDR) <- FDR$txID
out <- DTU.HepG2$meta_table_tx[DTU.HepG2$meta_table_tx$tx %in%
                                 DTU.HepG2$sig_tx, ]
out$gene_FDR_value <- FDR[out$tx, ]$gene
out$tx_FDR_value <- FDR[out$tx, ]$transcript
out$Transcript_Type <- trans.info.hg38[out$tx, ]$`Transcript type`
write.table(out,
            "DTU.HepG2.tsv", sep = "\t")
