setwd("E:/LabWork/Project/SCANSeq2/merge/Part1/")

library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(data.table)

source("E:/LabWork/code/MyFunction.R")
ReadGeneInfo("hg38", classify.protein = T)
ReadGeneInfo("mm10")

# load data
rc_gene_9CL <- readRDS("../../batch/1 9CL/data/rc_gene.RDS")
rc_gene_4CL <- readRDS("../../batch/3 4CL/data/rc_gene.RDS")
rc_gene_100 <- readRDS("../../batch/4 UMI_100/data/rc_gene.RDS")
rc_gene_200 <- readRDS("../../batch/5 UMI_200/data/rc_gene.RDS")
rc_mES <- readRDS("../../Test/SCANSeq_test/rc_gene_mES.RDS")

sample.info.9CL <- read.table("../../batch/1 9CL/sample.info.txt")
sample.info.9Mix <-
  read.table("../../batch/2 9CL_Mix/sample.info.txt")
sample.info.4CL <- read.table("../../batch/3 4CL/sample.info.txt")
sample.info.100 <-
  read.table("../../batch/4 UMI_100/sample.info.txt")
sample.info.200 <-
  read.table("../../batch/5 UMI_200/sample.info.txt")

sample.info.9CL$Cell_Line <- factor(
  sample.info.9CL$Cell_Line,
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
sample.info.9Mix$Cell_Line <- factor(
  sample.info.9Mix$Cell_Line,
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
sample.info.4CL$Cell_Line <- factor(
  sample.info.4CL$Cell_Line,
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
sample.info.100$Cell_Line <- factor(
  sample.info.100$Cell_Line,
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
sample.info.200$Cell_Line <- factor(
  sample.info.200$Cell_Line,
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

sample.info.9Mix2 <-
  read.table("../../archive/20211026_1/sample.info.txt")
sample.info.9Mix2 <-
  sample.info.9Mix2[sample.info.9Mix2$Pass_QC == 1, ]
sample.info.9Mix2$Cell_Line <- factor(
  sample.info.9Mix2$Cell_Line,
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
sample.info.9Mix2$Library <- "9CL_Mix2"
sample.info.9Mix3 <-
  read.table("../../archive/20210904/sample.info.txt")
sample.info.9Mix3 <-
  sample.info.9Mix3[sample.info.9Mix3$Pass_QC == 1, ]
sample.info.9Mix3$Cell_Line <- factor(
  sample.info.9Mix3$Cell_Line,
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
sample.info.9Mix3$Library <- "9CL_Mix3"

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


# 1. demultiplexing -------------------------------------------------------
TotalReads <- 84801426
plot.data <- aggregate(sample.info.9CL$Raw_Reads,
                       list(sample.info.9CL$Cell_Line),
                       FUN = sum)
colnames(plot.data) <- c("Cell_Line", "Reads")
plot.data$Cell_Line <- as.character(plot.data$Cell_Line)
plot.data <- plot.data[order(plot.data$Reads, decreasing = T), ]
plot.data$Cell_Line <-
  factor(plot.data$Cell_Line, levels = plot.data$Cell_Line)
temp <- t(as.data.frame(c(
  "unclassified",
  TotalReads - sum(plot.data$Reads)
)))
colnames(temp) <- c("Cell_Line", "Reads")
plot.data <- rbind(plot.data, temp)
rownames(plot.data) <- plot.data$Cell_Line
plot.data$Freq <- as.numeric(plot.data$Reads) / TotalReads
plot.data$label <- paste(plot.data$Cell_Line,
                         "\n",
                         round(plot.data$Freq, 3) * 100,
                         "%",
                         sep = "")

pdf("Pie.demultiplex.pdf", 6, 6)
ggpubr::ggdonutchart(
  plot.data,
  "Freq",
  label = "label",
  fill = "Cell_Line",
  color = "white",
  lab.pos = "in",
  palette = c(color.cell, "unclassified" = "gray75")
)
dev.off()

sample.info.9CL <- sample.info.9CL[sample.info.9CL$Pass_QC == 1,]
sample.info.9Mix <-
  sample.info.9Mix[sample.info.9Mix$Pass_QC == 1,]
sample.info.4CL <- sample.info.4CL[sample.info.4CL$Pass_QC == 1,]
sample.info.100 <- sample.info.100[sample.info.100$Pass_QC == 1,]
sample.info.200 <- sample.info.200[sample.info.200$Pass_QC == 1,]
sample.info.9Mix2 <-
  sample.info.9Mix2[sample.info.9Mix2$Pass_QC == 1,]
sample.info.9Mix3 <-
  sample.info.9Mix3[sample.info.9Mix3$Pass_QC == 1,]

# 2. Sensitivity ----------------------------------------------------------
dir.create("Sensitivity")

## 1.1 Number of Detected Genes and Isoforms---------------
pdf("Sensitivity/Vln_Gene_Trans.9CL.pdf", 6, 5)
temp <-
  aggregate(sample.info.9CL$Mapped_percent,
            list(sample.info.9CL$Cell_Line),
            median)
p1 <- ggplot(sample.info.9CL) +
  geom_violin(
    aes(x = Cell_Line, y = Mapped_percent * 100, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Mapped_percent * 100, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  ylab("Mapping Rate (%)") +
  scale_color_manual(values = color.cell) +
  scale_fill_manual(values = color.cell) +
  scale_x_discrete(breaks = temp$Group.1,
                   labels = paste0(temp$Group.1, "\n", signif(temp$x * 100, 2), "%")) +
  theme_lzy() + xlab(NULL)
temp <-
  aggregate(sample.info.9CL$Gene_Detected,
            list(sample.info.9CL$Cell_Line),
            median)
p2 <- ggplot(sample.info.9CL) +
  geom_violin(
    aes(x = Cell_Line, y = Gene_Detected / 1000, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Gene_Detected / 1000, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  ylab("Getected Genes (10^3)") +
  scale_color_manual(values = color.cell) +
  scale_fill_manual(values = color.cell) +
  scale_x_discrete(breaks = temp$Group.1,
                   labels = paste0(temp$Group.1, "\n", signif(temp$x, 4))) +
  theme_lzy() + xlab(NULL)
temp <-
  aggregate(sample.info.9CL$Trans_Detected,
            list(sample.info.9CL$Cell_Line),
            median)
p3 <- ggplot(sample.info.9CL) +
  geom_violin(
    aes(x = Cell_Line, y = Trans_Detected / 1000, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Trans_Detected / 1000, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  ylab("Getected Isoforms (10^3)") +
  scale_color_manual(values = color.cell) +
  scale_fill_manual(values = color.cell) +
  scale_x_discrete(breaks = temp$Group.1,
                   labels = paste0(temp$Group.1, "\n", signif(temp$x, 4))) +
  theme_lzy() + xlab(NULL)
multiplot(p1, p2, p3, cols = 1)
dev.off()

pdf("Sensitivity/Vln_Gene_Trans.9CL.pdf", 6, 3)
multiplot(p2, p3, cols = 1)
dev.off()
pdf("Sensitivity/Vln_Mappint.9CL.pdf", 6, 1.5)
plot(p1)
dev.off()

pdf("Sensitivity/Vln_Gene_Trans.4CL.pdf", 4, 5)
p1 <- ggplot(sample.info.4CL) +
  geom_violin(
    aes(x = Cell_Line, y = Mapped_percent * 100, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Mapped_percent * 100, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  ylab("Mapping Rate (%)") +
  scale_color_manual(values = color.cell) +
  scale_fill_manual(values = color.cell) +
  theme_lzy()
p2 <- ggplot(sample.info.4CL) +
  geom_violin(
    aes(x = Cell_Line, y = Gene_Detected / 1000, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Gene_Detected / 1000, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  ylab("Getected Genes (10^3)") +
  scale_color_manual(values = color.cell) +
  scale_fill_manual(values = color.cell) +
  theme_lzy()
p3 <- ggplot(sample.info.4CL) +
  geom_violin(
    aes(x = Cell_Line, y = Trans_Detected / 1000, fill = Cell_Line),
    alpha = 0.7,
    show.legend = F
  ) +
  geom_boxplot(
    aes(x = Cell_Line, y = Trans_Detected / 1000, color = Cell_Line),
    alpha = 0,
    show.legend = F
  ) +
  ylab("Getected Isoforms (10^3)") +
  scale_color_manual(values = color.cell) +
  scale_fill_manual(values = color.cell) +
  theme_lzy()
multiplot(p1, p2, p3, cols = 1)
dev.off()

aggregate(
  sample.info.9CL$Gene_Detected,
  by = list(Cell_Line = sample.info.9CL$Cell_Line),
  FUN = mean
)
aggregate(
  sample.info.9CL$Trans_Detected,
  by = list(Cell_Line = sample.info.9CL$Cell_Line),
  FUN = mean
)
aggregate(
  sample.info.9CL$Gene_Detected,
  by = list(Cell_Line = sample.info.9CL$Cell_Line),
  FUN = median
)
aggregate(
  sample.info.9CL$Trans_Detected,
  by = list(Cell_Line = sample.info.9CL$Cell_Line),
  FUN = median
)

aggregate(
  sample.info.4CL$Gene_Detected,
  by = list(Cell_Line = sample.info.4CL$Cell_Line),
  FUN = mean
)
aggregate(
  sample.info.4CL$Trans_Detected,
  by = list(Cell_Line = sample.info.4CL$Cell_Line),
  FUN = mean
)
aggregate(
  sample.info.4CL$Gene_Detected,
  by = list(Cell_Line = sample.info.4CL$Cell_Line),
  FUN = median
)
aggregate(
  sample.info.4CL$Trans_Detected,
  by = list(Cell_Line = sample.info.4CL$Cell_Line),
  FUN = median
)

## 1.2 accumulated frequency along expression level----------------
sensi.9CL <-
  myGeneSensitivity(rc_gene_9CL[, rownames(sample.info.9CL)])
sensi.4CL <-
  myGeneSensitivity(rc_gene_4CL[, rownames(sample.info.4CL)])
sensi.100 <-
  myGeneSensitivity(rc_gene_100[, rownames(sample.info.100)])
sensi.200 <-
  myGeneSensitivity(rc_gene_200[, rownames(sample.info.200)])
sensi.mES <- myGeneSensitivity(rc_mES)
sensi.9CL$library <- "9CL"
sensi.4CL$library <- "4CL"
sensi.100$library <- "UMI_100"
sensi.200$library <- "UMI_200"
sensi.mES$library = "SCAN-Seq"
sensi <- rbind(sensi.9CL, sensi.4CL, sensi.100, sensi.200, sensi.mES)
rm(sensi.9CL, sensi.4CL, sensi.100, sensi.200, sensi.mES)

pdf("Sensitivity/Percent_Detected.pdf", 5, 3)
ggplot(sensi, aes(x = breaks, y = frequency * 100, color = library)) +
  geom_line(size = 1) +
  geom_errorbar(aes(
    ymin = (frequency - sd) * 100,
    ymax = (frequency + sd) * 100
  ),
  width = .1,
  size = 0.5) +
  scale_color_manual(values = brewer.pal(5, "Set2")) +
  xlim(c(0, 3.5)) +
  xlab("Normalized Expression Level") +
  ylab("Gene Detection (%)") +
  theme_classic()
dev.off()

## 1.3 saturation  ----------------
sample.info.merge <- rbind(
  sample.info.100,
  sample.info.200,
  sample.info.4CL,
  sample.info.9CL,
  sample.info.9Mix
)
pdf("Sensitivity/LOESS.saturation.pdf", 4, 2.4)
ggplot(sample.info.merge) +
  #geom_point(aes(x = Mapped_Reads, y = Gene_Detected, color = "Gene"))+
  #geom_point(aes(x = Mapped_Reads, y = Trans_Detected, color = "Isoform"))+
  geom_smooth(
    aes(x = Mapped_Reads / 1e5, y = Gene_Detected, color = "Gene"),
    method = "loess",
    se = T,
    span = 0.4,
    size = 0.7
  ) +
  geom_smooth(
    aes(x = Mapped_Reads / 1e5, y = Trans_Detected, color = "Isoform"),
    method = "loess",
    se = T,
    span = 0.7,
    size = 0.7
  ) +
  xlab("Mapped Reads (10^5)") +
  theme_classic()
dev.off()

pdf("Sensitivity/LOESS.saturation.subset.pdf", 3, 2.4)
ggplot(sample.info.merge) +
  #geom_point(aes(x = Mapped_Reads, y = Gene_Detected, color = "Gene"))+
  #geom_point(aes(x = Mapped_Reads, y = Trans_Detected, color = "Isoform"))+
  geom_smooth(
    aes(x = Mapped_Reads / 1e5, y = Gene_Detected, color = "Gene"),
    method = "loess",
    se = T,
    span = 0.8,
    size = 0.7
  ) +
  geom_smooth(
    aes(x = Mapped_Reads / 1e5, y = Trans_Detected, color = "Isoform"),
    method = "loess",
    se = T,
    span = 0.8,
    size = 0.7
  ) +
  # scale_x_continuous(breaks = 0:4 * 5e4,
  #                    labels = c("0","5e4","1e5","3e4","4e4"),
  #                    limits = c(0,2e5))+
  xlim(c(0, 2)) +
  xlab("Mapped Reads (10^5)") +
  theme_classic()
dev.off()

# 3. cross contamination --------------------------------------------------
dir.create("CrossContamination")
rc_trans_9CL <- readRDS("../../batch/1 9CL/data/rc_trans.RDS")
rc_trans_4CL <- readRDS("../../batch/3 4CL/data/rc_trans.RDS")

mouse.trans <-
  rownames(rc_trans_4CL)[grep("ENSMUST", rownames(rc_trans_4CL))]
human.trans <-
  rownames(rc_trans_4CL)[-grep("ENSMUST", rownames(rc_trans_4CL))]
nrow(rc_trans_4CL) - length(mouse.trans) - length(human.trans)

Barnyard.4CL <- data.frame(
  row.names = colnames(rc_trans_4CL),
  "cell" = colnames(rc_trans_4CL),
  "human_transcripts" = colSums(rc_trans_4CL[human.trans, ]),
  "mouse_transcripts" = colSums(rc_trans_4CL[mouse.trans, ])
)
Barnyard.9CL <- data.frame(
  row.names = colnames(rc_trans_9CL),
  "cell" = colnames(rc_trans_9CL),
  "human_transcripts" = colSums(rc_trans_9CL[human.trans, ]),
  "mouse_transcripts" = colSums(rc_trans_9CL[mouse.trans, ])
)

Barnyard.4CL <- Barnyard.4CL[rownames(sample.info.4CL), ]
Barnyard.9CL <- Barnyard.9CL[rownames(sample.info.9CL), ]
Barnyard.4CL$Cell_Line <- sample.info.4CL$Cell_Line
Barnyard.9CL$Cell_Line <- sample.info.9CL$Cell_Line
Barnyard.4CL$Library <- "4CL"
Barnyard.9CL$Library <- "9CL"

Barnyard <- rbind(Barnyard.4CL, Barnyard.9CL)
Barnyard$ratio_human <-
  Barnyard$human_transcripts / (Barnyard$human_transcripts +
                                  Barnyard$mouse_transcripts)
Barnyard$ratio_mouse <-
  Barnyard$mouse_transcripts / (Barnyard$human_transcripts +
                                  Barnyard$mouse_transcripts)

Barnyard$ratio_human %>% hist()

Barnyard$Organism <- "Mixed"
Barnyard[Barnyard$ratio_human > 0.9, ]$Organism <- "Human"
Barnyard[Barnyard$ratio_mouse > 0.9, ]$Organism <- "Mouse"
table(Barnyard$Organism)
table(Barnyard$Organism) / nrow(Barnyard)

Barnyard.4CL <- Barnyard[Barnyard$Library == "4CL", ]
Barnyard.9CL <- Barnyard[Barnyard$Library == "9CL", ]

pdf("CrossContamination/Barnyard.4CL.pdf", 4, 3)
ggplot(Barnyard.4CL) +
  geom_point(aes(
    x = human_transcripts / 1000,
    y = mouse_transcripts / 1000,
    color = Organism
  ),
  size = 1) +
  theme_classic() +
  xlab("Human UMI counts (k)") +
  ylab("Mouse UMI counts (k)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal(1.7) +
  scale_color_manual(values = c("#d8756e", "#808183", "#6dbbc1")) +
  geom_abline(slope = 1 / 9, lty = 2) + geom_abline(slope = 9, lty = 2)
ggplot(Barnyard.4CL) +
  geom_point(aes(
    x = human_transcripts / 1000,
    y = mouse_transcripts / 1000,
    color = Cell_Line
  ),
  size = 0.3) +
  theme_classic() +
  scale_color_manual(values = color.cell) +
  xlab("Human Transcripts (10^3)") +
  ylab("Mouse Transcripts (10^3)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal(1.7) +
  geom_abline(slope = 1 / 9, lty = 2) + geom_abline(slope = 9, lty = 2)
dev.off()

pdf("CrossContamination/Barnyard.9CL.pdf", 4, 3)
ggplot(Barnyard.9CL) +
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
ggplot(Barnyard.9CL) +
  geom_point(aes(
    x = human_transcripts / 1000,
    y = mouse_transcripts / 1000,
    color = Cell_Line
  ),
  size = 0.3) +
  theme_bw() +
  scale_color_manual(values = color.cell) +
  xlab("Human Transcripts (10^3)") +
  ylab("Mouse Transcripts (10^3)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal(2) +
  geom_abline(slope = 1 / 9, lty = 2) + geom_abline(slope = 9, lty = 2)
dev.off()

table(Barnyard.4CL$Organism)
table(Barnyard.9CL$Organism)

n1 <-
  aggregate(Barnyard.4CL$human_transcripts,
            by = list(sample.info.4CL$Organism),
            sum)$x
n2 <-
  aggregate(Barnyard.4CL$mouse_transcripts,
            by = list(sample.info.4CL$Organism),
            sum)$x
n1 / (n1 + n2) # human 0.98089715
n2 / (n1 + n2) # mouse 0.96737273

n1 <-
  aggregate(Barnyard.9CL$human_transcripts,
            by = list(sample.info.9CL$Organism),
            sum)$x
n2 <-
  aggregate(Barnyard.9CL$mouse_transcripts,
            by = list(sample.info.9CL$Organism),
            sum)$x
n1 / (n1 + n2) # human 0.9925047
n2 / (n1 + n2) # mouse 0.974472796

# 4. Quantification -------------------------------------------------------
dir.create("Quantification")

## 4.1 correlation with NGS------------------------
rc_Smart3 <-
  read.table("../../Test/Smart-Seq3/Smartseq3.HEK.fwdprimer.readcounts.txt")
uc_Smart3 <-
  read.table("../../Test/Smart-Seq3/Smartseq3.HEK.fwdprimer.UMIcounts.txt")
sample.info.Smart3 <-
  read.table("../../Test/Smart-Seq3/Smartseq3.HEK.fwdprimer.sample_annotation.txt",
             header = T)
table(rownames(rc_Smart3) %in% gene.info.hg38$`Gene stable ID`)

rc_SCAN2.1 <-
  rc_gene_4CL[, rownames(sample.info.4CL)[sample.info.4CL$Cell_Line == "293T"]]
rc_SCAN2.2 <-
  rc_gene_9CL[, rownames(sample.info.9CL)[sample.info.9CL$Cell_Line == "293T"]]
dim(rc_SCAN2.1)
dim(rc_SCAN2.2)
rc_SCAN2 <- cbind(rc_SCAN2.1, rc_SCAN2.2)
rm(rc_SCAN2.1, rc_SCAN2.2)

common.gene <- intersect(rownames(rc_Smart3), rownames(uc_Smart3))
common.gene <-
  intersect(common.gene, gene.info.hg38$`Gene stable ID`)
names(common.gene) <- gene.info.hg38[common.gene, ]$`Gene name`
common.gene <- common.gene[!duplicated(names(common.gene))]
common.gene <-
  common.gene[names(common.gene) %in% rownames(rc_SCAN2)]
rc_Smart3 <- rc_Smart3[common.gene, ]
uc_Smart3 <- uc_Smart3[common.gene, ]
rownames(rc_Smart3) <- names(common.gene)
rownames(uc_Smart3) <- names(common.gene)
rc_SCAN2 <- rc_SCAN2[names(common.gene), ]

colnames(rc_Smart3) <- paste0(colnames(rc_Smart3), "_Read")
colnames(uc_Smart3) <- paste0(colnames(uc_Smart3), "_UMI")

# dot plot for correlation
plot.data <- data.frame(
  "Smart3_Reads" = rowSums(rc_Smart3),
  "Smart3_UMI" = rowSums(uc_Smart3),
  "SCANSeq2" = rowSums(rc_SCAN2)
)

cor(log10(plot.data$Smart3_Reads + 1),
    log10(plot.data$Smart3_UMI + 1)) # 0.9525138
cor(log10(plot.data$Smart3_Reads + 1),
    log10(plot.data$SCANSeq2 + 1)) # 0.8309494
cor(log10(plot.data$Smart3_UMI + 1), log10(plot.data$SCANSeq2 + 1)) # 0.8224188

model <-
  lm(log10(plot.data$SCANSeq2 + 1) ~ log10(plot.data$Smart3_Reads + 1))
summary(model)

pdf("Quantification/Dot.Cor.Methods.pdf", 5, 4)
ggplot(plot.data) +
  geom_point(aes(x = log10(Smart3_Reads + 1), y = log10(SCANSeq2 + 1)), size = 0.5) +
  geom_abline(slope = model$coefficients[2], color = "blue") +
  xlab("log10(Smart-Seq3 Reads + 1)") +
  ylab("log10(SCAN-Seq2 UMIs + 1)") +
  theme_classic() +
  coord_equal(1)
dev.off()
saveRDS(plot.data, "NGS_Correlation.rds")

# heatmap for correlation
rc_cmp <- cbind(rc_SCAN2, rc_Smart3, uc_Smart3)

anno <- data.frame(
  row.names = colnames(rc_cmp),
  "Group" = c(
    rep("SCANSeq2_4CL", 170),
    rep("SCANSeq2_9CL", 80),
    rep("SmartSeq3_Reads", 117),
    rep("SmartSeq3_UMI", 117)
  ),
  "temp" = 1
)
anno$temp <- NULL

identical(rownames(anno), colnames(rc_cmp))

cor.gene.cmp <- matrix(nrow = nrow(anno),
                       ncol = nrow(anno))
rownames(cor.gene.cmp) <- rownames(anno)
colnames(cor.gene.cmp) <- rownames(anno)

for (i in 1:nrow(anno)) {
  for (j in 1:nrow(anno)) {
    cor.gene.cmp[i, j] = cor(rc_cmp[, i], rc_cmp[, j])
  }
}
rm(rc_cmp, rc_SCAN2, rc_Smart3, uc_Smart3)

pdf("Quantification/Heatmap.cor.Methods.pdf", 7, 5)
pheatmap(
  cor.gene.cmp,
  main = "Compare with Smart-Seq3, HEK293T",
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = anno,
  annotation_col = anno,
  annotation_names_row = F,
  annotation_names_col = F,
  # annotation_colors = list(cell_line = mycolor),
  show_rownames = F,
  show_colnames = F
)
dev.off()

median(cor.gene.cmp[colnames(rc_SCAN2), colnames(rc_Smart3)]) # 0.7585603
median(cor.gene.cmp[colnames(rc_SCAN2), colnames(uc_Smart3)]) # 0.6603563
median(cor.gene.cmp[colnames(uc_Smart3), colnames(rc_Smart3)]) # 0.867181

## 4.2. UMI-ERCC correlation--------------------------
sample.info.merge <- rbind(
  sample.info.100,
  sample.info.200,
  sample.info.4CL,
  sample.info.9CL,
  sample.info.9Mix
)
sample.info.merge$Library <- factor(sample.info.merge$Library,
                                    levels = c("UMI_100", "UMI_200",
                                               "9CL", "9CL_Mix", "4CL"))
sample.info.merge <- rbind(
  sample.info.100,
  sample.info.200,
  sample.info.4CL,
  sample.info.9CL,
  sample.info.9Mix,
  sample.info.9Mix2[, colnames(sample.info.9Mix)],
  sample.info.9Mix3[, colnames(sample.info.9Mix)]
)
sample.info.merge$Library <- factor(
  sample.info.merge$Library,
  levels = c(
    "UMI_100",
    "UMI_200",
    "9CL",
    "9CL_Mix",
    "9CL_Mix2",
    "9CL_Mix3",
    "4CL"
  )
)

pdf("Quantification/Box.Cor.ERCC.UMI.2.pdf", 3.2, 3)
ggplot(sample.info.merge) +
  geom_boxplot(aes(x = Library, y = cor_ERCC_rc, fill = Library), outlier.size = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  )) +
  scale_fill_manual(values = brewer.pal(7, "Set2")) +
  ylab("Correlation Coefficient (r)")
dev.off()

## 4.3. Correlation among cells------------------------
# 9CL
anno <- data.frame(
  row.names = rownames(sample.info.9CL),
  "cell_line" = sample.info.9CL$Cell_Line,
  "temp" = 1
)
anno$temp <- NULL
rc_gene_qc <- rc_gene_9CL[, rownames(anno)]
identical(rownames(anno), colnames(rc_gene_qc))
cor.gene.9CL <- matrix(nrow = nrow(anno),
                       ncol = nrow(anno))
rownames(cor.gene.9CL) <- rownames(anno)
colnames(cor.gene.9CL) <- rownames(anno)
for (i in 1:nrow(anno)) {
  for (j in 1:nrow(anno)) {
    cor.gene.9CL[i, j] = cor(rc_gene_qc[, i], rc_gene_qc[, j])
  }
}

pdf("Quantification/Heatmap.cor.9CL.pdf", 7, 6)
pheatmap(
  cor.gene.9CL,
  main = "Correlation at Gene Level",
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = anno,
  annotation_col = anno,
  annotation_names_row = F,
  annotation_names_col = F,
  annotation_colors = list(cell_line = color.cell),
  show_rownames = F,
  show_colnames = F
)
dev.off()

cor.gene.9CL.melt <- data.frame()
for (one in unique(sample.info.9CL$Cell_Line)) {
  temp <- data.frame("cell_line" = one,
                     "Cor" = as.numeric(cor.gene.9CL[sample.info.9CL[sample.info.9CL$Cell_Line ==
                                                                       one, ]$Rename,
                                                     sample.info.9CL[sample.info.9CL$Cell_Line ==
                                                                       one, ]$Rename]))
  cor.gene.9CL.melt <- rbind(cor.gene.9CL.melt, temp)
}
table(cor.gene.9CL.melt$cell_line)
cor.gene.9CL.melt <- cor.gene.9CL.melt[cor.gene.9CL.melt$Cor < 1, ]

aggregate(cor.gene.9CL.melt$Cor,
          list(cor.gene.9CL.melt$cell_line),
          median)

pdf("Quantification/Density.Cor.9CL.pdf", 4, 3)
ggplot(cor.gene.9CL.melt) +
  geom_density(aes(x = Cor, color = cell_line), ) +
  scale_color_manual(values = color.cell) +
  xlab("Cell-to-cell correlation (r)") +
  theme_lzy()
dev.off()

# 4CL
anno <- data.frame(
  row.names = rownames(sample.info.4CL),
  "cell_line" = sample.info.4CL$Cell_Line,
  "temp" = 1
)
anno <- anno[order(anno$cell_line), ]
anno$temp <- NULL
rc_gene_qc <- rc_gene_4CL[, rownames(anno)]
identical(rownames(anno), colnames(rc_gene_qc))
cor.gene.4CL <- matrix(nrow = nrow(anno),
                       ncol = nrow(anno))
rownames(cor.gene.4CL) <- rownames(anno)
colnames(cor.gene.4CL) <- rownames(anno)
for (i in 1:nrow(anno)) {
  for (j in 1:nrow(anno)) {
    cor.gene.4CL[i, j] = cor(rc_gene_qc[, i], rc_gene_qc[, j])
  }
}

pdf("Quantification/Heatmap.cor.4CL.pdf", 7, 6)
pheatmap(
  cor.gene.4CL,
  main = "Correlation at Gene Level",
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = anno,
  annotation_col = anno,
  annotation_names_row = F,
  annotation_names_col = F,
  annotation_colors = list(cell_line = color.cell),
  show_rownames = F,
  show_colnames = F
)
dev.off()

cor.gene.4CL.melt <- data.frame()
for (one in unique(sample.info.4CL$Cell_Line)) {
  temp <- data.frame("cell_line" = one,
                     "Cor" = as.numeric(cor.gene.4CL[sample.info.4CL[sample.info.4CL$Cell_Line ==
                                                                       one, ]$Rename,
                                                     sample.info.4CL[sample.info.4CL$Cell_Line ==
                                                                       one, ]$Rename]))
  cor.gene.4CL.melt <- rbind(cor.gene.4CL.melt, temp)
}
table(cor.gene.4CL.melt$cell_line)
cor.gene.4CL.melt <- cor.gene.4CL.melt[cor.gene.4CL.melt$Cor < 1, ]

aggregate(cor.gene.4CL.melt$Cor,
          list(cor.gene.4CL.melt$cell_line),
          median)

pdf("Quantification/Density.Cor.4CL.pdf", 4, 3)
ggplot(cor.gene.4CL.melt) +
  geom_density(aes(x = Cor, color = cell_line), ) +
  scale_color_manual(values = color.cell) +
  xlab("Cell-to-cell correlation (r)") +
  theme_lzy()
dev.off()

# 293T Smart-Seq3 rc
rc_Smart3 <-
  read.table("../../Test/Smart-Seq3/Smartseq3.HEK.fwdprimer.readcounts.txt")
cor.gene.Smart3 <- matrix(nrow = ncol(rc_Smart3),
                          ncol = ncol(rc_Smart3))
rownames(cor.gene.Smart3) <- colnames(rc_Smart3)
colnames(cor.gene.Smart3) <- colnames(rc_Smart3)
for (i in 1:ncol(rc_Smart3)) {
  for (j in 1:ncol(rc_Smart3)) {
    cor.gene.Smart3[i, j] = cor(rc_Smart3[, i], rc_Smart3[, j])
  }
}
cor.gene.Smart3.melt <- data.frame("cell_line" = "293T",
                                   "Cor" = as.numeric(cor.gene.Smart3))
table(cor.gene.Smart3.melt$cell_line)
cor.gene.Smart3.melt <-
  cor.gene.Smart3.melt[cor.gene.Smart3.melt$Cor < 1, ]
aggregate(cor.gene.Smart3.melt$Cor,
          list(cor.gene.Smart3.melt$cell_line),
          median)
# 293T Smart-Seq3 umi
uc_Smart3 <-
  read.table("../../Test/Smart-Seq3/Smartseq3.HEK.fwdprimer.UMIcounts.txt")
cor.gene.Smart3.umi <- matrix(nrow = ncol(uc_Smart3),
                              ncol = ncol(uc_Smart3))
rownames(cor.gene.Smart3.umi) <- colnames(uc_Smart3)
colnames(cor.gene.Smart3.umi) <- colnames(uc_Smart3)
for (i in 1:ncol(uc_Smart3)) {
  for (j in 1:ncol(uc_Smart3)) {
    cor.gene.Smart3.umi[i, j] = cor(uc_Smart3[, i], uc_Smart3[, j])
  }
}
cor.gene.Smart3.umi.melt <- data.frame("cell_line" = "293T",
                                       "Cor" = as.numeric(cor.gene.Smart3.umi))
table(cor.gene.Smart3.umi.melt$cell_line)
cor.gene.Smart3.umi.melt <-
  cor.gene.Smart3.umi.melt[cor.gene.Smart3.umi.melt$Cor < 1, ]
aggregate(cor.gene.Smart3.umi.melt$Cor,
          list(cor.gene.Smart3.umi.melt$cell_line),
          median)


# compare
cor.gene.9CL.293T <-
  cor.gene.9CL.melt[cor.gene.9CL.melt$cell_line == "293T", ]
cor.gene.9CL.293T$Method <- "SCAN-Seq2"
cor.gene.Smart3.melt$Method <- "Smart-Seq3_Reads"
cor.gene.Smart3.umi.melt$Method <- "Smart-Seq3_UMIs"
plot.data <-
  rbind(cor.gene.9CL.293T,
        cor.gene.Smart3.melt,
        cor.gene.Smart3.umi.melt)

pdf("Quantification/Vln.Cell2Cell.cor.pdf", 4, 3.3)
ggplot(plot.data) +
  geom_violin(aes(x = Method, y = Cor, fill = Method)) +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 60,
    hjust = 1,
    vjust = 1
  )) +
  scale_fill_manual(values = brewer.pal(4, "Set2")) +
  ylab("Cell-to-cell correlation (r)")
dev.off()

# 5. distinguish isoform --------------------------------------------------
dir.create("IsoformDistinction")
rc_trans_9CL <- readRDS("../../batch/1 9CL/data/rc_trans.RDS")
rc_trans_9Mix <- readRDS("../../batch/2 9CL_Mix/data/rc_trans.RDS")

obj.9CL <-
  readRDS("E:/LabWork/Project/SCANSeq2/batch/1 9CL/obj.9CL.rds")
obj.9Mix <-
  readRDS("E:/LabWork/Project/SCANSeq2/batch/2 9CL_Mix/obj.9CL_Mix.rds")
obj.9CL[["Isoform"]] <-
  CreateAssayObject(rc_trans_9CL[, colnames(obj.9CL)])
obj.9Mix[["Isoform"]] <-
  CreateAssayObject(rc_trans_9Mix[, colnames(obj.9Mix)])

obj.combined <- merge(
  obj.9CL,
  y = obj.9Mix,
  add.cell.ids = c("9CL", "9MIX"),
  project = "SCAN_Seq"
)
table(Idents(obj.combined))
rm(obj.9CL, obj.9Mix, rc_trans_9CL, rc_trans_9Mix)

DefaultAssay(obj.combined) <- "Isoform"
obj.combined <- NormalizeData(obj.combined)
obj.combined <- ScaleData(obj.combined,
                          features =  protien.info.hg38[protien.info.hg38$`Gene name` %in%
                                                          obj.combined@assays$Gene@var.features, ]$`Transcript name`)

PTPRC.iso <-
  grep("PTPRC-", rownames(obj.combined), value = T) # 17 isoform
PTPRC.iso <-
  PTPRC.iso[trans.info.hg38[PTPRC.iso, ]$`Transcript type` == "protein_coding"]
length(PTPRC.iso)  # 10 coding

PTPRC.iso <- sort(PTPRC.iso)

pdf("IsoformDistinction/Vln.iso.PTPRC.all.pdf", 8, 8.5)
VlnPlot(
  obj.combined,
  c("PTPRC", PTPRC.iso),
  idents = c("GM12878", "H9", "K562"),
  cols = color.cell,
  ncol = 3,
  same.y.lims = T
)
dev.off()

PTPRC.iso.selected <-
  PTPRC.iso[!PTPRC.iso %in% c("PTPRC-215", "PTPRC-216")]

pdf("IsoformDistinction/Vln.iso.PTPRC.selected.pdf", 8, 7)
VlnPlot(
  obj.combined,
  PTPRC.iso.selected,
  idents = c("GM12878", "H9", "K562"),
  cols = color.cell,
  ncol = 3,
  same.y.lims = T
)
dev.off()
pdf("IsoformDistinction/Vln.iso.PTPRC.selected.9CL.pdf", 7.5, 7)
VlnPlot(
  obj.9CL,
  c("PTPRC", PTPRC.iso.selected),
  idents = c("GM12878", "H9", "K562"),
  cols = color.cell,
  ncol = 3,
  same.y.lims = T
)
dev.off()
rm(obj.combined)
rm(
  rc_gene_qc,
  rc_trans_4CL,
  rc_trans_9CL,
  rc_gene_100,
  rc_gene_200,
  rc_gene_4CL,
  rc_gene_9CL
)
rm(rc_mES, rc_Smart3, temp, uc_Smart3)

# plot CD45 trancript
library(magrittr)
library(ggtranscript)
library(ggplot2)

gtf <-
  rtracklayer::import("IsoformDistinction/Homo_sapiens.GRCh38.101.PTPRC.gtf")
gtf = as.data.frame(gtf)

PTPRC_exons <- gtf %>% dplyr::filter(type == "exon")
temp <-
  PTPRC_exons[PTPRC_exons$transcript_biotype == "protein_coding", ]
temp <- temp[!duplicated(temp$exon_id), ]
temp$transcript_name <- "gene_PTPRC"
PTPRC_exons <- rbind(PTPRC_exons, temp)
PTPRC_exons.select <-
  PTPRC_exons[PTPRC_exons$transcript_name %in% c(c("gene_PTPRC",
                                                   "PTPRC-201",
                                                   "PTPRC-203",
                                                   "PTPRC-209")), ]
PTPRC_exons.select$transcript_name <-
  factor(PTPRC_exons.select$transcript_name,
         levels = rev(c(
           "gene_PTPRC",
           "PTPRC-201",
           "PTPRC-203",
           "PTPRC-209"
         )))

pdf("IsoformDistinction/CD45.Gene.pdf", 6, 1)
temp %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
  geom_range(aes(fill = transcript_biotype),
             show.legend = F,
             size = 0) +
  geom_intron(data = to_intron(temp, "transcript_name"), arrow = NULL) +
  scale_fill_manual(values = "#0000f7") +
  geom_vline(xintercept = c(198692700, 198703000)) +
  theme_classic()
dev.off()

pdf("IsoformDistinction/CD45.Isoform.pdf", 6, 4)
PTPRC_exons.select %>%
  ggplot(aes(xstart = start, xend = end, y = transcript_name)) +
  geom_range(aes(fill = transcript_biotype)) +
  geom_intron(data = to_intron(PTPRC_exons.select, "transcript_name"),
              aes(strand = strand)) +
  coord_cartesian(xlim = c(198692700, 198703000)) +
  scale_fill_manual(values = "#0000f7") +
  theme_classic()
dev.off()


# 6. Supplementary tables -------------------------------------------------
sample.info.merge <- rbind(
  sample.info.9CL,
  sample.info.4CL,
  sample.info.9Mix,
  sample.info.100,
  sample.info.200
)
table(sample.info.merge$Library)
write.csv(sample.info.merge, "sample.info.merge.csv")
