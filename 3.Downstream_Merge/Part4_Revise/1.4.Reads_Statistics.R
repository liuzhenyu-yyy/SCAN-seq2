library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

setwd("E:/LabWork/Project/SCANSeq2/revise/1.4.Reads_Statistics/")
source("E:/LabWork/code/MyFunction.R")


# 1. load data ------------------------------------------------------------
obj.9CL <- readRDS("../../batch/1 9CL/obj.9CL.rds")
sample.info <- obj.9CL@meta.data
rownames(sample.info) <- sample.info$Name

cell.list <- dir("data/") %>% gsub("_stat.summary.txt", "", .)
reads.stat <- data.frame()

for (cell in rownames(sample.info)) {
  counts <-
    read.table(paste0("data/", cell, "_stat.summary.txt"), row.names = 1) %>%
    t %>% as.data.frame
  counts$Cell <- cell
  counts$Stage <-
    c("Demultiplexing", "Quality_control", "Deduplication")
  reads.stat <- rbind(reads.stat, counts)
}

colnames(reads.stat)

reads.stat$log10_Total_bases <- log10(reads.stat$Total_bases)
reads.stat$log10_Number_of_reads <- log10(reads.stat$Number_of_reads)
reads.stat$Stage <- factor(reads.stat$Stage,
                           levels = c("Demultiplexing", "Quality_control", "Deduplication"))

# 2. Visualization --------------------------------------------------------
color.stage <-  brewer.pal(4, "Set3")
names(color.stage) <-
  c("Total", "Demultiplexing", "Quality_control", "Deduplication")

## 2.1. Boxplot of single-cell statistics -----

colnames(reads.stat)
plot.data <-
  melt(reads.stat[reads.stat$Cell %in% rownames(sample.info), c(1:4, 6, 8:11)],
       id.vars = c("Cell", "Stage"))
plot.data$Stage <- factor(plot.data$Stage,
                          levels = c("Demultiplexing", "Quality_control", "Deduplication"))
plot.data$variable <- factor(
  plot.data$variable,
  levels = c(
    "log10_Number_of_reads",
    "log10_Total_bases",
    "Read_length_N50",
    "Mean_read_length",
    "Median_read_length",
    "Mean_read_quality" ,
    "Median_read_quality"
  )
)

pdf("Box.Reads.Statitics.pdf", 8, 5)
ggplot(plot.data, aes(x = Stage, y = value)) +
  geom_boxplot(aes(fill = Stage), alpha = 1, outlier.size = 0.7) +
  facet_wrap( ~ variable, scales = "free", ncol = 3) +
  theme_bw() + scale_fill_manual(values = color.stage) +
  theme(axis.text.x = element_blank())
dev.off()

aggregate(reads.stat$Number_of_reads, by = list(reads.stat$Stage), mean)
aggregate(reads.stat$Number_of_reads, by = list(reads.stat$Stage), median)

## 2.1. circle plot of overall -----
names(color.stage) <-
  c("PromethION_output",
    "Demultiplexing",
    "Quality_control",
    "Deduplication")

plot.data <-
  data.frame(
    "Stage" = c(
      "PromethION_output",
      "Demultiplexing",
      "Quality_control",
      "Deduplication"
    )
  )
plot.data$Number_of_Reads <-
  c(84801426, aggregate(reads.stat$Number_of_reads,
                        list(reads.stat$Stage),
                        sum)$x)
plot.data$Total_bases <-
  c(106795866640, aggregate(reads.stat$Total_bases,
                            list(reads.stat$Stage),
                            sum)$x)
test <- data.frame(
  "Stage" = paste("test", 1:3),
  "Number_of_Reads" = NA,
  "Total_bases" = NA
)
plot.data <- rbind(plot.data, test)

plot.data$Stage <-
  factor(plot.data$Stage, levels = rev(
    c(
      "PromethION_output",
      "Demultiplexing",
      "Quality_control",
      "Deduplication",
      test$Stage
    )
  ))

pdf("Circle.Reads.Bases.pdf", 8, 6)
multiplot(
  ggplot(plot.data, aes(
    x = Stage, y = Number_of_Reads / 1e6, fill = Stage
  )) +
    geom_bar(stat = "identity") +
    coord_polar(
      theta = "y",
      start = 0,
      direction = 1
    ) +
    ylim(c(0, max(
      plot.data$Number_of_Reads
    ))) +
    scale_fill_manual(values = color.stage) +
    ggtitle("Number of reads (million)") +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(breaks = seq(0, 100, 10)),
  ggplot(plot.data, aes(
    x = Stage, y = Total_bases / 1e9, fill = Stage
  )) +
    geom_bar(stat = "identity") +
    coord_polar(
      theta = "y",
      start = 0,
      direction = 1
    ) +
    ylim(c(0, max(
      plot.data$Total_bases
    ))) +
    scale_fill_manual(values = color.stage) +
    ggtitle("Total bases (billion)") +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(breaks = seq(0, 1200, 10)),
  cols = 1
)
dev.off()
