library(Biostrings)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(corrplot)

source("E:/LabWork/code/MyFunction.R")


# 1. Editing distance of 96 barcodes --------------------------------------

barcode <-
  readDNAStringSet(
    "../1.1.Error_Rate/data/Barcode01_96.fa",
    format = "fasta",
    seek.first.rec = FALSE,
    use.names = T
  )

adist(barcode$Barcode01,
      barcode$Barcode02,
      partial = F,
      ignore.case = T)

qdis.mat <- matrix(0, 96, 96)
rownames(dis.mat) <- names(barcode)
colnames(dis.mat) <- names(barcode)

for (i in names(barcode)) {
  for (j in names(barcode)) {
    dis.mat[i, j] <- adist(barcode[[i]], barcode[[j]])
  }
}

pdf("Heatmap.96barcode.pdf")
pheatmap(dis.mat)
dev.off()

min(dis.mat[dis.mat != 0])

# 2. expand barcode set ---------------------------------------------------
n_ideal <- 999
cutoff.dist <- 11
barcode_plus <- barcode
id <- 97
times <- 0

while (length(barcode_plus) < n_ideal) {
  times <- times + 1

  seq <- DNA_ALPHABET[1:4] %>%
    sample(size = 24, replace = TRUE) %>%
    paste(collapse = "") %>%
    DNAString()

  min_dist <-
    lapply(barcode_plus, function(x) {
      return(as.numeric(adist(x, seq)))
    }) %>%
    unlist %>% min

  if (min_dist >= cutoff.dist) {
    barcode_plus[[paste("Barcode", id, sep = "")]] <- seq
    id <- id + 1
    message(
      paste0(
        Sys.time(),
        ": ",
        id,
        " barcodes identified, ",
        times,
        " barcodes tried in this round."
      )
    )
    times <- 0
  }
}
