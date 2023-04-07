#!/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/miniconda3/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

sample.list <- read.table(args[1])
sample.list <- sample.list$V1

cell1 <- sample.list[1]

gene.quant <- read.table(paste("02.process/", cell1, "/Salmon_Quant/", cell1, "_gene_quant.txt", sep = ""),
  head = T, stringsAsFactors = F
)
trans.quant <- read.table(paste("02.process/", cell1, "/Salmon_Quant/", cell1, "_trans_quant.txt", sep = ""),
  head = T, stringsAsFactors = F
)

gene.count <- gene.quant[, c(2, 1)]
trans.count <- trans.quant[, c(2, 1)]

rm(gene.quant, trans.quant, cell1)

for (cell in sample.list) {
  gene.quant <- read.table(paste("02.process/", cell, "/Salmon_Quant/", cell, "_gene_quant.txt", sep = ""),
    head = T, stringsAsFactors = F, row.names = 1
  )
  trans.quant <- read.table(paste("02.process/", cell, "/Salmon_Quant/", cell, "_trans_quant.txt", sep = ""),
    head = T, stringsAsFactors = F, row.names = 1
  )
  gene.count[, paste(cell)] <- gene.quant$NumReads
  trans.count[, paste(cell)] <- trans.quant$NumReads
}

gene.count <- gene.count[, 2:ncol(gene.count)]
trans.count <- trans.count[, 2:ncol(trans.count)]


write.table(gene.count, "03.summary/Gene.Count.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(trans.count, "03.summary/Trans.Count.txt", quote = F, row.names = F, col.names = T, sep = "\t")
