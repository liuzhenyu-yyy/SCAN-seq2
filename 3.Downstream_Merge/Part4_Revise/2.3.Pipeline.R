setwd("E:/LabWork/Project/SCANSeq2/revise/2.3.Pipeline/")

library(ggplot2)


load("../1.4.Reads_Statistics/1.4.Reads_Statistics.RData")

aggregate(reads.stat$Number_of_reads, list(reads.stat$Stage),"mean")

mean(sample.info$Mapped_Reads)

Isoform.raw <- read.table("Isoform.category.raw.tsv", sep = "\t", row.names = 1)
Isoform.qc <- read.table("Isoform.category.qc.tsv", sep = "\t", row.names = 1)
Isoform.raw <- Isoform.raw[sample.info$Name,]
Isoform.qc <- Isoform.qc[sample.info$Name,]
colnames(Isoform.raw) <- c("Total_Gene","Total_Trans","FSM","ISM","NIC","NNC","Genic_Genomic","Antisense","Fusion","Intergenic","Genic_intron","NIC_CJ","NIC_CS","NIC_IR" )
colnames(Isoform.qc) <- c("Total_Gene","Total_Trans","FSM","ISM","NIC","NNC","Genic_Genomic","Antisense","Fusion","Intergenic","Genic_intron","NIC_CJ","NIC_CS","NIC_IR" )

mean(Isoform.raw$Total_Gene)
mean(Isoform.raw$Total_Trans)
mean(Isoform.qc$Total_Gene)
mean(Isoform.qc$Total_Trans)
