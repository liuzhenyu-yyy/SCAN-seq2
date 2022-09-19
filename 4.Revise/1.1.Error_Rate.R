library(Biostrings)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("E:/LabWork/Project/SCANSeq2/revise/1.1.Error_Rate/")
source("E:/LabWork/code/MyFunction.R")

# 1. Barcode error rate ---------------------------------------------------

## 1.1. load data ----
obj.9CL <- readRDS("../../batch/1 9CL/obj.9CL.rds")
sample.info <- obj.9CL@meta.data
rownames(sample.info) <- sample.info$Name
sample.info$Barcode.3 <- paste0("Barcode", sample.info$Barcode.3)
sample.info$Barcode.5 <- paste0("Barcode", sample.info$Barcode.5)

barcode <- readDNAStringSet("data/Barcode01_96.fa", format="fasta",
                               seek.first.rec=FALSE, use.names=T)

## 1.2. Calculate average for each single cell ----
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)

sample.info$Mean_Mismatch_5p <- NA
sample.info$Mean_Mismatch_3p <- NA
sample.info$Mean_Indel_5p <- NA
sample.info$Mean_Indel_3p <- NA

for (cell in rownames(sample.info)) {
  message(paste0(Sys.time(),": Calculating barcode mismatch in cell ", cell, "..."))
  
  reads <- readDNAStringSet(paste0("data/subset.fq/",cell,"_sub5000_clean.fq.gz"),
                            format = "fastq", seek.first.rec=FALSE, use.names=F)
  
  # initiate statistics
  nMismatch_5p <- c()
  nMismatch_3p <- c()
  nIndel_5p <- c()
  nIndel_3p <- c()
  
  for (one in 1:1000) {
    agrep.5p <- adist(barcode[[sample.info[cell,]$Barcode.5]],
                      reads[[one]],
                      partial = T,
                      ignore.case = T)
    agrep.3p <- adist(barcode[[sample.info[cell,]$Barcode.3]],
                      partial = T,
                      reads[[one]],counts = T,
                      ignore.case = T)
    
    align.5p <- pairwiseAlignment(barcode[[sample.info[cell,]$Barcode.5]], reads[[one]],
                                  substitutionMatrix = sigma,
                                  gapOpening = 4,
                                  gapExtension = 1,
                                  scoreOnly = FALSE,
                                  type="global-local")
    align.3p <- pairwiseAlignment(barcode[[sample.info[cell,]$Barcode.3]], reads[[one]],
                                  substitutionMatrix = sigma,
                                  gapOpening = 4,
                                  gapExtension = 1,
                                  scoreOnly = FALSE,
                                  type="global-local")
    
    if (align.5p@score > 10) {
      score_5p <- c(score_5p, align.5p@score)
      nMismatch_5p <- c(nMismatch_5p, nmismatch(align.5p))
      nIndel_5p <- c(nIndel_5p, nindel(align.5p)@insertion[1,1] +  nindel(align.5p)@deletion[1,1])
    }
    
    if (align.3p@score > 10) {
      score_3p <- c(score_3p, align.3p@score)
      nMismatch_3p <- c(nMismatch_3p, nmismatch(align.3p))
      nIndel_3p <- c(nIndel_3p, nindel(align.3p)@insertion[1,1] +  nindel(align.3p)@deletion[1,1])
    }
  }
  
  sample.info[cell,]$Mean_Mismatch_5p <- mean(nMismatch_5p)
  sample.info[cell,]$Mean_Mismatch_3p <- mean(nMismatch_3p)
  sample.info[cell,]$Mean_Indel_5p <- mean(nIndel_5p)
  sample.info[cell,]$Mean_Indel_3p <- mean(nIndel_3p)
}

for (cell in rownames(sample.info)) {
  message(paste0(Sys.time(),": Calculating barcode mismatch in cell ", cell, "..."))
  
  reads <- readDNAStringSet(paste0("data/subset.fq/",cell,"_sub5000_clean.fq.gz"),
                            format = "fastq", seek.first.rec=FALSE, use.names=F)
  
  # initiate statistics
  score_5p <- c()
  score_3p <- c()
  nMismatch_5p <- c()
  nMismatch_3p <- c()
  nIndel_5p <- c()
  nIndel_3p <- c()
  
  for (one in 1:5000) {
    align.5p <- pairwiseAlignment(barcode[[sample.info[cell,]$Barcode.5]], reads[[one]],
                                  substitutionMatrix = sigma,
                                  gapOpening = 4,
                                  gapExtension = 1,
                                  scoreOnly = FALSE,
                                  type="global-local")
    align.3p <- pairwiseAlignment(barcode[[sample.info[cell,]$Barcode.3]], reads[[one]],
                                  substitutionMatrix = sigma,
                                  gapOpening = 4,
                                  gapExtension = 1,
                                  scoreOnly = FALSE,
                                  type="global-local")
    
    if (align.5p@score > 10) {
      score_5p <- c(score_5p, align.5p@score)
      nMismatch_5p <- c(nMismatch_5p, nmismatch(align.5p))
      nIndel_5p <- c(nIndel_5p, nindel(align.5p)@insertion[1,1] +  nindel(align.5p)@deletion[1,1])
    }
    
    if (align.3p@score > 10) {
      score_3p <- c(score_3p, align.3p@score)
      nMismatch_3p <- c(nMismatch_3p, nmismatch(align.3p))
      nIndel_3p <- c(nIndel_3p, nindel(align.3p)@insertion[1,1] +  nindel(align.3p)@deletion[1,1])
    }
  }
  
  sample.info[cell,]$Mean_Mismatch_5p <- mean(nMismatch_5p)
  sample.info[cell,]$Mean_Mismatch_3p <- mean(nMismatch_3p)
  sample.info[cell,]$Mean_Indel_5p <- mean(nIndel_5p)
  sample.info[cell,]$Mean_Indel_3p <- mean(nIndel_3p)
}

rm(nMismatch_5p, nMismatch_3p, nIndel_5p, nIndel_3p, align.5p, align.3p, reads)

mean(sample.info$Mean_Mismatch_5p)
mean(sample.info$Mean_Indel_5p)

mean(sample.info$Mean_Mismatch_5p)+
  mean(sample.info$Mean_Indel_5p)

## 1.3. Statistics of barcode error ----
pdf("Box.Barcode.Error.3p.pdf",8,6)
multiplot(
  ggplot(sample.info, aes(x = gsub("Barcode","",Barcode.3), y = Mean_Mismatch_3p))+
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~"Mismatch") + theme_lzy() +
    theme(axis.title.x = element_blank()) + ylab("nMismatch"),
  ggplot(sample.info, aes(x = gsub("Barcode","",Barcode.3), y = Mean_Indel_3p))+
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~"Indel") + theme_lzy()+
    theme(axis.title.x = element_blank()) + ylab("nIndel"),
  ggplot(sample.info, aes(x = gsub("Barcode","",Barcode.3), y = Mean_Mismatch_3p + Mean_Indel_3p))+
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~"Overall Error")+ theme_lzy()+
    xlab("3 prime barcode") + ylab("Total Error"),
  cols = 1)
dev.off()
pdf("Box.Barcode.Error.5p.pdf",8,6)
multiplot(
  ggplot(sample.info, aes(x = gsub("Barcode","",Barcode.5), y = Mean_Mismatch_5p))+
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~"Mismatch") + theme_lzy() +
    theme(axis.title.x = element_blank()) + ylab("nMismatch"),
  ggplot(sample.info, aes(x = gsub("Barcode","",Barcode.5), y = Mean_Indel_5p))+
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~"Indel") + theme_lzy()+
    theme(axis.title.x = element_blank()) + ylab("nIndel"),
  ggplot(sample.info, aes(x = gsub("Barcode","",Barcode.5), y = Mean_Mismatch_5p + Mean_Indel_5p))+
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~"Overall Error")+ theme_lzy()+
    xlab("5 prime barcode") + ylab("Total Error"),
  cols = 1)
dev.off()


# 2. UMI NGS TGS ----------------------------------------------------------
## 2.1. load data -----

gene2trans <- read.table("data/TranscriptID_GeneID_from_cDNA_fasta.txt")
colnames(gene2trans) <- c("Trans_id","Gene_id")
rownames(gene2trans) <- gene2trans$Trans_id

UMI_NGS <- as.data.frame(data.table::fread("data/NGS-UMI-NGS-UMI-1.NGS.UMIs.txt",
                                           header = F))
UMI_NGS$Method <- "NGS"
UMI_TGS <- as.data.frame(data.table::fread("data/M_4CL.TGS.UMIs.txt",
                                           header = F))
UMI_TGS$Method <- "TGS"

colnames(UMI_NGS)[1] <- "ID"
UMI_NGS$Trans_id <- gsub("_.+?$","",UMI_NGS$ID)
UMI_NGS$Seq <- gsub("^.+?_","",UMI_NGS$ID)
UMI_NGS$Gene_id <- gene2trans[UMI_NGS$Trans_id,]$Gene_id
UMI_NGS$ID <- paste(UMI_NGS$Gene_id, UMI_NGS$Seq, sep = "_")
UMI_NGS <- UMI_NGS[!duplicated(UMI_NGS$ID),] # 7521407 of 10318388

colnames(UMI_TGS)[1] <- "ID"
UMI_TGS$Trans_id <- gsub("_.+?$","",UMI_TGS$ID)
UMI_TGS$Seq <- gsub("^.+?_","",UMI_TGS$ID)
UMI_TGS$Gene_id <- gene2trans[UMI_TGS$Trans_id,]$Gene_id
UMI_TGS$ID <- paste(UMI_TGS$Gene_id, UMI_TGS$Seq, sep = "_")
UMI_TGS <- UMI_TGS[!duplicated(UMI_TGS$ID),] # 26479443 of 27434757

table(UMI_TGS$ID %in% UMI_NGS$ID) #F:24992853  T:1486590

## 2.2. ERCC no error ------
library(Vennerable)
temp <- table(UMI_NGS$Gene_id)
temp <- sort(temp, decreasing = T)
head(temp)

ERCC <- grep("ERCC",UMI_NGS$Gene_id, value = T)
ERCC <- unique(ERCC)

UMI_TGS.ERCC <- UMI_TGS[UMI_TGS$Gene_id %in% ERCC,]
UMI_NGS.ERCC <- UMI_NGS[UMI_NGS$Gene_id %in% ERCC,]
rownames(UMI_TGS.ERCC) <- UMI_TGS.ERCC$ID
rownames(UMI_NGS.ERCC) <- UMI_NGS.ERCC$ID

x <- UMI_TGS.ERCC$ID
y <- UMI_NGS.ERCC$ID

pdf("Venn.ERCC_UMI.E0.pdf", 7, 7)
plot(Venn(SetNames = c("Illumina","ONT"),
          Weight =  c(`01` = sum(!x %in% y),`11` = sum(x %in% y), `10` = sum(!y %in% x))),
          doWeights = TRUE,
          show = list( SetLabels = T, Faces = FALSE))
dev.off()
128117 / length(x) # 0.2822706

## 2.3. allow edit distance 1 -----
UMI_TGS.ERCC$Min_Dist <- 100000
UMI_TGS.ERCC$Min_UMI <- "none"

sort(table(UMI_TGS.ERCC$Gene_id))
sort(table(UMI_NGS.ERCC$Gene_id))

i = 1
for (one in ERCC) {
  message(paste0(Sys.time(),": processing ", one,", the ",i," th of 78 ..."))
  ID_NGS <- UMI_NGS.ERCC[UMI_NGS.ERCC$Gene_id == one,]$ID
  seq_NGS <- gsub("^.+?_","",ID_NGS)
  ID_TGS <- UMI_TGS.ERCC[UMI_TGS.ERCC$Gene_id == one,]$ID
  seq_TGS <- gsub("^.+?_","",ID_TGS)
  
  
  seq_TGS_split <- split(seq_TGS, ceiling(seq_along(seq_TGS)/10000))
  min.dist <- c()
  min.umi <- c()
  message(paste0(paste0(Sys.time(), ": ",length(seq_TGS)," UMIs detected by TGS, spliting into ",
                        length(seq_TGS_split), " trunks.")))
          
  for (trunk in names(seq_TGS_split)) {
    message(paste0(paste0(Sys.time(), ": Processing trunk ", trunk,"...")))
    dis.mat <- adist(seq_TGS_split[[trunk]], seq_NGS)
    min.dist_ <- apply(dis.mat, 1, min)
    min.umi_ <- apply(dis.mat, 1, which.min)
    min.dist <- c(min.dist, min.dist_)
    min.umi <- c(min.umi, min.umi_)
    gc()
  }
  
  UMI_TGS.ERCC[ID_TGS,]$Min_Dist <- min.dist
  UMI_TGS.ERCC[ID_TGS,]$Min_UMI <- seq_NGS[min.umi]
  
  i <- i + 1
}

rm(dis.mat, min.dist_, min.umi_, min.dist, min.umi,
   ID_NGS, ID_TGS, seq_NGS, seq_TGS, seq_TGS_split)

# overlap
UMI_TGS.ERCC$ID_corrected <- paste(UMI_TGS.ERCC$Gene_id, UMI_TGS.ERCC$Min_UMI, sep = "_")
selcted <- which(UMI_TGS.ERCC$Min_Dist > 1)
UMI_TGS.ERCC[selcted,]$ID_corrected <- UMI_TGS.ERCC[selcted,]$ID

x <- UMI_TGS.ERCC$ID_corrected
y <- UMI_NGS.ERCC$ID

pdf("Venn.ERCC_UMI.E1.pdf", 7, 7)
plot(Venn(SetNames = c("Illumina","ONT"),
          Weight =  c(`01` = sum(!x %in% y),`11` = sum(x %in% y), `10` = sum(!y %in% x))),
     doWeights = TRUE,
     show = list( SetLabels = T, Faces = FALSE))
dev.off()
