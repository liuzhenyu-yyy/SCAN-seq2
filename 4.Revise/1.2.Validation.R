setwd("E:/LabWork/Project/SCANSeq2/revise/1.2")
library(dplyr)
library(Seurat)
library(DTUrtle)
library(ggplot2)
library(Vennerable)

source("E:/LabWork/code/MyFunction.R")
ReadGeneInfo("hg38", classify.protein = T)


# 1. select novel isoform -------------------------------------------------
## GM12878
GM12878.bed <- read.table("data/GM12878.c10.bed")
classification <- read.table("data/GM12878.c10.SQANTI3_classification.txt", header = T)
report <- read.table("data/GM12878_trans_report.txt", header = T)

rownames(classification) <- classification$isoform
rownames(GM12878.bed) <- gsub("^.+;","",GM12878.bed$V4)
rownames(report) <- report$transcript_id

NIC <- classification$isoform[classification$structural_category == "novel_in_catalog"]
GM12878.bed <- GM12878.bed[NIC,]
classification <- classification[NIC,]
report <- report[NIC,]

table(classification$subcategory)
classification$subcategory <- gsub("combination_of_known_junctions", "CJ",
                                           classification$subcategory)
classification$subcategory <- gsub("combination_of_known_splicesites", "CS",
                                           classification$subcategory)
classification$subcategory <- gsub("intron_retention", "IR",
                                           classification$subcategory)

report$nCell <- unlist(lapply(strsplit(report$sources,','), FUN = length))

GM12878.bed.out <- GM12878.bed
GM12878.bed.out$V4 <- paste0(GM12878.bed.out$V4, "_",
                        classification$subcategory, "_",
                        report$nCell)
write.table(GM12878.bed.out, "GM12878.NIC.bed",
            row.names = F, col.names = F,
            sep = "\t", quote = F)

GM12878.bed$subcategory <- classification$subcategory
GM12878.bed$nCell <- report$nCell
write.table(GM12878.bed, "GM12878.NIC.tsv",
            row.names = F, col.names = T,
            sep = "\t", quote = F)

## H9
H9.bed <- read.table("data/H9.c10.bed")
classification <- read.table("data/H9.SQANTI3_classification.txt", header = T)
report <- read.table("data/H9_trans_report.txt", header = T)

rownames(classification) <- classification$isoform
rownames(H9.bed) <- gsub("^.+;","",H9.bed$V4)
rownames(report) <- report$transcript_id

NIC <- classification$isoform[classification$structural_category == "novel_in_catalog"]
H9.bed <- H9.bed[NIC,]
classification <- classification[NIC,]
report <- report[NIC,]

table(classification$subcategory)
classification$subcategory <- gsub("combination_of_known_junctions", "CJ",
                                   classification$subcategory)
classification$subcategory <- gsub("combination_of_known_splicesites", "CS",
                                   classification$subcategory)
classification$subcategory <- gsub("intron_retention", "IR",
                                   classification$subcategory)

report$nCell <- unlist(lapply(strsplit(report$sources,','), FUN = length))

H9.bed.out <- H9.bed
H9.bed.out$V4 <- paste0(H9.bed.out$V4, "_",
                             classification$subcategory, "_",
                             report$nCell)
write.table(H9.bed.out, "H9.NIC.bed",
            row.names = F, col.names = F,
            sep = "\t", quote = F)

H9.bed$subcategory <- classification$subcategory
H9.bed$nCell <- report$nCell
write.table(H9.bed, "H9.NIC.tsv",
            row.names = F, col.names = T,
            sep = "\t", quote = F)

validated <- c("HNRNPA2B1","ATP5F1B","CMC2","RPS5","ITGB1BP1",
               "DYRK4","KRT7","RPLP2","FAM216A","SLC38A5","ZNF749",
               "LAMA5","PCBP2","PPP6R1")

rownames(gene.info.hg38) 
validated <- gene.info.hg38[validated,]


# 2. pseudogene expression ------------------------------------------------
load("E:/LabWork/Project/SCANSeq2/merge/Part1/Pseudogene/Pseudogene.RData")
bulk.express <- Seurat::AverageExpression(obj.human)
bulk.express <- as.data.frame(bulk.express$Gene)

gene.selcted <- c("AL355375.1", "HSPA8P13", "AC017037.2", "AL359643.1", "AL031864.2")
write.table(bulk.express[gene.selcted,],"pseudogene.selected.tsv", sep = "\t")

gene.selected <- pseudogenes$`Gene name`[pseudogenes$`Gene name` %in% rownames(bulk.express)
                                         & pseudogenes$Parent %in% rownames(bulk.express)]

## 2.1. GM12878 -----
pseudogenes.GM <- pseudogenes[pseudogenes$`Gene name` %in% gene.selected,]
rownames(pseudogenes.GM) <- pseudogenes.GM$`Gene name`
pseudogenes.GM$Expressed_pgene <- NULL
pseudogenes.GM$Expressed_parent <- NULL
pseudogenes.GM$Expression_pgene <- bulk.express[pseudogenes.GM$`Gene name`,]$GM12878
pseudogenes.GM$Expression_parent <- bulk.express[pseudogenes.GM$Parent,]$GM12878

pseudogenes.GM <- pseudogenes.GM[pseudogenes.GM$Expression_pgene +
                                   pseudogenes.GM$Expression_parent > 0,]
pseudogenes.GM$log2_fc <- log2((pseudogenes.GM$Expression_pgene + 0.01) /
                                 (pseudogenes.GM$Expression_parent + 0.01))
write.table(pseudogenes.GM, "pseudogenes.GM12878.tsv", sep = "\t",
          row.names = F, col.names = T, quote = F)

## 2.2. H9 ------
pseudogenes.H9 <- pseudogenes[pseudogenes$`Gene name` %in% gene.selected,]
rownames(pseudogenes.H9) <- pseudogenes.H9$`Gene name`
pseudogenes.H9$Expressed_pgene <- NULL
pseudogenes.H9$Expressed_parent <- NULL
pseudogenes.H9$Expression_pgene <- bulk.express[pseudogenes.H9$`Gene name`,]$H9
pseudogenes.H9$Expression_parent <- bulk.express[pseudogenes.H9$Parent,]$H9

pseudogenes.H9 <- pseudogenes.H9[pseudogenes.H9$Expression_pgene +
                                   pseudogenes.H9$Expression_parent > 0,]
pseudogenes.H9$log2_fc <- log2((pseudogenes.H9$Expression_pgene + 0.01) /
                                 (pseudogenes.H9$Expression_parent + 0.01))
write.table(pseudogenes.H9, "pseudogenes.H9.tsv", sep = "\t",
            row.names = F, col.names = T, quote = F)

# 3. DTU of GM12878 and H9 ------------------------------------------------
## 3.1. load data -----------
rc_trans_9CL <- readRDS("../../batch/1 9CL/data/rc_trans.RDS")
rc_trans_9Mix <- readRDS("../../batch/2 9CL_Mix/data/rc_trans.RDS")
identical(rownames(rc_trans_9CL),rownames(rc_trans_9Mix))
rc_trans <- cbind(rc_trans_9CL, rc_trans_9Mix)
rm(rc_trans_9CL, rc_trans_9Mix)

obj.9CL <- readRDS("E:/LabWork/Project/SCANSeq2/batch/1 9CL/obj.9CL.rds")
obj.9Mix <- readRDS("E:/LabWork/Project/SCANSeq2/batch/2 9CL_Mix/obj.9CL_Mix.rds")

obj.9CL[["Isoform"]] <- CreateAssayObject(rc_trans[,colnames(obj.9CL)])
obj.9Mix[["Isoform"]] <- CreateAssayObject(rc_trans[,colnames(obj.9Mix)])

obj.combined <- merge(obj.9CL, y = obj.9Mix, 
                      add.cell.ids = c("9CL", "9MIX"), project = "SCAN_Seq")
table(Idents(obj.combined))

obj.combined <- subset(obj.combined, subset = Cell_Line %in% c("GM12878","H9"))
table(obj.combined$Cell_Line)

rc_trans <- rc_trans[,obj.combined$Rename]

# sacle data
DefaultAssay(obj.combined) <- "Gene"
obj.combined <- NormalizeData(obj.combined)
obj.combined <- ScaleData(obj.combined)
DefaultAssay(obj.combined) <- "Isoform"
obj.combined <- NormalizeData(obj.combined)
obj.combined <- ScaleData(obj.combined)

## 3.2. DEG analysis ---------------
DefaultAssay(obj.combined) <- "Gene"
DEGs <-  FindMarkers(object = obj.combined,
                     ident.1 = "H9", ident.2 = "GM12878",
                     min.pct = 0.25,
                     logfc.threshold = 0.25,
                     only.pos = F)
DEGs <- DEGs[DEGs$p_val_adj<0.05,]

gene2trans <- gene2trans[,c(4,3,1,2)]
rownames(gene2trans) <- gene2trans$trans.symbol

trans <- intersect(rownames(gene2trans), rownames(rc_trans))
rc_trans <- rc_trans[trans,]
dim(rc_trans)

table(rownames(rc_trans) %in% gene2trans$trans.symbol)

## 3.3 DTU analysis --------------

sample.info <- obj.combined@meta.data
rownames(sample.info) <- sample.info$Rename
table(sample.info$Cell_Line)
table(colnames(rc_trans) %in% sample.info$Rename)

DTU <- run_drimseq(counts = as.matrix(rc_trans),
                   tx2gene = gene2trans,
                   pd=sample.info, id_col = "Rename",
                   cond_col = "Cell_Line",
                   cond_levels = c("H9", "GM12878"),
                   filtering_strategy = "own",
                   filter_only = F,
                   BPPARAM = BiocParallel::SnowParam(6),
                   min_samps_gene_expr = 25,
                   min_gene_expr = 5,
                   min_samps_feature_expr = 25,
                   min_feature_expr = 5)

DTU <- posthoc_and_stager(dturtle = DTU, ofdr = 0.05, posthoc = 0.1)
DTU$meta_table_gene[1:5,]
DTU$meta_table_tx[1:5,]
length(DTU$sig_gene)
length(DTU$sig_tx)

pdf("DGE_DTU.pdf",4,3)
plot(Venn(list("DGE" = rownames(DEGs),
               "DTU" = DTU$sig_gene)),
     doWeights = TRUE,
     show = list( SetLabels = T, Faces = FALSE))
dev.off()

DTU <- create_dtu_table(dturtle = DTU,
                        add_tx_metadata = list(tx_expr_in_max = c("exp_in", max),
                                               gene.id = c("gene.id", function(x){
                                                 return(names(which.max(table(x))))})))
head(DTU$dtu_table, n=20)

# temp <- plot_proportion_barplot(dturtle = DTU, genes = "PTPRC")
# temp$TFAP2A
# temp <- plot_proportion_pheatmap(dturtle = DTU, genes = "TFAP2A", 
#                                  sample_meta_table_columns = c("sample_id","seurat_clusters_2"),
#                                  include_expression = TRUE, treeheight_col=20,
#                                  BPPARAM = BiocParallel::SnowParam(1))
# temp$TFAP2A
# plot_transcripts_view(dturtle = DTU, 
#                       genes = "TFAP2A", 
#                       gtf = gtf, 
#                       genome = NULL, reduce_introns_min_size = 100,
#                       arrow_colors = c("#f46d43","#74add1"),
#                       one_to_one = TRUE,
#                       savepath = "./Hela",filename_ext  = "_transcripts.pdf",
#                       width = 5, height = 3)


## 3.4. output table ---------------------------------------------------------
bulk.gene <- AverageExpression(object = obj.combined, assays = "Gene")
bulk.isoform <- AverageExpression(object = obj.combined, assays = "Isoform")

FDR <- DTU$FDR_table
rownames(FDR) <- FDR$txID
out <- DTU$meta_table_tx[DTU$meta_table_tx$tx %in% 
                                DTU$sig_tx,]
out$gene_FDR_value <- FDR[out$tx,]$gene
out$tx_FDR_value <- FDR[out$tx,]$transcript
out$Transcript_Type <- trans.info.hg38[out$tx,]$`Transcript type`
out$Expression_H9 <- bulk.isoform$Isoform[out$tx,2]
out$Expression_GM12878 <- bulk.isoform$Isoform[out$tx,1]

gene <- gsub("-.+?$","",out$tx)
table(gene %in% rownames(bulk.gene$Gene))
which(!gene %in% rownames(bulk.gene$Gene))
out$tx[80]
gene[80] <- "H3-3A"
gene[!gene %in% rownames(bulk.gene$Gene)]

out$Expression_Gene_H9 <- bulk.gene$Gene[gene,2]
out$Expression_Gene_GM12878 <- bulk.gene$Gene[gene,1]
out$Percent_H9 <- out$Expression_H9 / out$Expression_Gene_H9
out$Percent_GM12878 <- out$Expression_GM12878 / out$Expression_Gene_GM12878

write.table(out,
            "DTU.trans.tsv", sep = "\t")


out2 <- DTU$dtu_table
out2$Expression_H9 <- bulk.gene$Gene[out2$gene_ID,2]
out2$Expression_GM12878 <- bulk.gene$Gene[out2$gene_ID,1]
write.table(out2,
            "DTU.gene.tsv", sep = "\t")
View(DTU$dtu_table)

rm(rc_trans, obj.9CL, obj.9Mix, obj.combined)

# GAPDH
isoforms <- grep("GAPDH-", rownames(bulk.isoform$Isoform))
GAPDH <- bulk.isoform$Isoform[isoforms,]
GAPDH <- GAPDH[rowSums(GAPDH)>0,]
bulk.gene$Gene["GAPDH",]


# RANGRF
isoforms <- grep("RANGRF-", rownames(bulk.isoform$Isoform))
RANGRF <- bulk.isoform$Isoform[isoforms,]
RANGRF <- RANGRF[rowSums(RANGRF)>0,]
bulk.gene$Gene["GAPDH",]


# House Keeping
temp <- as.data.frame(bulk.gene$Gene %>% as.data.frame %>%.[.$H9>0 & .$GM12878>0, ])
temp$Mean_Expression <- rowMeans(temp)
temp$Ratio_GM_vs_H9 <- temp$GM12878 / temp$H9
temp <- temp[abs(temp$Ratio_GM_vs_H9 - 1) < 0.05,]
temp <- temp[order(temp$Mean_Expression, decreasing = T),]
temp$Gene <- rownames(temp)
write.csv(temp, "Control_Genes.csv", row.names = F, col.names = T)
