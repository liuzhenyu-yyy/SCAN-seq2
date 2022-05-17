setwd("E:/LabWork/Project/SCANSeq2/merge/Part2/")

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

sample.info.9CL <- read.table("../../batch/1 9CL/sample.info.txt")
sample.info.9Mix <- read.table("../../batch/2 9CL_Mix/sample.info.txt")
sample.info.4CL <- read.table("../../batch/3 4CL/sample.info.txt")
sample.info.9CL$Cell_Line <- factor(sample.info.9CL$Cell_Line,
                                levels = c("GM12878","HepG2","Hela","H9","K562","293T",
                                           "AtT20","MEF","3T3"))
sample.info.9Mix$Cell_Line <- factor(sample.info.9Mix$Cell_Line,
                                    levels = c("GM12878","HepG2","Hela","H9","K562","293T",
                                               "AtT20","MEF","3T3"))
sample.info.4CL$Cell_Line <- factor(sample.info.4CL$Cell_Line,
                                    levels = c("GM12878","HepG2","Hela","H9","K562","293T",
                                               "AtT20","MEF","3T3"))
sample.info.9CL <- sample.info.9CL[sample.info.9CL$Pass_QC > 0,]
sample.info.9Mix <- sample.info.9Mix[sample.info.9Mix$Pass_QC > 0,]
sample.info.4CL <- sample.info.4CL[sample.info.4CL$Pass_QC > 0,]

sample.info <- rbind(sample.info.9CL, sample.info.9Mix, sample.info.4CL)
sample.info$Library <- factor(sample.info$Library, levels = c("9CL","9CL_Mix","4CL"))
# color for cell lines
color.cell <- brewer.pal(11,"Set3")[c(1:8,10)]
names(color.cell) <- c("GM12878","HepG2","Hela","H9","K562","293T",
                    "AtT20","MEF","3T3")
color.cell["HepG2"] <- "#ffed6f"


# 1. single cell assembly -------------------------------------------------
dir.create("SingleCell")
# add assembly results
isoform.count.9CL <- read.table("../../batch/1 9CL/data/Isoform.category.tsv")
isoform.count.9Mix <- read.table("../../batch/2 9CL_Mix/data/Isoform.category.tsv")
isoform.count.4CL <- read.table("../../batch/3 4CL/data/Isoform.category.tsv")
colnames(isoform.count.9CL) <- c("Name","Isoforms","FSM","ISM","NIC", "NNC","Genic_Genome",
                                 "Antisense", "Fusion" ,"Intergenic","Genic_intron",
                                 "NIC_CJ","NIC_CS", "NIC_IR")
colnames(isoform.count.9Mix) <- c("Name","Isoforms","FSM","ISM","NIC", "NNC","Genic_Genome",
                                  "Antisense", "Fusion" ,"Intergenic","Genic_intron",
                                  "NIC_CJ","NIC_CS", "NIC_IR")
colnames(isoform.count.4CL) <- c("Name","Isoforms","FSM","ISM","NIC", "NNC","Genic_Genome",
                                 "Antisense", "Fusion" ,"Intergenic","Genic_intron",
                                 "NIC_CJ","NIC_CS", "NIC_IR")
rownames(isoform.count.9CL) <- isoform.count.9CL$Name
rownames(isoform.count.9Mix) <- isoform.count.9Mix$Name
rownames(isoform.count.4CL) <- isoform.count.4CL$Name
isoform.count.9CL <- isoform.count.9CL[sample.info.9CL$Name,]
isoform.count.9Mix <- isoform.count.9Mix[sample.info.9Mix$Name,]
isoform.count.4CL <- isoform.count.4CL[sample.info.4CL$Name,]
rownames(isoform.count.9CL) <- sample.info.9CL$Rename
rownames(isoform.count.9Mix) <- sample.info.9Mix$Rename
rownames(isoform.count.4CL) <- sample.info.4CL$Rename

isoform.count <- rbind(isoform.count.9CL, isoform.count.9Mix, isoform.count.4CL)
identical(rownames(isoform.count), rownames(sample.info))
sample.info <- cbind(sample.info, isoform.count[,2:14])

rm(isoform.count.9CL, isoform.count.9Mix,isoform.count.4CL)
rm(sample.info.9CL, sample.info.9Mix, sample.info.4CL)

## 1.1. Number of isoform in each cell--------------
sample.info$percent.FSM <- sample.info$FSM / sample.info$Isoforms
sample.info$percent.ISM <- sample.info$ISM / sample.info$Isoforms
sample.info$percent.NIC <- sample.info$NIC / sample.info$Isoforms
sample.info$percent.NNC <- sample.info$NNC / sample.info$Isoforms
sample.info$percent.CJ <- sample.info$NIC_CJ / sample.info$Isoforms
sample.info$percent.CS <- sample.info$NIC_CS / sample.info$Isoforms
sample.info$percent.IR <- sample.info$NIC_IR / sample.info$Isoforms
sample.info$percent.NIC.CJ <- sample.info$NIC_CJ / sample.info$NIC
sample.info$percent.NIC.CS <- sample.info$NIC_CS / sample.info$NIC
sample.info$percent.NIC.IR <- sample.info$NIC_IR / sample.info$NIC

pdf("SingleCell/Vln_nGene_Isoform.pdf",6,2)
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Gene_Detected/1000,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Gene_Detected/1000,color=Cell_Line),alpha=0,show.legend = F)+
  ylab("Getected Genes (10^3)")+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Trans_Detected/1000,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Trans_Detected/1000,color=Cell_Line),alpha=0,show.legend = F)+
  ylab("Getected Isoforms (10^3)")+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Isoforms/1000,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Isoforms/1000,color=Cell_Line),alpha=0,show.legend = F)+
  ylab("Assembled Isoforms (10^3)")+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  theme_lzy()
dev.off()

pdf("SingleCell/Vln_nGene_Isoform.split.pdf",6,5)
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Gene_Detected/1000,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Gene_Detected/1000,color=Cell_Line),alpha=0,show.legend = F)+
  ylab("Getected Genes (10^3)")+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  facet_wrap(~Library,ncol = 1)+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Trans_Detected/1000,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Trans_Detected/1000,color=Cell_Line),alpha=0,show.legend = F)+
  ylab("Getected Isoforms (10^3)")+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  facet_wrap(~Library,ncol = 1)+
  theme_lzy()
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Isoforms/1000,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Isoforms/1000,color=Cell_Line),alpha=0,show.legend = F)+
  ylab("Assembled Isoforms (10^3)")+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  facet_wrap(~Library,ncol = 1)+
  theme_lzy()
dev.off()

pdf("SingleCell/Vln_nGene_Isoform.merge.pdf",6,1.5)
temp <- aggregate(sample.info$Isoforms,list(sample.info$Cell_Line),median)
ggplot(sample.info)+
  geom_violin(aes(x=Cell_Line,y=Isoforms/1000,fill=Cell_Line),alpha=0.7,show.legend = F)+
  geom_boxplot(aes(x=Cell_Line,y=Isoforms/1000,color=Cell_Line),alpha=0,show.legend = F)+
  ylab("Assembled Isoforms (10^3)")+
  xlab(NULL)+
  scale_color_manual(values = color.cell)+
  scale_fill_manual(values = color.cell)+
  scale_x_discrete(breaks = temp$Group.1,
                   labels = paste0(temp$Group.1,"\n",signif(temp$x,4)))+
  theme_lzy()
dev.off()

aggregate(sample.info$Isoforms, list(sample.info$Cell_Line), mean)
aggregate(sample.info$Isoforms, list(sample.info$Cell_Line), median)
table(sample.info$Cell_Line)

## 1.2. percent of each category-------------------------------------------------
features <- grep("percent.",colnames(sample.info))
pdf("SingleCell/Vln_percent.category.pdf",6,2)
for (one in features) {
  p <- ggplot(sample.info)+
    geom_violin(aes(x=Cell_Line,y=sample.info[,one]*100,fill=Cell_Line),alpha=0.7,show.legend = F)+
    geom_boxplot(aes(x=Cell_Line,y=sample.info[,one]*100,color=Cell_Line),alpha=0,show.legend = F)+
    ylab(colnames(sample.info)[one])+
    scale_color_manual(values = color.cell)+
    scale_fill_manual(values = color.cell)+
    theme_lzy()
  plot(p)
}
dev.off()


# 2. cell-line assembly ---------------------------------------------------
dir.create("CellLine")

# merge cell line assemblys
assembly.CL.bed <- data.frame()

for (one in names(color.cell)) {
  temp <- as.data.frame(fread(paste0("data/cell_line/",one,"_trans_report.txt")))
  temp$Cell_Line <- one
  assembly.CL.bed <- rbind(assembly.CL.bed, temp)
  message(paste0("Finish reading assembly of ", one, "."))
}
table(assembly.CL.bed$Cell_Line)

temp <- unlist(lapply(strsplit(assembly.CL.bed$sources,','), FUN = length))
table(assembly.CL.bed$num_clusters - temp)
assembly.CL.bed$num_cells <- temp

## 2.1. threshold for cell number------------------------------------------
table(sample.info$Cell_Line)
table(sample.info$Library)

breaks <- 1:50
nIsoform <- data.frame()
for (one in names(color.cell)) {
  nCells <- assembly.CL.bed[assembly.CL.bed$Cell_Line==one,]$num_cells
  temp <- unlist(lapply(breaks, function(x){return(sum(nCells >= x))}))
  nIsoform_one <- data.frame("Threshold"=breaks,
                             "nIsoform"=temp,
                             "Cell_Line"=rep(one))
  nIsoform <- rbind(nIsoform,nIsoform_one)
  rm(nCells, temp, nIsoform_one)
  message(paste0("Finish filter assembly of ", one, "."))
}

pdf("CellLine/Threshold.nIsoforms.pdf",4,3)
ggplot(nIsoform[nIsoform$Threshold <= 20,], aes(x = Threshold, y = nIsoform / 1000, color = Cell_Line))+
  geom_point()+  geom_line()+
  scale_color_manual(values = color.cell)+
  ylab("Number of Isoforms (10^3)")+
  xlab("Number of Cells")+
  theme_classic()
ggplot(nIsoform[nIsoform$Threshold <= 40,], aes(x = Threshold, y = log10(nIsoform), color = Cell_Line))+
  geom_point(size = 0.5)+  geom_line()+
  scale_color_manual(values = color.cell)+
  ylab("log10(Number of Isoforms)")+
  xlab("Number of Cells")+
  theme_classic()
dev.off()

nIsoform[nIsoform$Threshold==10,]

assembly.CL.filtered.bed <- assembly.CL.bed[assembly.CL.bed$num_cells >= 10,]
saveRDS(assembly.CL.bed,"assembly.CL.bed.rds")
rm(assembly.CL.bed)
gc()

for (one in names(color.cell)) {
  bed <- as.data.frame(fread(paste0("data/cell_line/",one,".bed")))
  trans.id <- strsplit(bed$V4,";")
  trans.id <- unlist(lapply(trans.id, function(x){return(x[2])}))
  if(length(trans.id) != nrow(bed)){
    stop("error in extracting trans.id")   
  }
  
  trans.select <- assembly.CL.filtered.bed[assembly.CL.filtered.bed$Cell_Line==one,]$transcript_id
  bed <- bed[trans.id %in% trans.select,]
  write.table(bed, paste0("data/cell_line/",one,".c10.bed"),
              row.names = F, col.names = F,
              sep = "\t", quote = F)
  message(paste0("Finish writing bed file of ", one, "."))
}
rm(bed, trans.select, trans.id)


# 3. merged assembly ------------------------------------------------------
dir.create("Merge")

# load data
isoform.mouse <- read.table("data/merge/merge_mouse_3CL.SQANTI3_classification.txt",
                            header = T, stringsAsFactors = F)
temp <- read.table("data/merge/merge_mouse_3CL_trans_report.txt",
                                    header = T, stringsAsFactors = F)
table(temp$transcript_id%in%isoform.mouse$isoform)
rownames(temp) <- temp$transcript_id
temp <- temp[isoform.mouse$isoform,]
identical(temp$transcript_id,isoform.mouse$isoform)
isoform.mouse$source <- temp$sources
isoform.mouse$Organism <- "Mouse"
rm(temp)

isoform.human <- read.table("data/merge/merge_human_6CL.SQANTI3_classification.txt",
                            header = T, stringsAsFactors = F)
temp <- read.table("data/merge/merge_human_6CL_trans_report.txt",
                   header = T, stringsAsFactors = F)
table(temp$transcript_id%in%isoform.human$isoform)
rownames(temp) <- temp$transcript_id
temp <- temp[isoform.human$isoform,]
identical(temp$transcript_id,isoform.human$isoform)
isoform.human$source <- temp$sources
isoform.human$Organism <- "Human"
rm(temp)


## 3.1. basic statistics for gene, isoform and junctions -----------------------
isoform.merge <- rbind(isoform.mouse, isoform.human)
table(isoform.merge$structural_category, isoform.merge$Organism)
temp <- isoform.merge$structural_category
temp <- gsub("antisense","Antisense",temp)
temp <- gsub("full-splice_match","FSM",temp)
temp <- gsub("fusion","Fusion",temp)
temp <- gsub("genic","Genic",temp)
temp <- gsub("incomplete-splice_match","ISM",temp)
temp <- gsub("interGenic","Intergenic",temp)
temp <- gsub("novel_in_catalog","NIC",temp)
temp <- gsub("novel_not_in_catalog","NNC",temp)
isoform.merge$structural_category <- factor(temp,
                                            levels = c("FSM","ISM","NIC","NNC",
                                                       "Genic","Antisense","Fusion",
                                                       "Intergenic"))
table(isoform.merge$structural_category, isoform.merge$Organism)
isoform.merge$Organism <- factor(isoform.merge$Organism,
                                 levels = c("Human","Mouse"))

isoform.merge$nCell_Line <- as.numeric(unlist(lapply(strsplit(isoform.merge$source,","),length)))

color.isoform <- c("#6baed6","#fc8d59","#78c679","#ee6a50",
                   "#969696","#66c2a4","#ffc125","#e9967a")
names(color.isoform) <- levels(isoform.merge$structural_category)

plot.data <- as.data.frame(table(isoform.merge$structural_category,
                                 isoform.merge$Organism,
                                 isoform.merge$nCell_Line))
colnames(plot.data) <- c("structural_category", "Organism", "nCell_Line","Frequency")
plot.data$nCell_Line <- factor(as.numeric(plot.data$nCell_Line), levels = 6:1)

pdf("Merge/Bar.percent.isoform.pdf",5,3)
ggplot(plot.data[plot.data$Organism == "Human",])+
  geom_bar(aes(x = structural_category,y =  Frequency / sum(Frequency) *100,
               fill = structural_category,
               alpha = nCell_Line),
           stat = "identity", position="stack", show.legend = F)+
  scale_fill_manual(values = color.isoform)+
  scale_alpha_manual(values =  1 - 0.13 * 0:5)+
  ylab("Percent of Isoforms")+
  theme_classic()
ggplot(plot.data[plot.data$Organism == "Mouse",])+
  geom_bar(aes(x = structural_category,y =  Frequency / sum(Frequency) *100,
               fill = structural_category,
               alpha = nCell_Line),
           stat = "identity", position="stack", show.legend = F)+
  scale_fill_manual(values = color.isoform)+
  scale_alpha_manual(values =  1 - 0.25 * c(0,0,0,0,1,2))+
  ylab("Percent of Isoforms")+
  theme_classic()
dev.off()

table(isoform.merge$Organism,
      isoform.merge$structural_category) / as.numeric(table(isoform.merge$Organism))

# junctions
table(junctions.human$canonical,junctions.human$junction_category)


## 3.2. detailed classification of NIC -------------------------------------
isoform.merge.NIC <- isoform.merge[isoform.merge$structural_category == "NIC", ]
table(isoform.merge.NIC$Organism)
table(isoform.merge.NIC$subcategory, isoform.merge.NIC$Organism)
temp <- isoform.merge.NIC$subcategory
temp <- gsub("combination_of_known_junctions","CJ",temp)
temp <- gsub("combination_of_known_splicesites","CS",temp)
temp <- gsub("intron_retention","IR",temp)
isoform.merge.NIC$subcategory <- factor(temp,
                                        levels = c("CJ","CS","IR"))

color.NIC <- brewer.pal(3,"Set2")
names(color.NIC) <- c("CJ","CS","IR")

pdf("Merge/Pie.NIC.percent.pdf",4,3)
plot.data <- as.data.frame(table(isoform.merge.NIC[isoform.merge.NIC$Organism == "Human",]$subcategory))
plot.data$label <- paste0(plot.data$Var1, "\n",
                         signif(plot.data$Freq / sum(plot.data$Freq)*100, 3), "%")
ggpubr::ggdonutchart(plot.data, "Freq", label = "label",                               
                     fill = "Var1", color = "white", lab.pos = "in",
                     palette = color.NIC)

plot.data <- as.data.frame(table(isoform.merge.NIC[isoform.merge.NIC$Organism == "Mouse",]$subcategory))
plot.data$label <- paste0(plot.data$Var1, "\n",
                          signif(plot.data$Freq / sum(plot.data$Freq)*100, 3), "%")
ggpubr::ggdonutchart(plot.data, "Freq", label = "label",                               
                     fill = "Var1", color = "white", lab.pos = "in",
                     palette = color.NIC)
dev.off()

rm(isoform.human, isoform.mouse)


# 4. B cell receptor ---------------------------------------------------
dir.create("BCR")

## 4.1. constant region ----------------------------------------------------
obj.9CL <- readRDS("E:/LabWork/Project/SCANSeq2/batch/1 9CL/obj.9CL.rds")
obj.9Mix <- readRDS("E:/LabWork/Project/SCANSeq2/batch/2 9CL_Mix/obj.9CL_Mix.rds")
obj.combined <- merge(obj.9CL, y = obj.9Mix, 
                      add.cell.ids = c("9CL", "9MIX"), project = "SCAN_Seq")
table(sample.info$Cell_Line)

IG_C_genes <- sort(gene.info.hg38[gene.info.hg38$`Gene type`=="IG_C_gene",]$`Gene name`)

pdf("BCR/Vln.IGH.constant.pdf",12,4.5)
VlnPlot(obj.combined,
        grep("IGH",IG_C_genes,value = T),
        idents = c("GM12878", "H9", "K562"),
        cols = color.cell,ncol = 5,same.y.lims = T)
dev.off()
pdf("BCR/Vln.IGL.constant.pdf",12,2.3)
VlnPlot(obj.combined,
        grep("IGH",IG_C_genes,value = T,invert = T),
        idents = c("GM12878", "H9", "K562"),
        cols = color.cell,ncol = 5,same.y.lims = T)
dev.off()

rc <- as.data.frame(obj.combined@assays$Gene@counts)
rc <- rc[,obj.combined$Cell_Line=="GM12878"]

rowSums(rc > 1)["IGHD"] # 88 of 187
rowSums(rc > 1)["IGHM"] # 186 of 187

rm(rc)

## 4.2. V(D)J Recombination ------------------------------------------------
# VDJFilter <- function(temp, min.reads = 10, E.cutoff = 1){
#   if(nrow(temp)< min.reads){
#     temp$v_call <- "none"
#     temp$d_call <- "none"
#     temp$j_call <- "none"
#     return(temp[1,])
#   }
#   temp[is.na(temp)] <- 100
#   
#   res <- temp[1,]
#   res$v_support <- temp$v_support[which.min(temp$v_support)]
#   res$v_score <- temp$v_score[which.min(temp$v_support)]
#   res$v_call <- temp$v_call[which.min(temp$v_support)]
#   if(res$v_support>E.cutoff){res$v_call <- "none"}
#   res$d_support <- temp$d_support[which.min(temp$d_support)]
#   res$d_score <- temp$d_score[which.min(temp$d_support)]
#   res$d_call <- temp$d_call[which.min(temp$d_support)]
#   if(res$d_support>E.cutoff){res$d_call <- "none"}
#   res$d_support <- temp$d_support[which.min(temp$j_support)]
#   res$j_score <- temp$j_score[which.min(temp$j_support)]
#   res$j_call <- temp$j_call[which.min(temp$j_support)]
#   if(res$j_support>E.cutoff){res$j_call <- "none"}
#   return(res)
# }

SimplifyVDJCall <- function(x){
  x <- strsplit(x, ",")
  x <- lapply(x, function(x){return(gsub("_.+?$", "",
                                         gsub("*","_",x, fixed = T)))})
  x <- lapply(x, function(x){return(names(sort(table(x), decreasing = T))[1])})
  return(unlist(x))
}

PlotVDJDonut <- function(data,col){
  plot.data <- as.data.frame(table(data[,col]))
  plot.data$Var1 <- as.character(plot.data$Var1)
  plot.data <- plot.data[order(plot.data$Freq, decreasing = T),]
  plot.data$Label <- paste(plot.data$Var1,plot.data$Freq,sep = "\n")
  plot.data$Color <- as.character(plot.data$Var1)
  if (min(plot.data$Freq )<= 5) {
    plot.data[plot.data$Freq <= 5 ,]$Label <- NA
  }
  plot.data$Color <- factor(plot.data$Color,
                            levels = rev(unique(rev(c(plot.data$Color,"none")))))
  
  colors <- rep(brewer.pal(8, "Pastel2"),3)[1:length(plot.data$Color)]
  names(colors) <- as.character(plot.data$Color)
  colors["none"] <- "gray70"
  p <- ggpubr::ggdonutchart(plot.data, x = "Freq", label = "Label",                               
                            fill = "Color", color = "white", lab.pos = "in",
                            palette = colors)
  return(p)
}

cell.list <- rownames(sample.info)[sample.info$Cell_Line=="GM12878"]
cell.list <- gsub("9CL_Mix","9Mix",cell.list)
VDJ.IGH <- data.frame()
VDJ.IGL <- data.frame()

for (one in cell.list) {
  if(file.exists(paste0("data/bcr/",one,".IgBlast.IGH.tsv"))){
    temp <- read.table(paste0("data/bcr/",one,".IgBlast.IGH.tsv"), header = T, sep = "\t")
    #temp <- VDJFilter(temp)
    temp$Cell <- one
    VDJ.IGH <- rbind(VDJ.IGH, temp)
  }
  
  if (file.exists(paste0("data/bcr/",one,".IgBlast.IGL.tsv"))) {
    temp <- read.table(paste0("data/bcr/",one,".IgBlast.IGL.tsv"), header = T, sep = "\t")
    #temp <- VDJFilter(temp)
    temp$Cell <- one
    VDJ.IGL <- rbind(VDJ.IGL, temp)
  }
  
  message(paste0("Filtering VDJ for ",one,"..."))
}

empty <- cell.list[!cell.list%in%VDJ.IGH$Cell]
VDJ.IGH[(nrow(VDJ.IGH)+1):(nrow(VDJ.IGH)+length(empty)),] <- "none"
VDJ.IGH[(nrow(VDJ.IGH)-length(empty)+1):(nrow(VDJ.IGH)),]$Cell <- empty
empty <- cell.list[!cell.list%in%VDJ.IGL$Cell]
VDJ.IGL[(nrow(VDJ.IGL)+1):(nrow(VDJ.IGL)+length(empty)),] <- "none"
VDJ.IGL[(nrow(VDJ.IGL)-length(empty)+1):(nrow(VDJ.IGL)),]$Cell <- empty

rownames(VDJ.IGH) <- VDJ.IGH$Cell
rownames(VDJ.IGL) <- VDJ.IGL$Cell
VDJ.IGH <- VDJ.IGH[cell.list,]
VDJ.IGL <- VDJ.IGL[cell.list,]
VDJ.IGH[is.na(VDJ.IGH)] <- "none"
VDJ.IGH[VDJ.IGH==""] <- "none"
VDJ.IGL[is.na(VDJ.IGL)] <- "none"
VDJ.IGL[VDJ.IGL==""] <- "none"

BCR.stat <- read.table("data/BCR.VDJ.stat.tsv")
colnames(BCR.stat) <- c("Cell","IGL_Reads","IGH_Reads","IGL_Cluster","IGH_Cluster")
rownames(BCR.stat) <- BCR.stat$Cell
identical(BCR.stat$Cell,VDJ.IGL$Cell)
identical(BCR.stat$Cell,VDJ.IGH$Cell)
VDJ.IGH[BCR.stat$IGH_Reads < 5,c("v_call","d_call","j_call")] <- "none"
VDJ.IGL[BCR.stat$IGL_Reads < 5,c("v_call","d_call","j_call")] <- "none"

VDJ.IGH$v_call_sim <- SimplifyVDJCall(VDJ.IGH$v_call)
VDJ.IGH$d_call_sim <- SimplifyVDJCall(VDJ.IGH$d_call)
VDJ.IGH$j_call_sim <- SimplifyVDJCall(VDJ.IGH$j_call)
VDJ.IGL$v_call_sim <- SimplifyVDJCall(VDJ.IGL$v_call)
VDJ.IGL$j_call_sim <- SimplifyVDJCall(VDJ.IGL$j_call)


p1 <- PlotVDJDonut(VDJ.IGH,"v_call_sim")+NoLegend()
p2 <- PlotVDJDonut(VDJ.IGH,"d_call_sim")+NoLegend()
p3 <- PlotVDJDonut(VDJ.IGH,"j_call_sim")+NoLegend()
p4 <- PlotVDJDonut(VDJ.IGL,"v_call_sim")+NoLegend()
p5 <- PlotVDJDonut(VDJ.IGL,"j_call_sim")+NoLegend()

pdf("BCR/Pie.IGH.VDJ.pdf",3,3)
plot(p1)
plot(p2)
plot(p3)
dev.off()

pdf("BCR/Pie.IGL.VDJ.pdf",3,3)
plot(p4)
plot(p5)
dev.off()

# 4.3. define subclone relationship ---------------------------------------
VDJ.IGH$Freq <- 1
VDJ.IGH$subclone <- paste(VDJ.IGH$v_call_sim, VDJ.IGH$d_call_sim, VDJ.IGH$j_call_sim,
                          sep = "_")

VDJ.IGH.sub <- VDJ.IGH[VDJ.IGH$v_call_sim != "none" &
                         VDJ.IGH$d_call_sim != "none" &
                         VDJ.IGH$j_call_sim != "none" ,]
library(ggsankey)
subclones <- sort(table(VDJ.IGH.sub$subclone),decreasing = T)
subclones <- names(subclones)[subclones>=3]
names(subclones) <- paste("Subclone",1:3,sep="_")

VDJ.IGH.sub[VDJ.IGH.sub$subclone == subclones[1],]$subclone <- names(subclones)[1]
VDJ.IGH.sub[VDJ.IGH.sub$subclone == subclones[2],]$subclone <- names(subclones)[2]
VDJ.IGH.sub[VDJ.IGH.sub$subclone == subclones[3],]$subclone <- names(subclones)[3]
VDJ.IGH.sub[!VDJ.IGH.sub$subclone %in% names(subclones),]$subclone <- "others"
table(VDJ.IGH.sub$subclone)

VDJ.IGH.sub$v_call_sim <- factor(VDJ.IGH.sub$v_call_sim,
                                 levels = names(sort(table(VDJ.IGH.sub$v_call_sim),
                                                     decreasing = T)))
VDJ.IGH.sub$d_call_sim <- factor(VDJ.IGH.sub$d_call_sim,
                                 levels = names(sort(table(VDJ.IGH.sub$d_call_sim),
                                                     decreasing = T)))
VDJ.IGH.sub$j_call_sim <- factor(VDJ.IGH.sub$j_call_sim,
                                 levels = names(sort(table(VDJ.IGH.sub$j_call_sim),
                                                     decreasing = T)))

plot.data <- VDJ.IGH.sub %>%
  make_long(v_call_sim, d_call_sim, j_call_sim)
plot.data$node <- factor(plot.data$node, levels = rev(c(levels(VDJ.IGH.sub$v_call_sim),
                                                        levels(VDJ.IGH.sub$d_call_sim),
                                                        levels(VDJ.IGH.sub$j_call_sim))))

color.subclone <- c(brewer.pal(3,"Pastel2"),"gray75")
names(color.subclone) <- c(names(subclones),"others")
sort(table(plot.data$node),decreasing = T)

temp <- color.subclone[rep(c("Subclone_1", "Subclone_2", "Subclone_3","others"),time = c(3,3,3,4))]
names(temp) <-names(sort(table(plot.data$node),decreasing = T)) 

pdf("BCR/Sankey.subclone.GM12878.pdf",6,4)
ggplot(plot.data, aes(x = x, next_x = next_x,
                      node = node, next_node = next_node,
                      fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .8, smooth=8,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray30") +
  scale_x_discrete(breaks = c("v_call_sim","d_call_sim","j_call_sim"),
                   labels =c("IGHV", "IGHD", "IGHJ")) +
  # scale_fill_viridis_d(option = "D", alpha = .8)+
  scale_fill_manual(values = temp)+
  theme_sankey(base_size = 15) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("GM12878 subclones")
dev.off()

rm(p1,p2,p3,p4,p5,p,plot.data,obj.9CL,obj.9Mix,temp)

# 5. T cell receptor ---------------------------------------------------
dir.create("TCR")

## 5.1. constant region ----------------------------------------------------
TR_C_genes <- sort(gene.info.hg38[gene.info.hg38$`Gene type`=="TR_C_gene",]$`Gene name`)

pdf("TCR/Vln.TC.constant.pdf",7,4.5)
VlnPlot(obj.combined,
        grep("*",TR_C_genes,value = T),
        idents = c("GM12878", "H9", "K562"),
        cols = color.cell,ncol = 3,same.y.lims = T)
dev.off()

rc <- as.data.frame(obj.combined@assays$Gene@counts)
rc <- rc[,obj.combined$Cell_Line=="H9"]

rowSums(rc > 1)["TRBC1"] # 240  of 271
rowSums(rc > 1)["TRGC1"] # 38 of 271
table(as.numeric(rc["TRBC1",]>1), as.numeric(rc["TRGC1",]>1))

rm(rc)

## 4.2. V(D)J Recombination ------------------------------------------------
cell.list <- rownames(sample.info)[sample.info$Cell_Line=="H9"]
cell.list <- gsub("9CL_Mix","9Mix",cell.list)
VDJ.TRB <- data.frame()
VDJ.TRA <- data.frame()

for (one in cell.list) {
  if(file.exists(paste0("data/tcr/",one,".IgBlast.TRB.tsv"))){
    temp <- read.table(paste0("data/tcr/",one,".IgBlast.TRB.tsv"), header = T, sep = "\t")
    #temp <- VDJFilter(temp)
    temp$Cell <- one
    VDJ.TRB <- rbind(VDJ.TRB, temp)
  }

  if (file.exists(paste0("data/tcr/",one,".IgBlast.TRA.tsv"))) {
    temp <- read.table(paste0("data/tcr/",one,".IgBlast.TRA.tsv"), header = T, sep = "\t")
    #temp <- VDJFilter(temp)
    temp$Cell <- one
    VDJ.TRA <- rbind(VDJ.TRA, temp)
  }

  message(paste0("Filtering VDJ for ",one,"..."))
}

empty <- cell.list[!cell.list%in%VDJ.TRB$Cell]
VDJ.TRB[(nrow(VDJ.TRB)+1):(nrow(VDJ.TRB)+length(empty)),] <- "none"
VDJ.TRB[(nrow(VDJ.TRB)-length(empty)+1):(nrow(VDJ.TRB)),]$Cell <- empty
empty <- cell.list[!cell.list%in%VDJ.TRA$Cell]
VDJ.TRA[(nrow(VDJ.TRA)+1):(nrow(VDJ.TRA)+length(empty)),] <- "none"
VDJ.TRA[(nrow(VDJ.TRA)-length(empty)+1):(nrow(VDJ.TRA)),]$Cell <- empty

rownames(VDJ.TRB) <- VDJ.TRB$Cell
rownames(VDJ.TRA) <- VDJ.TRA$Cell
VDJ.TRB <- VDJ.TRB[cell.list,]
VDJ.TRA <- VDJ.TRA[cell.list,]
VDJ.TRB[is.na(VDJ.TRB)] <- "none"
VDJ.TRB[VDJ.TRB==""] <- "none"
VDJ.TRA[is.na(VDJ.TRA)] <- "none"
VDJ.TRA[VDJ.TRA==""] <- "none"

TCR.stat <- read.table("data/TCR.VDJ.stat.tsv")
colnames(TCR.stat) <- c("Cell","TRA_Reads","TRB_Reads","TRA_Cluster","TRB_Cluster")
rownames(TCR.stat) <- TCR.stat$Cell
identical(TCR.stat$Cell,VDJ.TRA$Cell)
identical(TCR.stat$Cell,VDJ.TRB$Cell)
VDJ.TRA[TCR.stat$TRA_Reads < 5,c("v_call","d_call","j_call")] <- "none"
VDJ.TRB[TCR.stat$TRB_Reads < 5,c("v_call","d_call","j_call")] <- "none"


VDJ.TRB$v_call_sim <- SimplifyVDJCall(VDJ.TRB$v_call)
VDJ.TRB$d_call_sim <- SimplifyVDJCall(VDJ.TRB$d_call)
VDJ.TRB$j_call_sim <- SimplifyVDJCall(VDJ.TRB$j_call)
VDJ.TRA$v_call_sim <- SimplifyVDJCall(VDJ.TRA$v_call)
VDJ.TRA$j_call_sim <- SimplifyVDJCall(VDJ.TRA$j_call)

p1 <- PlotVDJDonut(VDJ.TRB,"v_call_sim")+NoLegend()
p2 <- PlotVDJDonut(VDJ.TRB,"d_call_sim")+NoLegend()
p3 <- PlotVDJDonut(VDJ.TRB,"j_call_sim")+NoLegend()
p4 <- PlotVDJDonut(VDJ.TRA,"v_call_sim")+NoLegend()
p5 <- PlotVDJDonut(VDJ.TRA,"j_call_sim")+NoLegend()

pdf("TCR/Pie.TRB.VDJ.pdf",3,3)
plot(p1)
plot(p2)
plot(p3)
dev.off()

pdf("TCR/Pie.TRA.VDJ.pdf",3,3)
plot(p4)
plot(p5)
dev.off()

# 5.3. define subclone relationship ---------------------------------------
VDJ.TRA$Freq <- 1
VDJ.TRA$subclone <- paste(VDJ.TRA$v_call_sim, VDJ.TRA$j_call_sim,
                          sep = "_")

VDJ.TRA.sub <- VDJ.TRA[VDJ.TRA$v_call_sim != "none" &
                         VDJ.TRA$j_call_sim != "none" ,]
library(ggsankey)
subclones <- sort(table(VDJ.TRA.sub$subclone),decreasing = T)
subclones <- names(subclones)[subclones>=3]
names(subclones) <- paste("Subclone",1:2,sep="_")

VDJ.TRA.sub[VDJ.TRA.sub$subclone == subclones[1],]$subclone <- names(subclones)[1]
VDJ.TRA.sub[VDJ.TRA.sub$subclone == subclones[2],]$subclone <- names(subclones)[2]
#VDJ.TRA.sub[VDJ.TRA.sub$subclone == subclones[3],]$subclone <- names(subclones)[3]
VDJ.TRA.sub[!VDJ.TRA.sub$subclone %in% names(subclones),]$subclone <- "others"
table(VDJ.TRA.sub$subclone)

VDJ.TRA.sub$v_call_sim <- factor(VDJ.TRA.sub$v_call_sim,
                                 levels = names(sort(table(VDJ.TRA.sub$v_call_sim),
                                                     decreasing = T)))
VDJ.TRA.sub$j_call_sim <- factor(VDJ.TRA.sub$j_call_sim,
                                 levels = names(sort(table(VDJ.TRA.sub$j_call_sim),
                                                     decreasing = T)))

plot.data <- VDJ.TRA.sub %>%
  make_long(v_call_sim, j_call_sim)
plot.data$node <- factor(plot.data$node, levels = rev(c(levels(VDJ.TRA.sub$v_call_sim),
                                                        levels(VDJ.TRA.sub$j_call_sim))))

color.subclone <- c(brewer.pal(3,"Pastel2"),"gray75")
names(color.subclone) <- c(names(subclones),"others")
sort(table(plot.data$node),decreasing = T)

temp <- color.subclone[rep(c("Subclone_1", "Subclone_2", "Subclone_3","others"),time = c(3,3,3,4))]
names(temp) <-names(sort(table(plot.data$node),decreasing = T)) 

pdf("TCR/Sankey.subclone.H9.pdf",6,4)
ggplot(plot.data, aes(x = x, next_x = next_x,
                      node = node, next_node = next_node,
                      fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .8, smooth=8,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray30") +
  scale_x_discrete(breaks = c("v_call_sim","d_call_sim","j_call_sim"),
                   labels =c("TRAV", "TRAD", "TRAJ")) +
  # scale_fill_viridis_d(option = "D", alpha = .8)+
  scale_fill_manual(values = temp)+
  theme_sankey(base_size = 15) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("GM12878 subclones")
dev.off()

rm(p1,p2,p3,p4,p5,p,plot.data,obj.9CL,obj.9Mix,temp)

# 6. supplementary tables -------------------------------------------------
table(isoform.merge$Organism)
write.table(isoform.merge[isoform.merge$Organism == "Human",],
            "isoform.human.tsv", sep = "\t")
write.table(isoform.merge[isoform.merge$Organism == "Mouse",],
            "isoform.mouse.tsv", sep = "\t")

common <- intersect(colnames(VDJ.TRA),colnames(VDJ.TRB))
common <- intersect(common,colnames(VDJ.IGH))
common <- intersect(common,colnames(VDJ.IGL))

write.table(VDJ.TRA[,common],
            "VDJ.TRA.tsv", sep = "\t")
write.table(VDJ.TRB[,common],
            "VDJ.TRB.tsv", sep = "\t")
write.table(VDJ.IGH[,common],
            "VDJ.IGH.tsv", sep = "\t")
write.table(VDJ.IGL[,common],
            "VDJ.IGL.tsv", sep = "\t")
