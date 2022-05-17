setwd("E:/LabWork/Project/SCANSeq2/merge/Part1/Pseudogene/")

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
obj.9CL <- readRDS("E:/LabWork/Project/SCANSeq2/batch/1 9CL/obj.9CL.rds")
obj.9Mix <- readRDS("E:/LabWork/Project/SCANSeq2/batch/2 9CL_Mix/obj.9CL_Mix.rds")
obj.4CL <- readRDS("E:/LabWork/Project/SCANSeq2/batch/3 4CL/obj.4CL.rds")

obj.combined <- merge(obj.9CL, y = c(obj.9Mix, obj.4CL),
                      add.cell.ids = c("9CL", "9Mix", "4CL"),
                      project = "SCAN_Seq2")
obj.combined$Organism <- ifelse(obj.combined$Cell_Line %in% c("GM12878","HepG2","Hela","H9","K562","293T"),
                                "Human","Mouse")
rm(obj.9CL, obj.9Mix, obj.4CL)

# color for cell lines
color.cell <- brewer.pal(11,"Set3")[c(1:8,10)]
names(color.cell) <- c("GM12878","HepG2","Hela","H9","K562","293T",
                       "AtT20","MEF","3T3")
color.cell["HepG2"] <- "#ffed6f"


# 0. primary Seurat pipeline: DR and cluster ------------------------------
# all cells
obj.combined <- obj.combined %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors()

ElbowPlot(obj.combined)
obj.combined <- RunUMAP(obj.combined, dims = 1:10,
                        n.neighbors = 40, min.dist = 0.1)

pdf("UMAP.Sample.info.pdf",8,3)
multiplot(DimPlot(obj.combined, group.by = "Cell_Line", cols = color.cell)+
            theme_lzy()+coord_equal(0.7),
            DimPlot(obj.combined, group.by = "Library")+
            theme_lzy()+coord_equal(0.7),cols = 2)
dev.off()

# subset human
table(obj.combined$Organism, obj.combined$Cell_Line)
obj.human <- subset(obj.combined, subset = Organism == "Human")
obj.human$Cell_Line <- factor(obj.human$Cell_Line,
                              levels = c("GM12878","HepG2","Hela","H9","K562","293T",
                                         "AtT20","MEF","3T3"))
Idents(obj.human) <- obj.human$Cell_Line
obj.human <- obj.human %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors()

ElbowPlot(obj.human)
table(obj.human$Cell_Line)
obj.human <- RunUMAP(obj.human, dims = 1:6,
                        n.neighbors = 60, min.dist = 0.3)

pdf("UMAP.Sample.info.human.pdf",8,3)
multiplot(DimPlot(obj.human, group.by = "Cell_Line", cols = color.cell)+
            theme_lzy()+coord_equal(0.8),
          DimPlot(obj.human, group.by = "Library")+
            theme_lzy()+coord_equal(0.8),cols = 2)
dev.off()

# calculate tpm
rc <- as.data.frame(obj.human@assays$Gene@counts) # 1538 human cells
tpm <- as.data.frame(apply(rc, 2, function(x){return(x * 1e6 / sum(x))}))
colSums(tpm)
rm(rc)

# 1. Distinguish Pseudogene and Parent ---------------------------------------

## 1.1  Distribution of similarity ------------------
similarity <- read.table("E:/LabWork/genome/pseudogene/similarity.dat.txt", header = T)
table(similarity$parent %in% trans.info.hg38$`Transcript stable ID version`) 
table(similarity$pgene %in% trans.info.hg38$`Transcript stable ID version`) 
similarity$parent <- gsub("\\..+?$","",similarity$parent)
similarity$pgene <- gsub("\\..+?$","",similarity$pgene)
table(similarity$parent %in% trans.info.hg38$`Transcript stable ID`) 
table(similarity$pgene %in% trans.info.hg38$`Transcript stable ID`) 

plot.data <- melt(similarity)
pdf("Density.similarity.pdf",3.8,3)
ggplot(plot.data)+
  geom_histogram(aes(x = value, fill = variable), alpha = 0.5, binwidth = 0.02)+
  theme_lzy()+facet_wrap(~variable,ncol = 1)+
  geom_vline(xintercept = 0.96, lty=2)+
  xlab("Similarity")+ylab("Number of Genes")
dev.off()

table(similarity$body<0.96) / nrow(similarity) # 0.92520872  
table(similarity$utr<0.96) / nrow(similarity) # 0.8155148 

## 1.2 performance at high similarity--------------
similarity.sub <- similarity
similarity.sub <- similarity[similarity$parent %in% trans.info.hg38$`Transcript stable ID`,]
similarity.sub <- similarity.sub[similarity.sub$pgene %in% trans.info.hg38$`Transcript stable ID`,]

temp <- trans.info.hg38$`Gene name`
names(temp) <- trans.info.hg38$`Transcript stable ID`
similarity.sub$parent_gene <- temp[similarity.sub$parent]
similarity.sub$pgene_gene <- temp[similarity.sub$pgene]
temp <- trans.info.hg38$`Transcript name`
names(temp) <- trans.info.hg38$`Transcript stable ID`
similarity.sub$parent_trans <- temp[similarity.sub$parent]
similarity.sub$pgene_trans <- temp[similarity.sub$pgene]

tpm.pgene <- tpm[unique(similarity.sub$pgene_gene),]
tpm.parent <- tpm[unique(similarity.sub$parent_gene),]

head(sort(rowMeans(tpm.parent),decreasing = T), n = 20)
head(sort(rowMeans(tpm.pgene),decreasing = T), n = 20)

features <- c("RPS2", "RPS2P46",# 0.9874858
              "EEF1A1","EEF1A1P19", # 0.9695431
              "GAPDH","GAPDHP1", # 0.9620758
              "EMB", "EMBP1"  # 0.9777035
)

rowMeans(tpm)[features]

pdf("Vln.Parent.Pgene.pdf",5,10)
p1 <- VlnPlot(obj.human, features = features[c(1,2)],
              cols = color.cell, pt.size = 0.01, ncol = 2,same.y.lims = T)
p2 <- VlnPlot(obj.human, features = features[c(3,4)],
              cols = color.cell, pt.size = 0.01, ncol = 2,same.y.lims = T)
p3 <- VlnPlot(obj.human, features = features[c(5,6)],
              cols = color.cell, pt.size = 0.01, ncol = 2,same.y.lims = T)
p4 <- VlnPlot(obj.human, features = features[c(7,8)],
              cols = color.cell, pt.size = 0.01, ncol = 2,same.y.lims = T)
multiplot(p1, p2, p3, p4, cols = 1)
dev.off()
rm(tpm.pgene, tpm.parent)

# pairwise alignment
library(spiralize)
p.align = grid.grabExpr({
  
  lines = readLines("alignment/RPS2-RPS2P46_copy.txt")
  query = lines[seq(1, length(lines), by = 4)]
  subject = lines[seq(3, length(lines), by = 4)]
  
  query = gsub("^\\S+\\s+\\S+\\s+|\\s+\\S+$", "", query)
  query = paste(query, collapse = "")
  query = strsplit(query, "")[[1]]
  
  subject = gsub("^\\S+\\s+\\S+\\s+|\\s+\\S+$", "", subject)
  subject = paste(subject, collapse = "")
  subject = strsplit(subject, "")[[1]]
  
  # # only first 500 base-pairs because of the space in the final combined figure
  # query = query[1:500]
  # subject = subject[1:500]
  
  query_alt = query
  query_alt[query_alt == subject] <- "."
  subject_alt = subject
  subject_alt[subject_alt == query] <- "."
  
  n = length(query)
  
  col.query = c("A" = 2, "T" = 4, "C" = 3, "G" = 7, "-" = "red", "." = NULL)
  col.subject = c("A" = 2, "T" = 4, "C" = 3, "G" = 7, "-" = "red", "." = NULL)
  
  spiral_initialize(xlim = c(0, n), start = 180, end = 360*3, scale_by = "curve")
  spiral_track(height = 0.6)
  spiral_text(1:n - 0.5, 0.3, query_alt, facing = "inside", 
              gp = gpar(fontsize = 8, col = col.query[query_alt]), just = "top",)
  spiral_text(1:n - 0.5, 0.7, subject_alt, facing = "inside", 
              gp = gpar(fontsize = 8, col = col.subject[subject_alt]), just = "bottom")
  # x = which(query == subject)
  # spiral_segments(x - 0.5, 0.35, x - 0.5, 0.65, gp = gpar(col = "#444444"))
  
  spiral_track(height = 0.15, background = FALSE)
  spiral_lines(TRACK_META$xlim, 0.5, gp = gpar(col = "gray40", lty = 2))
  
  spiral_lines(TRACK_META$xlim, -1.3, gp = gpar(col = "#87B6A7", lwd = 2.5),)
  spiral_lines(TRACK_META$xlim, -3.8, gp = gpar(col = "#AEADF0", lwd = 2.5))
  at = seq(0, n, by = 100)
  spiral_points(at - 0.5, 0.5, pch = 16, size = unit(3, "pt"))
  l = at %% 100 == 0
  spiral_text(at[l] - 0.5, 1, paste0(at[l], "bp"), just = "bottom",
              facing = "inside", nice_facing = TRUE, gp = gpar(fontsize = 8, col = "#444444"))
  
}, width = 20/3, height = 12/2)

pdf("Spiral.RPS2.pdf",4,4)
grid.draw(p.align)
dev.off()

# 2. global expression level of pesudogenes -------------------------------
# get pseudogene annotation
pseudogenes <- gene.info.hg38[grep("pseudogene",gene.info.hg38$`Gene type`),]
table(pseudogenes$`Gene type`)
pseudogenes <- pseudogenes[grep("IG_|TR_|rRNA_",pseudogenes$`Gene type`,invert = T),]
table(pseudogenes$`Gene type`)
pseudogenes$pseudo_type <- gsub("transcribed_|translated_","",pseudogenes$`Gene type`)
table(pseudogenes$pseudo_type)

View(pseudogenes[pseudogenes$pseudo_type=="pseudogene",])
pseudogenes <- pseudogenes[pseudogenes$pseudo_type != "pseudogene",] # 14501 ps
table(pseudogenes$pseudo_type)

table(pseudogenes$`Gene name` %in% rownames(obj.human))
table(duplicated(pseudogenes$`Gene name`))
pseudogenes <- pseudogenes[!duplicated(pseudogenes$`Gene name`),]

tpm.ps <- tpm[pseudogenes$`Gene name`,]

## 2.1. identify expressed pseudogenes -----------
breaks.cell <- c(3,5,10,20,50)
breaks.tpm <- c(0,1,10,100)
nGene.ps <- data.frame()

for (i in breaks.tpm) {
  for (j in breaks.cell) {
    nCells <- rowSums(tpm.ps > i)
    nCells <- aggregate(nCells, by = list("Type"=pseudogenes$pseudo_type),
                        FUN = function(x){return(sum(x > j))})
    colnames(nCells)[2] <- "Count"
    nCells$Cutoff <- i
    nCells$nCell <- j
    nGene.ps <- rbind(nGene.ps, nCells)
    rm(nCells)
  }
}

nGene.ps.1 <- nGene.ps[nGene.ps$Cutoff==breaks.tpm[1],]
nGene.ps.2 <- nGene.ps[nGene.ps$Cutoff==breaks.tpm[2],]
nGene.ps.3 <- nGene.ps[nGene.ps$Cutoff==breaks.tpm[3],]
nGene.ps.4 <- nGene.ps[nGene.ps$Cutoff==breaks.tpm[4],]
nGene.ps.1$Count <- nGene.ps.1$Count - nGene.ps.2$Count
nGene.ps.2$Count <- nGene.ps.2$Count - nGene.ps.3$Count
nGene.ps.3$Count <- nGene.ps.3$Count - nGene.ps.4$Count
nGene.ps <- rbind(nGene.ps.1,nGene.ps.2,nGene.ps.3,nGene.ps.4)
rm(nGene.ps.1,nGene.ps.2,nGene.ps.3,nGene.ps.4)
nGene.ps$Cutoff <- factor(nGene.ps$Cutoff, levels = breaks.tpm)
nGene.ps$nCell <- factor(nGene.ps$nCell, levels = breaks.cell)
nGene.ps$Type <- factor(nGene.ps$Type,
                        levels = c("processed_pseudogene","unprocessed_pseudogene",
                                   "unitary_pseudogene", "polymorphic_pseudogene"))
pdf("Bar.nGenes.pdf",6,3)
ggplot(nGene.ps)+
  geom_bar(aes(x = nCell, y = Count, fill = Cutoff),stat = "identity")+
  facet_wrap(~Type,scales = "free_y")+
  scale_fill_manual(name = "TPM Cutoff",
                    values = brewer.pal(4,"Set2")[c(4,2,1,3)])+
  theme_lzy()+
  xlab("Number of Cells")+ ylab("Number of Expressed Genes")
dev.off()

signif(c(742, 625, 66, 11) / c(10663, 3552, 235, 49)*100, 3) 

# expressed pseudogenes
identical(pseudogenes$`Gene name`,rownames(tpm.ps))
ps.expressed <- pseudogenes[rowSums(tpm.ps > 0) > 10,] # 1444
table(ps.expressed$pseudo_type)

identical(colnames(tpm.ps), colnames(obj.human))

temp <- apply(tpm.ps, 1, FUN = function(x){
  return(aggregate(x, list(obj.human$Cell_Line), mean)$x)
  })
temp <- t(as.data.frame(temp))
colnames(temp) <- levels(obj.human$Cell_Line)[1:6]
head(temp)
colSums(temp > 1)
mean(colSums(temp > 1))

## 2.2. expression level of pseudogenes -----------
table(obj.combined$Cell_Line)
ps.expressed$Mean_TPM <- apply(tpm.ps[ps.expressed$`Gene name`, ],
                               1,
                               function(x){return(mean(head(sort(x, decreasing = T),100)))})
ps.expressed$Expression_Level <- "none"
ps.expressed[ps.expressed$Mean_TPM >= 0,]$Expression_Level <- "0-10"
ps.expressed[ps.expressed$Mean_TPM >= 10,]$Expression_Level <- "10-50"
ps.expressed[ps.expressed$Mean_TPM >= 50,]$Expression_Level <- "50-100"
ps.expressed[ps.expressed$Mean_TPM >= 100,]$Expression_Level <- "100-500"
ps.expressed[ps.expressed$Mean_TPM >= 500,]$Expression_Level <- ">500"
table(ps.expressed$Expression_Level)
ps.expressed$Expression_Level <- factor(ps.expressed$Expression_Level,
                                        levels = rev(c("0-10","10-50","50-100","100-500",">500")))

plot.data <- as.data.frame(table(ps.expressed$pseudo_type, ps.expressed$Expression_Level))
colnames(plot.data) <- c("Type","TPM","Count")
plot.data$Type <- gsub("_pseudogene","",plot.data$Type)
plot.data$Type <- factor(plot.data$Type, levels = c("processed","unprocessed",
                                                    "unitary","polymorphic"))
plot.data$Percent <- plot.data$Count / rep(table(ps.expressed$pseudo_type)) *100

pdf("Bar.Expression.level.pdf",8,4)
multiplot(
  ggplot(plot.data)+
    geom_bar(aes(x = Type, y = Count, fill = TPM),stat = "identity")+
    theme_lzy()+theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    scale_fill_manual(values =rev(brewer.pal(6,"Reds")[1:5])),
  ggplot(plot.data)+
    geom_bar(aes(x = Type, y = Percent, fill = TPM),stat = "identity")+
    theme_lzy()+theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    scale_fill_manual(values =rev(brewer.pal(6,"Reds")[1:5])),
  cols=2)
dev.off()

temp <- ps.expressed$Mean_TPM
names(temp) <- ps.expressed$`Gene name`
head(sort(temp, decreasing = T),n = 10)

pdf("Vln.Pseudogene.Marker.pdf",10,2.5)
VlnPlot(object = obj.combined, c("AC092053.1","TPTEP1","C19orf48","GPX1"),
        cols = color.cell,ncol = 4)
dev.off()

# 3. pseudogene and parents -----------------------------------------------

## 3.1. identify co-expressed pairs ----------
ps.parent <- as.data.frame(fread("E:/LabWork/genome/pseudogene/pseudogene.parent.ensembl101.tsv",
                                 sep = "\t"))
rownames(ps.parent) <- ps.parent$`Gene stable ID`
ps.parent <- ps.parent[pseudogenes$`Gene stable ID`,]
pseudogenes$Parent <- ps.parent$parent
pseudogenes[pseudogenes$Parent == pseudogenes$`Gene name`,]$Parent <- "none"
rm(ps.parent)
length(unique(pseudogenes$Parent)) # 3256 - 1

# pgene
pseudogenes$Expressed_pgene <- ifelse(pseudogenes$`Gene name` %in% ps.expressed$`Gene name`,
                                      1,0)
table(pseudogenes$Expressed_pgene)

# parent
parents <- unique(pseudogenes$Parent) # 3256 - 1
table(parents %in% rownames(tpm))
parents <- parents[parents %in% rownames(tpm)] # 3246
parents <- parents[rowSums(tpm[parents,] > 0) > 10] # 2997

pseudogenes$Expressed_parent <- ifelse(pseudogenes$Parent %in% parents,
                                       1,0)
table(pseudogenes$Expressed_parent)

table(pseudogenes$Expressed_pgene,pseudogenes$Expressed_parent)

ps.coexpress <- pseudogenes[pseudogenes$Expressed_pgene ==1 & pseudogenes$Expressed_parent ==1,]
table(ps.coexpress$pseudo_type)

## 3.2. compare expression level of pgene and parent
ps.coexpress$Mean_TPM_Pgene <- apply(tpm.ps[ps.coexpress$`Gene name`, ],
                                     1,
                                     function(x){return(mean(head(sort(x, decreasing = T),100)))})
ps.coexpress$Mean_TPM_Parent <- apply(tpm[ps.coexpress$Parent, ],
                                      1,
                                      function(x){return(mean(head(sort(x, decreasing = T),100)))})

pdf("Bar.Expression.Compare.pdf",4,3.5)
ggplot(ps.coexpress[ps.coexpress$pseudo_type %in% c("processed_pseudogene","unprocessed_pseudogene"),])+
  geom_histogram(aes(x = log10(Mean_TPM_Pgene / Mean_TPM_Parent),
                     fill = pseudo_type), 
                 show.legend = F, alpha = 0.5,binwidth = 0.15)+
  facet_wrap(~pseudo_type,nrow = 3)+
  xlim(c(-4,4))+
  theme_lzy()+
  xlab(("log10(TPM_Pseudogene / TPM_Parent)"))
dev.off()

sum(ps.coexpress$Mean_TPM_Parent) / sum(ps.coexpress$Mean_TPM_Pgene)

## 3.2. correlation of pseudogene and parents------------
one <- ps.coexpress$`Gene stable ID`[1]
ps.coexpress$Spearman_cor <- 2
ps.coexpress$p.value <- 2

for (one in 1:nrow(ps.coexpress)) {
  test <- cor.test(x = as.numeric(tpm[ps.coexpress[one,]$`Gene name`,]),
                   y = as.numeric(tpm[ps.coexpress[one,]$Parent,]),
                   method = "s",exact = F)
  ps.coexpress[one,]$Spearman_cor <- test$estimate
  ps.coexpress[one,]$p.value <- test$p.value
}

ps.coexpress$p.adjust <- p.adjust(p = ps.coexpress$p.value,
                                  method = "BH",
                                  n = nrow(ps.coexpress))
table(ps.coexpress$p.adjust < 0.01)

ps.coexpress$Correltaion <- "not significant"
ps.coexpress[ps.coexpress$Spearman_cor > 0.2,]$Correltaion <- "Positive"
ps.coexpress[ps.coexpress$Spearman_cor < (-0.2),]$Correltaion <- "Negative"
table(ps.coexpress$Correltaion)

ps.coexpress$label <- paste(ps.coexpress$Parent, ps.coexpress$`Gene name`, sep = "_")
ps.coexpress[ps.coexpress$Correltaion == "not significant",]$label <- NA

pdf("Dot.Parent.Correlation.pdf",4.5,4)
ggplot(ps.coexpress,aes(x = Spearman_cor, y = -log10(p.adjust)))+
  geom_point(aes(color = Correltaion))+
  ggrepel::geom_label_repel(aes(label = label,color = Correltaion),
                            box.padding = 0.1,
                            label.padding = 0.1,
                            label.r = 0.05,
                            label.size = 0.1,
                            min.segment.length = 10,
                            max.overlaps = 10,
                            size = 2,
                   segment.color = 'grey50') +
  scale_color_manual(values = c("#619cff","gray75","#f8766d"))+
  theme_classic()+
  ylim(c(0,150))+
  xlim(c(-0.6,0.6))
dev.off()

# dot plot for each cell
ps.pos <- ps.coexpress[ps.coexpress$Correltaion == "Positive",]$`Gene stable ID`
ps.pos <- ps.pos[order(ps.coexpress[ps.pos,]$Spearman_cor)]
ps.neg <- ps.coexpress[ps.coexpress$Correltaion == "Negative",]$`Gene stable ID`
ps.neg <- ps.neg[order(ps.coexpress[ps.neg,]$Spearman_cor)]

p <- list()
for (one in ps.pos) {
  plot.data <- data.frame("Parent"=as.numeric(tpm[ps.coexpress[one,]$Parent,]),
                          "Pgene"=as.numeric(tpm[ps.coexpress[one,]$`Gene name`,]),
                          "Cell_Line"=obj.human@meta.data[colnames(tpm),]$Cell_Line)
  p1 <- ggplot(plot.data)+
    geom_point(aes(x = Parent,
                   y = Pgene,
                   color = Cell_Line),show.legend = F)+
    theme_lzy()+
    scale_color_manual(values = color.cell)+
    xlab(ps.coexpress[one,]$Parent)+
    ylab(ps.coexpress[one,]$`Gene name`)+
    ggtitle(paste("Cor = ", signif(ps.coexpress[one,]$Spearman_cor,2)))
  p[[one]] <- p1
}

pdf("Dot.Positive.Correlation.pdf",9,9)
multiplot(plotlist = p, cols = 3)
dev.off()

p <- list()
for (one in ps.neg) {
  plot.data <- data.frame("Parent"=as.numeric(tpm[ps.coexpress[one,]$Parent,]),
                          "Pgene"=as.numeric(tpm[ps.coexpress[one,]$`Gene name`,]),
                          "Cell_Line"=obj.human@meta.data[colnames(tpm),]$Cell_Line)
  p1 <- ggplot(plot.data)+
    geom_point(aes(x = Parent,
                   y = Pgene,
                   color = Cell_Line),show.legend = F)+
    theme_lzy()+
    scale_color_manual(values = color.cell)+
    xlab(ps.coexpress[one,]$Parent)+
    ylab(ps.coexpress[one,]$`Gene name`)+
    ggtitle(paste("Cor = ", signif(ps.coexpress[one,]$Spearman_cor,2)))
  p[[one]] <- p1
}

pdf("Dot.Negative.Correlation.pdf",6,6)
multiplot(plotlist = p, cols = 2)
dev.off()

rm(plot.data, p1,p)

# 4. dr and cluster on pseudogenes ----------------------------------------
rc <- as.data.frame(obj.human@assays$Gene@counts)
rc_ps <- rc[pseudogenes$`Gene name`,]
rm(rc)

obj.human.ps <- CreateSeuratObject(rc_ps,meta.data = obj.human@meta.data)
obj.human.ps <- obj.human.ps %>% NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
ElbowPlot(obj.human.ps)

obj.human.ps <- RunUMAP(obj.human.ps, dims = 1:9,
                        n.neighbors = 40, min.dist = 0.1)

obj.human.ps <- FindNeighbors(obj.human.ps)
obj.human.ps <- FindClusters(obj.human.ps,resolution = 0.4)

pdf("PS.UMAP.Sample.info.pdf",8,4.5)
multiplot(DimPlot(obj.human.ps, group.by = "Cell_Line", cols = color.cell,label = T)+
            theme_lzy()+coord_equal(0.9)+NoLegend(),
          DimPlot(obj.human.ps, group.by = "RNA_snn_res.0.4",label = T,cols = unname(color.cell))+
            theme_lzy()+coord_equal(0.9)+NoLegend(),cols = 2)
dev.off()

table(obj.human.ps$Cell_Line, obj.human.ps$RNA_snn_res.0.4)
Idents(obj.human.ps) <- obj.human.ps$Cell_Line

markers.ps <- FindAllMarkers(obj.human.ps, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)

table(markers.ps$cluster)

top10 <- markers.ps %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
table(top10$cluster)

pdf("PS.Heatmap.Marker.all.pdf",6,5)
obj.human.ps <- ScaleData(obj.human.ps,features =  top10$gene)
DoHeatmap(obj.human.ps, features = top10$gene,
          draw.lines = F,
          group.colors = color.cell)+
  scale_fill_gradientn(colors = rev(brewer.pal(7,"RdYlBu")))
obj.human.ps <- ScaleData(obj.human.ps)
dev.off()


# 5. supplementary tables -------------------------------------------------
write.table(pseudogenes,"pseudogenes.tsv",sep ="\t")
write.table(ps.expressed,"ps.expressed.tsv",sep ="\t")
write.table(ps.coexpress,"ps.coexpress.tsv",sep ="\t")
write.table(markers.ps,"markers.ps.tsv", sep ="\t")
