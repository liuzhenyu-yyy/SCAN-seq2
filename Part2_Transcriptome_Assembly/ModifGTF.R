setwd("E:/LabWork/Project/SCANSeq2/merge/Part2/data/")

# 1.IGH and IGL GTF -------------------------------------------------------
BCR.gtf <- read.table("GTF/BCR.gtf",sep = "\t")

IGH.gtf <- BCR.gtf[grep("IGH",BCR.gtf$V9),]
IGH.gtf$symbols <- unlist(lapply(strsplit(IGH.gtf$V9,";"),
                                 function(x){return(gsub(" gene_name ", "",
                                                         grep("gene_name",x,value = T)))}))
table(IGH.gtf$symbols)
IGH.gtf <- IGH.gtf[IGH.gtf$V1 == 14,]
IGH.gtf <- IGH.gtf[grep("P",IGH.gtf$symbols, invert  = T),]
IGH.gtf$Type <- "none"
IGH.gtf[grep("IGHV",IGH.gtf$symbols),]$Type <- "V"
IGH.gtf[grep("IGHD",IGH.gtf$symbols),]$Type <- "D"
IGH.gtf[grep("IGHJ",IGH.gtf$symbols),]$Type <- "J"
IGH.gtf[grep("IGHM|IGHA|IGHG|IGHE",IGH.gtf$symbols),]$Type <- "Constant"
IGH.gtf[IGH.gtf$symbols=="IGHD",]$Type <- "Constant"
table(IGH.gtf$Type)

color.vdj <- c("#fed995","#c4574f","#0079bb","#9086b9")
names(color.vdj) <- c("V","D","J","Constant")
IGH.gtf$Color <- color.vdj[IGH.gtf$Type]
IGH.gtf$V9 <- paste0(IGH.gtf$V9," type ", IGH.gtf$Type, "; color ", IGH.gtf$Color,";")

write.table(IGH.gtf,"GTF/Hg38.IGH.gtf",row.names = F, col.names = F, quote = F,sep = "\t")

IGL.gtf <- BCR.gtf[grep("IGL",BCR.gtf$V9),]
IGL.gtf$symbols <- unlist(lapply(strsplit(IGL.gtf$V9,";"),
                                 function(x){return(gsub(" gene_name ", "",
                                                         grep("gene_name",x,value = T)))}))
table(IGL.gtf$symbols)
IGL.gtf <- IGL.gtf[IGL.gtf$V1 == 22,]
IGL.gtf <- IGL.gtf[grep("P",IGL.gtf$symbols, invert  = T),]
IGL.gtf <- IGL.gtf[grep("IGLL",IGL.gtf$symbols, invert  = T),]
IGL.gtf$Type <- "none"
IGL.gtf[grep("IGLV",IGL.gtf$symbols),]$Type <- "V"
# IGL.gtf[grep("IGLD",IGL.gtf$symbols),]$Type <- "D"
IGL.gtf[grep("IGLJ",IGL.gtf$symbols),]$Type <- "J"
IGL.gtf[grep("IGLC",IGL.gtf$symbols),]$Type <- "Constant"
table(IGL.gtf$Type)

IGL.gtf$Color <- color.vdj[IGL.gtf$Type]
IGL.gtf$V9 <- paste0(IGL.gtf$V9," type ", IGL.gtf$Type, "; color ", IGL.gtf$Color,";")

write.table(IGL.gtf,"GTF/Hg38.IGL.gtf",row.names = F, col.names = F, quote = F,sep = "\t")

IGK.gtf <- BCR.gtf[grep("IGK",BCR.gtf$V9),]
IGK.gtf$symbols <- unlist(lapply(strsplit(IGK.gtf$V9,";"),
                                 function(x){return(gsub(" gene_name ", "",
                                                         grep("gene_name",x,value = T)))}))
table(IGK.gtf$symbols)
IGK.gtf <- IGK.gtf[IGK.gtf$V1 == 2,]
IGK.gtf <- IGK.gtf[grep("P",IGK.gtf$symbols, invert  = T),]
IGK.gtf$Type <- "none"
IGK.gtf[grep("IGKV",IGK.gtf$symbols),]$Type <- "V"
IGK.gtf[grep("IGKD",IGK.gtf$symbols),]$Type <- "D"
IGK.gtf[grep("IGKJ",IGK.gtf$symbols),]$Type <- "J"
IGK.gtf[grep("IGKC",IGK.gtf$symbols),]$Type <- "Constant"
IGK.gtf[IGK.gtf$symbols=="IGKD",]$Type <- "Constant"
table(IGK.gtf$Type)

color.vdj <- c("#fed995","#c4574f","#0079bb","#9086b9")
names(color.vdj) <- c("V","D","J","Constant")
IGK.gtf$Color <- color.vdj[IGK.gtf$Type]
IGK.gtf$V9 <- paste0(IGK.gtf$V9," type ", IGK.gtf$Type, "; color ", IGK.gtf$Color,";")

write.table(IGK.gtf,"GTF/Hg38.IGK.gtf",row.names = F, col.names = F, quote = F,sep = "\t")


# 2.TRB and TRA GTF -------------------------------------------------------
TCR.gtf <- read.table("GTF/TCR.gtf",sep = "\t")


TRB.gtf <- TCR.gtf[grep("TRB",TCR.gtf$V9),]
TRB.gtf$symbols <- unlist(lapply(strsplit(TRB.gtf$V9,";"),
                                 function(x){return(gsub(" gene_name ", "",
                                                         grep("gene_name",x,value = T)))}))
table(TRB.gtf$symbols)
TRB.gtf <- TRB.gtf[TRB.gtf$V1 == 7,]
TRB.gtf <- TRB.gtf[grep("P",TRB.gtf$symbols, invert  = T),]
TRB.gtf$Type <- "none"
TRB.gtf[grep("TRBV",TRB.gtf$symbols),]$Type <- "V"
TRB.gtf[grep("TRBD",TRB.gtf$symbols),]$Type <- "D"
TRB.gtf[grep("TRBJ",TRB.gtf$symbols),]$Type <- "J"
TRB.gtf[grep("TRBC",TRB.gtf$symbols),]$Type <- "Constant"
table(TRB.gtf$Type)

color.vdj <- c("#fed995","#c4574f","#0079bb","#9086b9")
names(color.vdj) <- c("V","D","J","Constant")
TRB.gtf$Color <- color.vdj[TRB.gtf$Type]
TRB.gtf$V9 <- paste0(TRB.gtf$V9," type ", TRB.gtf$Type, "; color ", TRB.gtf$Color,";")

write.table(TRB.gtf,"GTF/Hg38.TRB.gtf",row.names = F, col.names = F, quote = F,sep = "\t")

TRA.gtf <- TCR.gtf[grep("TRA",TCR.gtf$V9),]
TRA.gtf$symbols <- unlist(lapply(strsplit(TRA.gtf$V9,";"),
                                 function(x){return(gsub(" gene_name ", "",
                                                         grep("gene_name",x,value = T)))}))
table(TRA.gtf$symbols)
TRA.gtf <- TRA.gtf[TRA.gtf$V1 == 14,]
TRA.gtf <- TRA.gtf[grep("P",TRA.gtf$symbols, invert  = T),]
TRA.gtf <- TRA.gtf[grep("TRAL",TRA.gtf$symbols, invert  = T),]
TRA.gtf$Type <- "none"
TRA.gtf[grep("TRAV",TRA.gtf$symbols),]$Type <- "V"
# TRA.gtf[grep("TRAD",TRA.gtf$symbols),]$Type <- "D"
TRA.gtf[grep("TRAJ",TRA.gtf$symbols),]$Type <- "J"
TRA.gtf[grep("TRAC",TRA.gtf$symbols),]$Type <- "Constant"
table(TRA.gtf$Type)

TRA.gtf$Color <- color.vdj[TRA.gtf$Type]
TRA.gtf$V9 <- paste0(TRA.gtf$V9," type ", TRA.gtf$Type, "; color ", TRA.gtf$Color,";")

write.table(TRA.gtf,"GTF/Hg38.TRA.gtf",row.names = F, col.names = F, quote = F,sep = "\t")

#### merge count
merge.gtf <- rbind(IGH.gtf, IGL.gtf, TRA.gtf, TRB.gtf)
merge.gtf <- merge.gtf[merge.gtf$V3=="gene",]
table(duplicated(merge.gtf$symbols))

merge.gtf$Locus <- substr(merge.gtf$symbols, 1,3)
table(merge.gtf$Locus, merge.gtf$Type)
