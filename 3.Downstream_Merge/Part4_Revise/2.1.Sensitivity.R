setwd("E:/LabWork/Project/SCANSeq2/revise/2.1.Sensitivity/")

library(ggplot2)


# 1. load SCAN-seq2 data --------------------------------------------------
rc_gene_9CL <- readRDS("../../batch/1 9CL/data/rc_gene.RDS")
rc_gene_9Mix <- readRDS("../../batch/2 9CL_Mix/data/rc_gene.RDS")
rc_gene_4CL <- readRDS("../../batch/3 4CL/data/rc_gene.RDS")
rc_gene_100 <- readRDS("../../batch/4 UMI_100/data/rc_gene.RDS")
rc_gene_200 <- readRDS("../../batch/5 UMI_200/data/rc_gene.RDS")

sample.info.9CL <- read.table("../../batch/1 9CL/sample.info.txt")
sample.info.9Mix <-
  read.table("../../batch/2 9CL_Mix/sample.info.txt")
sample.info.4CL <- read.table("../../batch/3 4CL/sample.info.txt")
sample.info.100 <-
  read.table("../../batch/4 UMI_100/sample.info.txt")
sample.info.200 <-
  read.table("../../batch/5 UMI_200/sample.info.txt")
sample.info.merge <- rbind(
  sample.info.100,
  sample.info.200,
  sample.info.4CL,
  sample.info.9CL,
  sample.info.9Mix
)
sample.info.merge <-
  sample.info.merge[sample.info.merge$Pass_QC == 1, ]

rm(
  sample.info.100,
  sample.info.200,
  sample.info.4CL,
  sample.info.9CL,
  sample.info.9Mix
)


# 2. downsample Smart-seq3 data ---------------------------------------------
rc_Smart3 <-
  read.table("../../Test/Smart-Seq3/Smartseq3.HEK.fwdprimer.readcounts.txt")

DownSample <- function(x, prob = 0.1) {
  pool <- rep(rownames(rc_Smart3), times = x)
  samples <- sample(pool, size = as.integer(sum(x) * prob))
  out <- rep(0, nrow(rc_Smart3))
  names(out) <- rownames(rc_Smart3)
  samples <- table(samples)
  out[names(samples)] <- samples
  return(out)
}

rc_Smart3_3 <- apply(rc_Smart3, 2, DownSample, prob = 1 / 3)
rc_Smart3_5 <- apply(rc_Smart3, 2, DownSample, prob = 0.2)
rc_Smart3_10 <- apply(rc_Smart3, 2, DownSample, prob = 0.1)
rc_Smart3_20 <- apply(rc_Smart3, 2, DownSample, prob = 0.05)
rc_Smart3_50 <- apply(rc_Smart3, 2, DownSample, prob = 0.02)
rc_Smart3_100 <- apply(rc_Smart3, 2, DownSample, prob = 0.01)
rc_Smart3_500 <- apply(rc_Smart3, 2, DownSample, prob = 0.002)
rc_Smart3_1000 <- apply(rc_Smart3, 2, DownSample, prob = 0.001)

rc_Smart3_merge <- cbind(
  rc_Smart3_3,
  rc_Smart3_5,
  rc_Smart3_10,
  rc_Smart3_20,
  rc_Smart3_50,
  rc_Smart3_100,
  rc_Smart3_500,
  rc_Smart3_1000
)
rm(
  rc_Smart3_3,
  rc_Smart3_5,
  rc_Smart3_10,
  rc_Smart3_20,
  rc_Smart3_50,
  rc_Smart3_100,
  rc_Smart3_500,
  rc_Smart3_1000
)

sample.info.Smart3 <- data.frame(
  "cell" = colnames(rc_Smart3_merge),
  "Mapped_Reads" = colSums(rc_Smart3_merge),
  "Gene_Detected" =  colSums(rc_Smart3_merge >= 1)
)



# 3. compare methods ------------------------------------------------------
sample.info.merge$Method <- "SCAN-seq2"
sample.info.Smart3$Method <- "Smart-seq3"

temp <- intersect(colnames(sample.info.Smart3),
                  colnames(sample.info.merge))
plot.data <-
  rbind(sample.info.Smart3[, temp], sample.info.merge[, temp])

pdf("Sensititivy.compare.pdf", 5, 3)
ggplot(plot.data) +
  geom_smooth(
    aes(
      x = Mapped_Reads / 1e5,
      y = Gene_Detected / 1000,
      color = Method
    ),
    method = "loess",
    se = T,
    span = 0.4,
    size = 0.7
  ) +
  xlab("Mapped Reads (10^5)") +
  theme_classic()
dev.off()
