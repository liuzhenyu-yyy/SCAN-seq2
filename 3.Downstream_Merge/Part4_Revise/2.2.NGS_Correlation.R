library(ggplot2)
library(dplyr)

setwd("E:/LabWork/Project/SCANSeq2/revise/2.2.NGS_Correlation/")
source("E:/LabWork/code/MyFunction.R")

# 1. overall correlations -------------------------------------------------
plot.data <- readRDS("../../merge/Part1/NGS_Correlation.rds")
plot.data$Density <- get_density(log10(plot.data$Smart3_Reads+1),
                                 log10(plot.data$SCANSeq2+1),
                                 nbins = 1000)

cor(log10(plot.data$Smart3_Reads+1),log10(plot.data$Smart3_UMI+1)) # 0.9525138
cor(log10(plot.data$Smart3_Reads+1),log10(plot.data$SCANSeq2+1)) # 0.8309494
cor(log10(plot.data$Smart3_UMI+1),log10(plot.data$SCANSeq2+1)) # 0.8224188

model <- lm(log10(plot.data$SCANSeq2+1)~log10(plot.data$Smart3_Reads+1))
summary(model)

cor.value <- cor(log10(plot.data$Smart3_Reads+1),log10(plot.data$SCANSeq2+1))

pdf("Dot.Cor.Methods.pdf",5,4)
ggplot(plot.data)+
  geom_point(aes(x = log10(Smart3_Reads+1), y = log10(SCANSeq2+1), color = Density),size = 0.5)+
  geom_abline(slope = model$coefficients[2],
              intercept = model$coefficients[1],
              color = "blue")+
  xlab("log10(Smart-Seq3 Reads + 1)")+
  ylab("log10(SCAN-Seq2 UMIs + 1)")+
  theme_classic()+
  scale_color_viridis_c()+
  coord_equal(1)+
  ggtitle(paste0("Correlation: ", signif(cor.value, 2)))
dev.off()

# 3. each quantile --------------------------------------------------------
breaks <- quantile(plot.data$SCANSeq2,c(0,0.25,0.5,0.75,1))

p.list <- list()

for (one in 1:4) {
  plot.data.sub <- plot.data[plot.data$SCANSeq2 >= breaks[one] &
                               plot.data$SCANSeq2 < breaks[one + 1], ]
  
  model <- lm(log10(plot.data.sub$SCANSeq2+1) ~ log10(plot.data.sub$Smart3_Reads+1))
  
  cor.value <- cor(log10(plot.data.sub$Smart3_Reads+1),log10(plot.data.sub$SCANSeq2+1))
  
  p <- ggplot(plot.data.sub)+
    geom_point(aes(x = log10(Smart3_Reads+1), y = log10(SCANSeq2+1), color = Density),size = 0.5)+
    geom_abline(slope = model$coefficients[2], intercept = model$coefficients[1], color = "blue")+
    xlab("log10(Smart-Seq3 Reads + 1)")+
    ylab("log10(SCAN-Seq2 UMIs + 1)")+
    theme_classic()+
    scale_color_viridis_c()+
    xlim(c(0,6.5)) +ylim(c(0,6.5))+
    coord_equal(1)+
    ggtitle(paste0("Correlation: ", signif(cor.value, 2)))
  p.list[[one]] <- p
}

pdf("Dot.Cor.Split4.pdf", 8,8)
multiplot(plotlist = p.list, cols = 2)
dev.off()
