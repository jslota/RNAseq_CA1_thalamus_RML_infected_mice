###Visualize allen gsea
#2021-08-26
#Jessy Slota

library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(reshape2)
library(ggplot2)

#load results
res <- read.csv("CA1 allen gsea/Full_gsea_res.csv")

#Use best selected gene sets
names <- read.delim("CA1 allen gsea/Best gene sets.txt", header = FALSE)$V1
hmp_dt <- reshape(res[,c(1,4,8)], idvar = "NAME", timevar = "Sample", direction = "wide", new.row.names = seq(1, length(unique(res$NAME))))
rownames(hmp_dt) <- hmp_dt$NAME
hmp_dt <- hmp_dt[names ,-1]
colnames(hmp_dt) <- gsub("NES.", "", colnames(hmp_dt))

plt_cols <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(100))
brks <- c(seq(min(na.omit(hmp_dt)), -0.75, length=50), seq(0.75, max(na.omit(hmp_dt)), length=50))

#plot heatmap
pheatmap(hmp_dt, plt_cols, breaks = brks, treeheight_row = 15, treeheight_col = 15, border_color = NA)

#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_CA1",]
raw_dat <- raw_dat[,rownames(smpl_dat)]

#Add in gsea data
smpl_dat$CA1 <- as.numeric(hmp_dt["FIELD_CA1",rownames(smpl_dat)])
smpl_dat$CA2 <- as.numeric(hmp_dt["FIELD_CA2_STRATUM_ORIENS",rownames(smpl_dat)])
smpl_dat$CA3 <- as.numeric(hmp_dt["FIELD_CA3_STRATUM_ORIENS",rownames(smpl_dat)])
smpl_dat$VAP <- as.numeric(hmp_dt["MANTLE_ZONE_OF_VAP",rownames(smpl_dat)])
smpl_dat$Thalamus <- as.numeric(hmp_dt["THALAMUS",rownames(smpl_dat)])
smpl_dat$P1B <- as.numeric(hmp_dt["SUPERFICIAL_STRATUM_OF_P1B",rownames(smpl_dat)])

#normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Prions)
dds <- DESeq(dds)
norm_dat <- vst(dds)

#make PCA
pcaData <- plotPCA(norm_dat, intgroup=c("Prions", "Time","CA1", "CA2", "CA3", "VAP", "P1B", "Thalamus" ), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData <- melt(pcaData, id.vars = c("PC1", "PC2", "group", "Prions", "Time", "name"))

ggplot(pcaData, aes(PC1, PC2, color=value, shape=Prions)) +
  geom_point(size=2) +
  scale_color_gradientn(colours = plt_cols) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  facet_wrap(~variable) +
  theme_bw() +
  theme(panel.grid = element_blank())

unique(pcaData[pcaData$PC1 >0 &pcaData$PC2 > -10,]$name)

#make PCA CA1 only
pcaData <- plotPCA(norm_dat, intgroup=c("Prions", "Time","CA1"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pcaData <- melt(pcaData, id.vars = c("PC1", "PC2", "group", "Prions", "Time", "name"))

ggplot(pcaData, aes(PC1, PC2, color=CA1, shape=Prions)) +
  geom_point(size=2) +
  scale_color_gradientn(colours = plt_cols) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/supplementary plots/PCA_plot_CA1_Allen_ES.svg", width = 3.5, height = 2.5)
