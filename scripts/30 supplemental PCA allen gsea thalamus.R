###Visualize allen gsea
#2021-08-25
#Jessy Slota

library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(reshape2)
library(ggplot2)

#load results
res <- read.csv("thalamus allen gsea/Full_gsea_res.csv")

#Use best selected gene sets
names <- read.delim("thalamus allen gsea/Best gene sets.txt", header = FALSE)$V1
hmp_dt <- reshape(res[,c(1,4,8)], idvar = "NAME", timevar = "Sample", direction = "wide", new.row.names = seq(1, length(unique(res$NAME))))
rownames(hmp_dt) <- hmp_dt$NAME
hmp_dt <- hmp_dt[names ,-1]
colnames(hmp_dt) <- gsub("NES.", "", colnames(hmp_dt))

plt_cols <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(100))
brks <- c(seq(min(na.omit(hmp_dt)), -0.75, length=50), seq(0.75, max(na.omit(hmp_dt)), length=50))

#plot heatmap
pheatmap(hmp_dt, plt_cols, breaks = brks, treeheight_row = 15, treeheight_col = 15)

#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_TH",]
raw_dat <- raw_dat[,rownames(smpl_dat)]

#Add in gsea data
smpl_dat$Thalamus <- as.numeric(hmp_dt["THALAMUS",rownames(smpl_dat)])
smpl_dat$CA2 <- as.numeric(hmp_dt["FIELD_CA2",rownames(smpl_dat)])
smpl_dat$DentateGyrus <- as.numeric(hmp_dt["DENTATE_GYRUS",rownames(smpl_dat)])
smpl_dat$PreOptic <- as.numeric(hmp_dt["PREOPTIC_AREA",rownames(smpl_dat)])
smpl_dat$SupraOptic <- as.numeric(hmp_dt["SUPRAOPTIC_NUCLEUS",rownames(smpl_dat)])
smpl_dat$ZonaIncerta <- as.numeric(hmp_dt["ZONA_INCERTA",rownames(smpl_dat)])

#normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Prions)
dds <- DESeq(dds)
norm_dat <- vst(dds)

#make PCA
pcaData <- plotPCA(norm_dat, intgroup=c("Prions", "Time", "Thalamus", "CA2", "DentateGyrus", "PreOptic", "SupraOptic", "ZonaIncerta"), returnData=TRUE)
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

smpl_dat[smpl_dat$Thalamus > 1 & smpl_dat$Dentate_Gyrus < 0 & smpl_dat$SupraOptic < 0,]

#make PCA TH only
pcaData <- plotPCA(norm_dat, intgroup=c("Prions", "Time","Thalamus"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pcaData <- melt(pcaData, id.vars = c("PC1", "PC2", "group", "Prions", "Time", "name"))

ggplot(pcaData, aes(PC1, PC2, color=Thalamus, shape=Prions)) +
  geom_point(size=2) +
  scale_color_gradientn(colours = plt_cols) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/supplementary plots/PCA_plot_thalamus_Allen_ES.svg", width = 3.5, height = 2.5)
