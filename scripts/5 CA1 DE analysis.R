###DE analysis
#2021-08-26
#Jessy Slota

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_CA1",]
raw_dat <- raw_dat[,rownames(smpl_dat)]

#Best samples
smpls <- c("RML_118_CA1", "RML_120_CA1", "RML_122_CA1", "RML_121_CA1", "RML_134_CA1", "RML_135_CA1", "RML_138_CA1", "RML_165_CA1", "RML_166_CA1", "RML_167_CA1", "RML_185_CA1", "RML_186_CA1", "RML_195_CA1",
           "PBS_22_CA1", "PBS_24_CA1", "PBS_26_CA1", "PBS_27_CA1", "PBS_71_CA1", "PBS_72_CA1", "PBS_88_CA1")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)

#PCA plot
norm_dat <- vst(dds)
pcaData <- plotPCA(norm_dat, intgroup=c("Prions", "Time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Prions, shape=Time, label=name)) +
  geom_point(size=2) +
  geom_label(size=2, nudge_y = 2) +
  scale_color_manual(values = c("Firebrick", "Navy")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())

#make plot for figure
ggplot(pcaData, aes(PC1, PC2, color=Prions, shape=Time)) +
  geom_point(size=2) +
  scale_color_manual(values = c("Navy", "Firebrick")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/PCA_plot_CA1_smpls.svg", width = 3.5, height = 2.5)

#70 dpi
smpls <- c("RML_120_CA1", "RML_122_CA1","RML_118_CA1", "RML_121_CA1",
           "PBS_22_CA1", "PBS_24_CA1", "PBS_88_CA1", "PBS_26_CA1", "PBS_27_CA1")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/CA1_DE_70_dpi.csv")

#90 dpi
smpls <- c("RML_134_CA1", "RML_135_CA1", "RML_138_CA1",
           "PBS_71_CA1", "PBS_72_CA1")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/CA1_DE_90_dpi.csv")

#130 dpi
smpls <- c("RML_165_CA1", "RML_166_CA1", "RML_167_CA1",
           "PBS_71_CA1", "PBS_72_CA1")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/CA1_DE_130_dpi.csv")


#150 dpi
smpls <- c("RML_185_CA1", "RML_186_CA1", "RML_195_CA1",
           "PBS_22_CA1", "PBS_24_CA1", "PBS_26_CA1", "PBS_27_CA1", "PBS_71_CA1", "PBS_72_CA1", "PBS_88_CA1")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/CA1_DE_150_dpi.csv")

