###DE analysis
#2021-08-19
#Jessy Slota

library(DESeq2)
library(ggplot2)

#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_TH",]

smpls <- c("RML_118_TH", "RML_120_TH", "RML_121_TH", "RML_134_TH", "RML_135_TH", "RML_138_TH", 
           "RML_165_TH", "RML_166_TH", "RML_167_TH", "RML_185_TH", "RML_186_TH", "RML_195_TH",
          "PBS_24_TH", "PBS_26_TH", "PBS_38_TH",
          "PBS_71_TH", "PBS_72_TH", "PBS_87_TH", "PBS_88_TH", "PBS_90_TH")
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
ggsave("plots/PCA_plot_thalamus_smpls.svg", width = 3.5, height = 2.5)

#70 dpi
smpls <- c("RML_118_TH", "RML_120_TH", "RML_121_TH",
           "PBS_22_TH", "PBS_24_TH", "PBS_26_TH", "PBS_71_TH", "PBS_90_TH")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/TH_DE_70_dpi.csv")

#90 dpi
smpls <- c("RML_134_TH", "RML_135_TH", "RML_138_TH",
           "PBS_37_TH", "PBS_38_TH", "PBS_41_TH", "PBS_71_TH", "PBS_72_TH", "PBS_87_TH", "PBS_88_TH", "PBS_90_TH")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/TH_DE_90_dpi.csv")


###130 dpi
smpls <- c("RML_165_TH", "RML_166_TH", "RML_167_TH",
          "PBS_22_TH", "PBS_24_TH", "PBS_26_TH", "PBS_41_TH", "PBS_69_TH",
          "PBS_71_TH", "PBS_72_TH", "PBS_88_TH", "PBS_90_TH")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/TH_DE_130_dpi.csv")

###150 dpi
smpls <- c("RML_185_TH", "RML_186_TH", "RML_195_TH",
           "PBS_71_TH", "PBS_72_TH", "PBS_87_TH", "PBS_88_TH", "PBS_90_TH")
dds <- DESeqDataSetFromMatrix(raw_dat[,smpls], smpl_dat[smpls,], ~Prions)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Prions", "pos", "neg"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
write.csv(res, "results/TH_DE_150_dpi.csv")
