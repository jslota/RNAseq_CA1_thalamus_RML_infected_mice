###Differential expression analysis of PBS samples only
#Jessy Slota
#2021-09-01

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)

###What is different between PBS samples?
#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_TH",]
smpl_dat <- smpl_dat[smpl_dat$Prions == "neg",]

#Add in thalamus enrichment scores
th <- read.csv("thalamus allen gsea/Full_gsea_res.csv")
th <- th[is.element(th$NAME, "THALAMUS") & is.element(th$Sample, rownames(smpl_dat)),]
smpl_dat <- smpl_dat[th$Sample,]
smpl_dat$th <- th$NES
rm(th)
raw_dat <- raw_dat[,rownames(smpl_dat)]

#Assign samples into groups
smpl_dat$Grp <- "A"
smpl_dat[c("PBS_22_TH", "PBS_26_TH"),]$Grp <- "B"
smpl_dat[c("PBS_69_TH", "PBS_37_TH"),]$Grp <- "C"
smpl_dat[c("PBS_41_TH", "PBS_38_TH"),]$Grp <- "D"

#Normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Grp)
dds <- DESeq(dds)

#PCA plot
norm_dat <- vst(dds)
pcaData <- plotPCA(norm_dat, intgroup=c("Grp", "Time", "th"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plt_cols <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(100))
ggplot(pcaData, aes(PC1, PC2, color=th, shape=Grp, label=name)) +
  geom_point(size=3) +
  #geom_label(size=2, nudge_y = 2, color="black") +
  scale_color_gradientn(colours = plt_cols) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/supplementary plots/PCA_plot_thalamus_mock.svg", width = 4, height = 3)

#DE results
res <- results(dds, contrast = c("Grp", "A", "D"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,])

#DE results
res <- results(dds, contrast = c("Grp", "A", "B"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- union(de_genes, rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,]))

#DE results
res <- results(dds, contrast = c("Grp", "A", "C"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- union(de_genes, rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,]))

#DE results
res <- results(dds, contrast = c("Grp", "B", "D"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- union(de_genes, rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,]))

#DE results
res <- results(dds, contrast = c("Grp", "B", "C"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- union(de_genes, rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,]))

#DE results
res <- results(dds, contrast = c("Grp", "C", "D"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- union(de_genes, rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,]))


#Only get functional genes
de_genes <- de_genes[-grep("RP2*", de_genes)]
de_genes <- de_genes[-grep("AC1*", de_genes)]
de_genes <- de_genes[-grep("Gm", de_genes)]
writeClipboard(de_genes)

#heatmap all samples
#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_TH",]

#Add in thalamus enrichment scores
th <- read.csv("thalamus allen gsea/Full_gsea_res.csv")
th <- th[is.element(th$NAME, "THALAMUS") & is.element(th$Sample, rownames(smpl_dat)),]
smpl_dat <- smpl_dat[th$Sample,]
smpl_dat$th <- th$NES
rm(th)
raw_dat <- raw_dat[,rownames(smpl_dat)]

#Normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Prions)
dds <- DESeq(dds)
norm_dat <- vst(dds)

hmap <- assay(norm_dat)[de_genes,]
hmap <- (hmap - rowMeans(hmap))/rowSds(hmap)

smpls <- smpl_dat %>% 
  mutate(region=case_when(Experiment=="RML_CA1"~"CA1",
                          Experiment=="RML_TH"~"thalamus")) %>%
  mutate(Treatment=case_when(Prions=="neg"~"Mock",
                             Prions=="pos"~"RML"))

z_scale = circlize::colorRamp2(c(seq(min(hmap), -0.75, length.out=5), 0, seq(0.75, -min(hmap), length.out=5)), rev(brewer.pal(11, "PuOr")))
NES = circlize::colorRamp2(c(seq(min(smpls$th), -0.1, length.out=5), 0, seq(0.1, max(smpls$th), length.out=5)), rev(brewer.pal(11, "PuOr")))

gns <- brewer.pal(9, "Greens")
#add annotation
ha <- HeatmapAnnotation(Treatment=smpls$Treatment,
                        dpi=as.factor(smpls$dpi),
                        thalamus_NES=smpls$th,
                        col=list(Treatment=c("Mock"="navy", "RML"="firebrick"),
                                 dpi=c(`70`=gns[2], `90`=gns[4], `130`=gns[7], `150`=gns[9]),
                                 thalamus_NES=NES))

pdf("plots/supplementary plots/thalamus_highly_variable_hmp.pdf", width=5, height = 7)
ht_list = Heatmap(hmap, col = z_scale, name="Z-score", show_column_names = FALSE, show_row_names = FALSE, top_annotation = ha)
draw(ht_list, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()
