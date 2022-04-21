###Differential expression analysis of PBS samples only
#Jessy Slota
#2021-09-01

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

###What is different between PBS samples?
#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_CA1",]
smpl_dat <- smpl_dat[smpl_dat$Prions == "neg",]

#Add in ca1 enrichment scores
ca1 <- read.csv("CA1 allen gsea/Full_gsea_res.csv")
ca1 <- ca1[is.element(ca1$NAME, "FIELD_CA1") & is.element(ca1$Sample, rownames(smpl_dat)),]
smpl_dat <- smpl_dat[ca1$Sample,]
smpl_dat$ca1 <- ca1$NES
rm(ca1)
raw_dat <- raw_dat[,rownames(smpl_dat)]

#Assign samples into groups
smpl_dat$Grp <- "A"
smpl_dat[c("PBS_41_CA1", "PBS_37_CA1"),]$Grp <- "B"
smpl_dat[c("PBS_71_CA1", "PBS_72_CA1"),]$Grp <- "C"
smpl_dat[c("PBS_26_CA1", "PBS_27_CA1"),]$Grp <- "D"
smpl_dat[c("PBS_90_CA1", "PBS_38_CA1", "PBS_87_CA1"),]$Grp <- "E"

#Normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Grp)
dds <- DESeq(dds)

#PCA plot
norm_dat <- vst(dds)
pcaData <- plotPCA(norm_dat, intgroup=c("Grp", "Time", "ca1"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plt_cols <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(100))
ggplot(pcaData, aes(PC1, PC2, color=ca1, shape=Grp, label=name)) +
  geom_point(size=3) +
  #geom_label(size=2, nudge_y = 2, color="black") +
  scale_color_gradientn(colours = plt_cols) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave("plots/supplementary plots/PCA_plt_CA1_mock.svg", width = 4, height = 3)


#DE results
res <- results(dds, contrast = c("Grp", "A", "D"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,])

#DE results
res <- results(dds, contrast = c("Grp", "A", "E"))
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
write.table(rownames(res[res$baseMean>15 & res$padj < 0.05 & abs(res$log2FoldChange) > 0.85,]), "gene lists/poor_sample_de_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#DE results
res <- results(dds, contrast = c("Grp", "B", "D"))
res <- na.omit(as.data.frame(res[order(res$padj),]))
summary(res$padj < 0.05)
writeClipboard(rownames(res[res$log2FoldChange > 0.85 & res$padj < 0.05,]))
writeClipboard(rownames(res[res$log2FoldChange < -0.85 & res$padj < 0.05,]))
de_genes <- union(de_genes, rownames(res[res$baseMean>25 & res$padj < 0.01 & abs(res$log2FoldChange) > 1,]))

#DE results
res <- results(dds, contrast = c("Grp", "B", "E"))
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
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_CA1",]

#Add in ca1 enrichment scores
ca1 <- read.csv("CA1 allen gsea/Full_gsea_res.csv")
ca1 <- ca1[is.element(ca1$NAME, "FIELD_CA1") & is.element(ca1$Sample, rownames(smpl_dat)),]
smpl_dat <- smpl_dat[ca1$Sample,]
smpl_dat$ca1 <- ca1$NES
rm(ca1)
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
NES = circlize::colorRamp2(c(seq(min(smpls$ca1), -0.1, length.out=5), 0, seq(0.1, max(smpls$ca1), length.out=5)), rev(brewer.pal(11, "PuOr")))

gns <- brewer.pal(9, "Greens")
#add annotation
ha <- HeatmapAnnotation(Treatment=smpls$Treatment,
                        dpi=as.factor(smpls$dpi),
                        CA1_NES=smpls$ca1,
                        col=list(Treatment=c("Mock"="navy", "RML"="firebrick"),
                                 dpi=c(`70`=gns[2], `90`=gns[4], `130`=gns[7], `150`=gns[9]),
                                 CA1_NES=NES))

pdf("plots/supplementary plots/CA1_highly_variable_hmp.pdf", width=5, height = 7)
ht_list = Heatmap(hmap, col = z_scale, name="Z-score", show_column_names = FALSE, show_row_names = FALSE, top_annotation = ha)
draw(ht_list, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

