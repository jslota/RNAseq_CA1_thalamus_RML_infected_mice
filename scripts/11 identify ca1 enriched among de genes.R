library(dplyr)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(enrichR)

#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
raw_dat <- raw_dat[,rownames(smpl_dat)]

#normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Experiment+Treatment+Time)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Experiment", "RML_CA1", "RML_TH"))
res <- na.omit(as.data.frame(res))
write.csv(res, "gene lists/DE_res_CA1_vs_TH.csv")
ca1_enriched_genes <- res %>% filter(baseMean > 15 & padj < 0.01 & log2FoldChange > 1) %>% rownames()
th_enriched_genes <- res %>% filter(baseMean > 15 & padj < 0.01 & log2FoldChange < -1) %>% rownames()

gene_list <- read.delim("cell type databases/final_cell_type_list.txt")
ca1_enriched_genes <- intersect(ca1_enriched_genes, gene_list$ms)
th_enriched_genes <- intersect(th_enriched_genes, gene_list$ms)

write.table(ca1_enriched_genes, "gene lists/ca1_enriched_compared_to_th_genes.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(th_enriched_genes, "gene lists/th_enriched_compared_to_ca1_genes.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
