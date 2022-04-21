###Comparing DE res RML CA1 and TH
#2021-09-08
#Jessy Slota

library(dplyr)

#load results
#CA1
ca1_70 <- read.csv("results/CA1_DE_70_dpi.csv")
ca1_90 <- read.csv("results/CA1_DE_90_dpi.csv")
ca1_130 <- read.csv("results/CA1_DE_130_dpi.csv")
ca1_150 <- read.csv("results/CA1_DE_150_dpi.csv")

#TH
th_70 <- read.csv("results/TH_DE_70_dpi.csv")
th_90 <- read.csv("results/TH_DE_90_dpi.csv")
th_130 <- read.csv("results/TH_DE_130_dpi.csv")
th_150 <- read.csv("results/TH_DE_150_dpi.csv")

ca1_genes_up <- ca1_150 %>% filter(baseMean > 15 & padj < 0.05 & log2FoldChange > 0.5) %>% pull(X)

ca1_genes_dn <- ca1_150 %>% filter(baseMean > 15 & padj < 0.05 & log2FoldChange < -0.5) %>% pull(X)

#n DE genes
ca1_150 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.5) %>% nrow()
th_130 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.5) %>% nrow()
th_150 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.5) %>% nrow()

#Best overall marker genes
th_genes_up <- union(th_150 %>% filter(baseMean > 15 & padj < 0.05 & log2FoldChange > 0.5) %>% pull(X),
                     th_130 %>% filter(baseMean > 15 & padj < 0.05 & log2FoldChange > 0.5) %>% pull(X))

th_genes_dn <- union(th_150 %>% filter(baseMean > 15 & padj < 0.05 & log2FoldChange < -0.5) %>% pull(X),
                     th_130 %>% filter(baseMean > 15 & padj < 0.05 & log2FoldChange < -0.5) %>% pull(X))


#filter out poor sample quality genes from results
tmp <- read.delim("gene lists/poor_sample_de_genes.txt", header = FALSE)$V1
ca1_genes_up <- setdiff(ca1_genes_up, tmp)
ca1_genes_dn <- setdiff(ca1_genes_dn, tmp)
th_genes_up <- setdiff(th_genes_up, tmp)
th_genes_dn <- setdiff(th_genes_dn, tmp)
rm(tmp)

#get common and all up and down genes
common_up <- intersect(ca1_genes_up, th_genes_up)
all_up <- union(ca1_genes_up, th_genes_up)
common_dn <- intersect(ca1_genes_dn, th_genes_dn)
all_dn <- union(ca1_genes_dn, th_genes_dn)

#remaining genes
remain_genes <- ca1_70 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.85) %>% pull(X)
remain_genes <- union(remain_genes, ca1_90 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.85) %>% pull(X))
remain_genes <- union(remain_genes, ca1_130 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.85) %>% pull(X))
remain_genes <- union(remain_genes, th_70 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.85) %>% pull(X))
remain_genes <- union(remain_genes, th_90 %>% filter(baseMean > 15 & padj < 0.05 & abs(log2FoldChange) > 0.85) %>% pull(X))

remain_genes <- setdiff(remain_genes, all_up)
remain_genes <- setdiff(remain_genes, all_dn)
write.table(remain_genes, "gene lists/de_genes_other_timepoints.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#save lists
write.table(ca1_genes_up, "gene lists/ca1_genes_up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(ca1_genes_dn, "gene lists/ca1_genes_dn.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(th_genes_up, "gene lists/th_genes_up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(th_genes_dn, "gene lists/th_genes_dn.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(common_up, "gene lists/common_up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(common_dn, "gene lists/common_dn.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(all_up, "gene lists/all_up.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(all_dn, "gene lists/all_dn.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#get genes for SYNGO conversion
genes <- union(all_up, all_dn)
genes <- genes[-grep("Rik", genes)]
genes <- genes[-grep("\\(", genes)]
writeClipboard(genes)

