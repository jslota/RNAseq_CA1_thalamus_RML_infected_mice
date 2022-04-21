#Heatmaps of genes in CA1 and TH
#2021-09-09
#Jessy Slota

library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)
library(ggplot2)

###Make heatmaps of functional gene categories
###lists of colors
up_cols <- brewer.pal(9 ,"Oranges")[4:9]
dn_cols <- brewer.pal(9 ,"Purples")[4:9]
cat_cols <- brewer.pal(8, "Dark2")

#get normalized data
ca1_dt <- as.matrix(read.csv("normalized data/CA1_all_normalized_counts.csv", row.names = 1))
th_dt <- as.matrix(read.csv("normalized data/TH_all_normalized_counts.csv", row.names = 1))

#calculate z-scores
ca1_dt <- (ca1_dt - rowMeans(ca1_dt))/rowSds(ca1_dt)
th_dt <- (th_dt - rowMeans(th_dt))/rowSds(th_dt)
hmp_dt <- cbind(ca1_dt, th_dt)
brks <- c(seq(min(na.omit(hmp_dt)), -0.1, length=50), seq(0.1, max(na.omit(hmp_dt)), length=50))

#sample info annotation
ann_dt <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1)[colnames(hmp_dt),c(1,2,4)]
plt_cols <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(100))
ann_colors <- list(`Prions` = c(`neg` = "navy", `pos` = "firebrick"),
                   `Experiment` = c(`RML_CA1` = "grey", `RML_TH` = "black"),
                   `DE` = c(`CA1` = cat_cols[1], `TH` = cat_cols[2], `Both` = cat_cols[3]))

#Make list of genes as either DE in CA1, TH or Both
tmp <- read.delim("gene lists final cell type/all_genes.txt", header = FALSE)
ca1 <- union(read.delim("gene lists/ca1_genes_dn.txt", header = FALSE)$V1,
             read.delim("gene lists/ca1_genes_up.txt", header = FALSE)$V1)
th <- union(read.delim("gene lists/th_genes_dn.txt", header = FALSE)$V1,
            read.delim("gene lists/th_genes_up.txt", header = FALSE)$V1)  
both <- intersect(ca1, th)
ca1 <- setdiff(ca1, both)
th <- setdiff(th, both)
tmp <- tmp %>% mutate(DE = case_when(V1 %in% both ~ "Both",
                                     V1 %in% ca1 ~ "CA1",
                                     V1 %in% th ~ "TH"))
gene_list <- data.frame(row.names = tmp$V1,
                        DE = tmp$DE)
rm(tmp, ca1, th, both)

#Make heatmap of every cell type group of genes
lists <- Sys.glob("final lists of genes/*")
lists <- lists[c(2,3,4,5,6,8)]
for (i in lists) {
  x <- gsub("final lists of genes/", "", gsub("_genes.txt", "", i))
  print(x)
  print(paste0("plots/cell type specific/", x, "_hmp.png"))
  #make heatmap
  tmp <- read.delim(i, header = FALSE)$V1
  hmp <- pheatmap(hmp_dt[tmp,], plt_cols, breaks = brks, treeheight_row = 15, cluster_cols = FALSE,
                  annotation_col = ann_dt, annotation_colors = ann_colors, border_color = NA, 
                  show_rownames = FALSE, show_colnames = FALSE, annotation_row = gene_list)
  ggsave(paste0("plots/cell type specific/", x, "_hmp.png"), hmp, width = 4.5, height = 7, dpi = 1200)
}



#Make enrichment plot for top pathways
lists <- Sys.glob("enrichr results/*")
lists <- lists[c(2,3,4,5,6,8)]
for (i in lists) {
  x <- gsub("enrichr results/", "", gsub("_genes.txt", "", i))
  print(x)
  tmp <- read.delim(i)
  rownames(tmp) <- seq(1, nrow(tmp))
  up <- tmp %>% filter(direction == "up" & set < 500) %>% arrange(Adjusted.P.value) %>% select(Term, list, set, Adjusted.P.value, Combined.Score) %>%
    dplyr::slice(n = 1:20) %>% mutate(Term = factor(Term, levels = Term))
  dn <- tmp %>% filter(direction == "dn" & set < 500) %>% arrange(Adjusted.P.value) %>% select(Term, list, set, Adjusted.P.value, Combined.Score) %>%
    dplyr::slice(n = 1:20) %>% mutate(Term = factor(Term, levels = Term))
  
  print(paste0("plots/cell type specific/", x, "_up_enrichplot.svg"))
  ggplot(up, aes(y=Term, x=-log10(`Adjusted.P.value`), color=`Combined.Score`)) +
    geom_point(size=3, alpha=0.66) +
    scale_color_gradientn(colors = up_cols) +
    scale_y_discrete(label = function(x) stringr::str_trunc(x, 40)) +
    theme_classic() +
    theme(legend.position = c(0.8,0.8))
  ggsave(paste0("plots/cell type specific/", x, "_up_enrichplot.svg"), width = 4.5, height = 3.5)
  
  print(paste0("plots/cell type specific/", x, "_dn_enrichplot.svg"))
  ggplot(dn, aes(y=Term, x=-log10(`Adjusted.P.value`), color=`Combined.Score`)) +
    geom_point(size=3, alpha=0.66) +
    scale_color_gradientn(colors = dn_cols) +
    scale_y_discrete(label = function(x) stringr::str_trunc(x, 40)) +
    theme_classic() +
    theme(legend.position = c(0.8,0.8))
  ggsave(paste0("plots/cell type specific/", x, "_dn_enrichplot.svg"), width = 4.5, height = 3.5)
}
