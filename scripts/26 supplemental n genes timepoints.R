library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

lvls <- c("CA1_DE_70_dpi", "CA1_DE_90_dpi", "CA1_DE_130_dpi", "CA1_DE_150_dpi",
          "TH_DE_70_dpi", "TH_DE_90_dpi", "TH_DE_130_dpi", "TH_DE_150_dpi")

smpl_qc_genes <- read.delim("gene lists/poor_sample_de_genes.txt", header = FALSE)$V1

#number of de genes
dt <- list()
for (i in Sys.glob("results/*")) {
  x <- gsub(".csv", "", gsub("results/*", "", i))
  print(x)
  dt[[x]]$up <- read.csv(i) %>% filter(baseMean > 15 & log2FoldChange > 0.5 & padj < 0.05) %>%
    pull(X) %>% setdiff(smpl_qc_genes) %>% length()
  dt[[x]]$dn <- read.csv(i) %>% filter(baseMean > 15 & log2FoldChange < -0.5 & padj < 0.05) %>%
    pull(X) %>% setdiff(smpl_qc_genes) %>% length()
  dt[[x]] <- as.data.frame(dt[[x]])
}
dt <- do.call(rbind, dt)
dt <- dt %>% mutate(comparison = rownames(dt)) %>% melt(id.var = "comparison")
dt <- dt %>% filter(comparison %in% lvls)
tmp <- read.csv("comparisons.csv")
dt <- left_join(dt, tmp, by="comparison")

ggplot(dt, aes(x=dpi, y=value, color=variable, shape=variable, label=value)) +
  geom_line() +
  geom_point(size=2, alpha=0.66) +
  geom_text(size=3, nudge_y = 30) +
  scale_color_manual(values = brewer.pal(11, "PuOr")[c(2,10)]) +
  scale_x_continuous(breaks = c(70, 90, 130, 150)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1500)) +
  facet_wrap(~region) +
  theme_classic()
ggsave("plots/supplementary plots/number_de_genes.svg", width = 4, height = 2.5)

#overlap of de genes
comparisons <- paste0("results/", lvls, ".csv")
dt_up <- matrix(nrow = length(comparisons), ncol = length(comparisons), dimnames = list(comparisons, comparisons))
dt_dn <- matrix(nrow = length(comparisons), ncol = length(comparisons), dimnames = list(comparisons, comparisons))

for (i in paste0("results/", lvls, ".csv")) {
    for (k in comparisons) {
    print(paste0(i, " vs ", k))
    a <- read.csv(i) %>% filter(baseMean > 15 & log2FoldChange > 0.5 & padj < 0.05) %>% pull(X) %>% setdiff(smpl_qc_genes)
    b <- read.csv(k) %>% filter(baseMean > 15 & log2FoldChange > 0.5 & padj < 0.05) %>% pull(X) %>% setdiff(smpl_qc_genes)
    dt_up[i,k] <- length(intersect(a, b))/length(union(a, b))
    
    a <- read.csv(i) %>% filter(baseMean > 15 & log2FoldChange < -0.5 & padj < 0.05) %>% pull(X) %>% setdiff(smpl_qc_genes)
    b <- read.csv(k) %>% filter(baseMean > 15 & log2FoldChange < -0.5 & padj < 0.05) %>% pull(X) %>% setdiff(smpl_qc_genes)
    dt_dn[i,k] <- length(intersect(a, b))/length(union(a, b))
    }
  comparisons <- comparisons[-grep(i, comparisons)]
}
colnames(dt_up) <-  gsub(".csv", "", gsub("results/*", "", colnames(dt_up)))
rownames(dt_up) <-  gsub(".csv", "", gsub("results/*", "", rownames(dt_up)))
colnames(dt_dn) <-  gsub(".csv", "", gsub("results/*", "", colnames(dt_dn)))
rownames(dt_dn) <-  gsub(".csv", "", gsub("results/*", "", rownames(dt_dn)))
colnames(dt_up) <-  gsub("_DE", "", gsub("all_ctrls", "dpi", colnames(dt_up)))
rownames(dt_up) <-  gsub("_DE", "", gsub("all_ctrls", "dpi", rownames(dt_up)))
colnames(dt_dn) <-  gsub("_DE", "", gsub("all_ctrls", "dpi", colnames(dt_dn)))
rownames(dt_dn) <-  gsub("_DE", "", gsub("all_ctrls", "dpi", rownames(dt_dn)))

#hierarchical clustering snippet
dt_up[,hclust(dist(dt_up))$order] %>% colnames()


#Heatmap up gene overlap
tmp <- dt_up %>% 
  melt() %>%
  na.omit()

plt_cls <- colorRampPalette(brewer.pal(9, "Oranges"))(100)
  
ggplot(tmp, aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradientn(colours = plt_cls) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0), 
        axis.title = element_blank())
ggsave("plots/supplementary plots/up_genes_overlap.svg", width = 3.5, height = 2.75)

#Heatmap down gene overlap
tmp <- dt_dn %>% 
  melt() %>%
  na.omit()

plt_cls <- colorRampPalette(brewer.pal(9, "Purples"))(100)

ggplot(tmp, aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradientn(colours = plt_cls) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0), 
        axis.title = element_blank())
ggsave("plots/supplementary plots/dn_genes_overlap.svg", width = 3.5, height = 2.75)
