library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(enrichR)

#load data
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
raw_dat <- raw_dat[,rownames(smpl_dat)]
smpl_dat <- smpl_dat %>% mutate(region=gsub(".*_", "", Experiment))

#normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~region+Prions)
dds <- DESeq(dds)
norm_dat <- vst(dds)

pcaData <- plotPCA(norm_dat, intgroup=c("region", "Prions"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=region, shape=Prions)) +
  geom_point(size=2, alpha=0.66) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,2)]) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_classic() +
  theme(text = element_text(size = 9), plot.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("plots/supplementary plots/PCA_ca1_th.svg", width = 5, height = 3.75)

#get genes for heatmap
resultsNames(dds)
res <- results(dds, c("region", "CA1", "TH"))
res <- res %>% as.data.frame() %>% arrange(padj) %>% na.omit()

genes <- res %>% filter(baseMean>25, abs(log2FoldChange) > 1, padj < 1e-08) %>% rownames()
genes <- genes[-grep("RP", genes)]
genes <- genes[-grep("AC", genes)]
genes <- genes[-grep("Gm", genes)]
genes <- genes[-grep("Rik", genes)]

hmp_dt <- assay(norm_dat)[genes,]
hmp_dt <- (hmp_dt-rowMeans(hmp_dt))/matrixStats::rowSds(hmp_dt)

zscore <- circlize::colorRamp2(c(-4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4), rev(brewer.pal(11, "PuOr")))

ha <- HeatmapAnnotation(region=smpl_dat$region, col = list(region=c("CA1"=brewer.pal(3, "Dark2")[1],
                                                                    "TH"=brewer.pal(3, "Dark2")[2])))
#draw heatmap
pdf("plots/supplementary plots/ca1_th_hmp.pdf", width = 6, height = 8)
ht <- Heatmap(hmp_dt, col = zscore, name="Z-score", show_row_names = FALSE, show_column_names = FALSE, top_annotation = ha)
draw(ht, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

#enrichment analysis
dbs <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "PanglaoDB_Augmented_2021") 
up_cols <- brewer.pal(9 ,"Oranges")[4:9]
dn_cols <- brewer.pal(9 ,"Purples")[4:9]

#Increased in CA1
tmp <- res %>% filter(baseMean>25, log2FoldChange > 1, padj < 1e-08) %>% rownames() %>% intersect(genes)
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, P.value, Combined.Score) %>%
  dplyr::slice(n = 1:20) %>% mutate(Term = factor(Term, levels = Term))

ggplot(plt_dt, aes(y=Term, x=-log10(`P.value`), color=`Combined.Score`)) +
  geom_point(size=3, alpha=0.66) +
  scale_color_gradientn(colors = up_cols) +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  theme_classic() +
  theme(legend.position = c(0.8,0.8))
ggsave("plots/supplementary plots/up_ca1_enrchplt.svg", width = 5, height = 4)

#Increased in Thalamus
tmp <- res %>% filter(baseMean>25, log2FoldChange < -1, padj < 1e-08) %>% rownames() %>% intersect(genes)
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, P.value, Combined.Score) %>%
  dplyr::slice(n = 1:20) %>% mutate(Term = factor(Term, levels = Term))

ggplot(plt_dt, aes(y=Term, x=-log10(`P.value`), color=`Combined.Score`)) +
  geom_point(size=3, alpha=0.66) +
  scale_color_gradientn(colors = dn_cols) +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  theme_classic() +
  theme(legend.position = c(0.8,0.8))
ggsave("plots/supplementary plots/dn_ca1_enrchplt.svg", width = 5, height = 4)
