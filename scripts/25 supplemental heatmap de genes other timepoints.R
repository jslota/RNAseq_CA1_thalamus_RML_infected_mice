#Heatmaps of genes in CA1 and TH
#2021-09-09
#Jessy Slota

library(ComplexHeatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(enrichR)

###Make heatmaps of functional gene categories
#get normalized data
ca1_dt <- as.matrix(read.csv("normalized data/CA1_all_normalized_counts.csv", row.names = 1))
th_dt <- as.matrix(read.csv("normalized data/TH_all_normalized_counts.csv", row.names = 1))

#calculate z-scores
ca1_dt <- (ca1_dt - rowMeans(ca1_dt))/rowSds(ca1_dt)
th_dt <- (th_dt - rowMeans(th_dt))/rowSds(th_dt)

hmp_dt <- cbind(ca1_dt, th_dt)
hmp_dt <- hmp_dt[read.delim("gene lists/de_genes_other_timepoints.txt", header = FALSE)$V1,]

#sample info annotation
ann_dt <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1)[colnames(hmp_dt),c(1,2,4)]

ann_dt <- ann_dt %>% 
  mutate(region=case_when(Experiment=="RML_CA1"~"CA1",
                          Experiment=="RML_TH"~"thalamus")) %>%
  mutate(Treatment=case_when(Prions=="neg"~"Mock",
                             Prions=="pos"~"RML"))

z_scale = circlize::colorRamp2(c(seq(min(hmp_dt), -0.75, length.out=5), 0, seq(0.75, max(hmp_dt), length.out=5)), rev(brewer.pal(11, "PuOr")))

gns <- brewer.pal(9, "Greens")
ha <- HeatmapAnnotation(Treatment=ann_dt$Treatment,
                        dpi=as.factor(ann_dt$dpi),
                        col=list(Treatment=c("Mock"="navy", "RML"="firebrick"),
                                 dpi=c(`70`=gns[2], `90`=gns[4], `130`=gns[7], `150`=gns[9])))

pdf("plots/supplementary plots/early_de_genes_hmp.pdf", width = 6, height = 7)
Heatmap(hmp_dt, col = z_scale, name = "Z-score", show_row_names = FALSE, show_column_names = FALSE, 
        column_split = ann_dt$region, top_annotation = ha)
dev.off()


#microglia associated genes
tmp <- read.delim("gene lists/de_genes_other_timepoints.txt", header = FALSE)$V1
hmp <- pheatmap(hmp_dt[tmp,], plt_cols, breaks = brks, treeheight_row = 0, cluster_cols = FALSE,
                annotation_col = ann_dt, annotation_colors = ann_colors, border_color = NA, 
                show_rownames = FALSE, show_colnames = FALSE)
ggsave("plots/de_genes_other_timepoints_heatmap.png", hmp, width = 5, height = 6, dpi = 1200)

#get pathways for enrichment analysis
dbs <- c("BioPlanet_2019", "GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 
erch <- enrichr(tmp, dbs)
erch <- do.call(rbind, erch)

erch <- erch %>% arrange(Adjusted.P.value) %>% select(Term, Adjusted.P.value, Combined.Score) %>%
  dplyr::slice(n = 1:20) %>% mutate(Term = factor(Term, levels = Term))

plt_cols <- brewer.pal(9 ,"Greens")[4:9]
ggplot(erch, aes(y=Term, x=-log10(`Adjusted.P.value`), color=`Combined.Score`)) +
  geom_point(size=3, alpha=0.66) +
  scale_color_gradientn(colors = plt_cols) +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  theme_classic() +
  theme(legend.position = c(0.8,0.8))
ggsave("plots/de_genes_other_timepoints_enrichplot.svg", width = 5, height = 4)
