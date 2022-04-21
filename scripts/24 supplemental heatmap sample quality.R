#heatmap of "poor sample quality" genes across all samples
#2021-10-01
#Jessy Slota

library(ComplexHeatmap)
library(RColorBrewer)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(enrichR)

#get normalized data
ca1_dt <- as.matrix(read.csv("normalized data/CA1_all_normalized_counts.csv", row.names = 1))
th_dt <- as.matrix(read.csv("normalized data/TH_all_normalized_counts.csv", row.names = 1))

#calculate z-scores
ca1_dt <- (ca1_dt - rowMeans(ca1_dt))/rowSds(ca1_dt)
th_dt <- (th_dt - rowMeans(th_dt))/rowSds(th_dt)

hmp_dt <- cbind(ca1_dt, th_dt)
genes <- read.delim("gene lists/poor_sample_de_genes.txt", header = FALSE)$V1
hmp_dt <- hmp_dt[genes,]
brks <- c(seq(min(na.omit(hmp_dt)), -0.1, length=50), seq(0.1, max(na.omit(hmp_dt)), length=50))

###Make heatmaps of functional gene categories
ann_dt <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1)[colnames(hmp_dt),c(1,4)]
ms_no <- ann_dt %>% rownames() %>% gsub("_CA1.*", "", .) %>% gsub("_TH.*", "", .)
ann_dt$ms_no <- ms_no
ann_dt <- ann_dt %>% 
  mutate(region=case_when(Experiment=="RML_CA1"~"CA1",
                          Experiment=="RML_TH"~"thalamus")) %>%
  mutate(mouse_number = gsub(".*_", "", ms_no))

z_scale = circlize::colorRamp2(c(seq(min(hmp_dt), -0.75, length.out=5), 0, seq(0.75, max(hmp_dt), length.out=5)), rev(brewer.pal(11, "PuOr")))

ha <- HeatmapAnnotation(region=ann_dt$region,
                        mouse_number = ann_dt$mouse_number,
                        col=list(region=c("CA1"="lightgrey", "thalamus"="black")))

pdf("plots/supplementary plots/sample_quality_hmp.pdf", width = 6, height = 7)
Heatmap(hmp_dt, col = z_scale, name = "Z-score", show_row_names = FALSE, show_column_names = FALSE, top_annotation = ha)
dev.off()

#get pathways for enrichment analysis
dbs <- c("BioPlanet_2019", "GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 
erch <- enrichr(genes, dbs)
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
ggsave("plots/smpl_qual_genes_enrichplot.svg", width = 5, height = 4)
