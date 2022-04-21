#Heatmaps of genes in CA1 and TH
#2021-09-09
#Jessy Slota

###Complex heatmap
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(enrichR)
library(ggplot2)

###Make heatmap
#get normalized data
ca1_dt <- as.matrix(read.csv("normalized data/CA1_normalized_counts.csv", row.names = 1))
th_dt <- as.matrix(read.csv("normalized data/TH_normalized_counts.csv", row.names = 1))

#calculate z-scores
ca1_dt <- (ca1_dt - rowMeans(ca1_dt))/matrixStats::rowSds(ca1_dt)
th_dt <- (th_dt - rowMeans(th_dt))/matrixStats::rowSds(th_dt)

#get gene info
genes <- read.delim("cell type databases/final_cell_type_list.txt")
rownames(genes) <- genes$ms
genes %>% pull(cell_type) %>% unique
genes <- genes %>% mutate(cell_short = factor(substr(cell_type, 1, 3), levels = c("mic", "ast", "end", "oli", "opc", "neu", NA)))

hmp_dt <- cbind(ca1_dt, th_dt)[genes$ms,]


#get sample info
smpls <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1)[colnames(hmp_dt),c(1,2,4)] %>% 
  mutate(region=case_when(Experiment=="RML_CA1"~"CA1",
                          Experiment=="RML_TH"~"thalamus")) %>%
  mutate(Treatment=case_when(Prions=="neg"~"Mock",
                             Prions=="pos"~"RML"))

z_scale = circlize::colorRamp2(c(seq(min(hmp_dt), -0.75, length.out=5), 0, seq(0.75, -min(hmp_dt), length.out=5)), rev(brewer.pal(11, "PuOr")))

gns <- brewer.pal(9, "Greens")
#add annotation
ha <- HeatmapAnnotation(Treatment=smpls$Treatment,
                        dpi=as.factor(smpls$dpi),
                        col=list(Treatment=c("Mock"="navy", "RML"="firebrick"),
                                 dpi=c(`70`=gns[2], `90`=gns[4], `130`=gns[7], `150`=gns[9])))
dk2 <- brewer.pal(8, "Dark2")

mkrs <- read.delim("cell type databases/gene markers for heatmap.txt", header=FALSE)$V1
mkrs_plot <- rownames(hmp_dt) %in% mkrs

pdf("plots/prion_altered_heatmap.pdf", width=6, height = 10)
ht_list = Heatmap(hmp_dt, name = "Z-score", col=z_scale, row_split = genes$cell_short, column_split = smpls$region, show_row_names = FALSE,
                  top_annotation = ha, show_column_names = FALSE, cluster_row_slices = FALSE, border = FALSE) +
  rowAnnotation(cell_type=genes$cell_type, show_annotation_name=FALSE,
                foo = anno_mark(at=which(mkrs_plot), labels = rownames(hmp_dt)[mkrs_plot], labels_gp = gpar(fontsize=8)),
                col=list(cell_type=c("astrocyte"=dk2[1], "neuron"=dk2[2], "opc"=dk2[3], "microglia"=dk2[4], "endothelial"=dk2[5], "oligodendrocyte"=dk2[6])),
                na_col = "white")
draw(ht_list, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

###Enrichr all increased/decreased genes
dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021")
up_cols <- brewer.pal(9 ,"Oranges")[4:9]
dn_cols <- brewer.pal(9 ,"Purples")[4:9]

tmp <- read.delim("gene lists/all_up.txt", header = FALSE)$V1
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, Overlap, P.value, Combined.Score) %>%
  mutate(Term = stringr::str_trunc(gsub(" \\(..*", "", Term), 40),
         n = as.numeric(gsub("/..*", "", Overlap))) %>%
  distinct(Term, .keep_all = TRUE) %>%
  dplyr::slice(n = 1:10)
plt_dt <- plt_dt %>% mutate(Term = factor(Term, levels = Term))

ggplot(plt_dt, aes(y=Term, x=-log10(`P.value`), color=`Combined.Score`, size=n)) +
  geom_point(alpha=0.66) +
  scale_color_gradientn(colors = up_cols) +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  theme_classic() +
  theme(legend.position = c(0.8,0.8))
ggsave("plots/all_genes_up_enrichr.svg", width = 4, height = 2.5)

tmp <- read.delim("gene lists/all_dn.txt", header = FALSE)$V1
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, Overlap, P.value, Combined.Score) %>%
  mutate(Term = stringr::str_trunc(gsub(" \\(..*", "", Term), 40),
         n = as.numeric(gsub("/..*", "", Overlap))) %>%
  distinct(Term, .keep_all = TRUE) %>%
  dplyr::slice(n = 1:10)
plt_dt <- plt_dt %>% mutate(Term = factor(Term, levels = Term))

ggplot(plt_dt, aes(y=Term, x=-log10(`P.value`), color=`Combined.Score`, size=n)) +
  geom_point(alpha=0.66) +
  scale_color_gradientn(colors = dn_cols) +
  scale_y_discrete(label = function(x) stringr::str_trunc(x, 40)) +
  theme_classic() +
  theme(legend.position = c(0.8,0.8))
ggsave("plots/all_genes_dn_enrichr.svg", width = 4, height = 2.5)