library(dplyr)
library(ggplot2)
library(enrichR)
library(ComplexHeatmap)

markers <- read.delim("gene lists region specific/common/microglia_up.txt", header = FALSE)$V1

#Heatmap of fold changes
fc_dt <- read.csv("all_fold_change_values.csv", row.names = 1) %>% 
  select(X, TH_DE_70_dpi, TH_DE_90_dpi, TH_DE_130_dpi, TH_DE_150_dpi,
         CA1_DE_70_dpi, CA1_DE_90_dpi, CA1_DE_130_dpi, CA1_DE_150_dpi) %>% 
  filter(X %in% markers)

hmp_dt <- fc_dt
rownames(hmp_dt) <- hmp_dt$X
hmp_dt <- hmp_dt[,-1]
hmp_dt <- as.matrix(hmp_dt)

fc_dt %>% arrange(desc(TH_DE_150_dpi))

gns <- brewer.pal(9, "Greens")
ann_dt <- data.frame(row.names = colnames(hmp_dt),
                     dpi = rep(c(70, 90, 130, 150), 2),
                     region=c(rep("thalamus", 4), rep("CA1", 4)))

ha <- HeatmapAnnotation(dpi=as.factor(ann_dt$dpi),
                        col=list(dpi=c(`70`=gns[2], `90`=gns[4], `130`=gns[7], `150`=gns[9])))

FC <- circlize::colorRamp2(c(min(hmp_dt), -3, -2.25, -1.5, -0.75, 0, 0.75, 1.5, 2.25, 3, max(hmp_dt)), rev(brewer.pal(11, "PuOr")))

mkrs <- read.delim("cell type databases/micro markers hmap.txt", header=FALSE)$V1
mkrs_plot <- rownames(hmp_dt) %in% mkrs

pdf("plots/microglia analysis/common_micro_hmp.pdf", width = 5, height = 6)
Heatmap(hmp_dt, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
        show_column_names = FALSE, column_split = ann_dt$region) +
  rowAnnotation(foo = anno_mark(at=which(mkrs_plot), labels = rownames(hmp_dt)[mkrs_plot], labels_gp = gpar(fontsize=8)))
dev.off()

###draw legends for ppi networks

#microglia PPI network
col_fun = circlize::colorRamp2(seq(0, 10, length.out=9), brewer.pal(9, "Oranges"))
lgd = Legend(col_fun = col_fun, title = "Log2(FC)")
pdf("plots/legend_microglia_ppi.pdf", width=1, height=2)
draw(lgd)
dev.off()

#other PPI netowkrs
col_fun = circlize::colorRamp2(seq(-3, 3, length.out=11), rev(brewer.pal(11, "PuOr")))
lgd = Legend(col_fun = col_fun, title = "Log2(FC)")
pdf("plots/legend_ppi.pdf", width=1, height=2)
draw(lgd)
dev.off()

#Enrichr analysis
up_cols <- brewer.pal(9 ,"Oranges")[4:9]
dn_cols <- brewer.pal(9 ,"Purples")[4:9]
dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 

#CA1 up
tmp <- read.delim("gene lists region specific/all/microglia_ca1_up.txt", header = FALSE)$V1
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
ggsave("plots/microglia analysis/ca1_up_enrichplt.svg", width = 4, height = 2.5)

#CA1 dn
tmp <- read.delim("gene lists region specific/all/microglia_ca1_dn.txt", header = FALSE)$V1
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
ggsave("plots/microglia analysis/ca1_dn_enrichplt.svg", width = 4, height = 2.5)

#th up
tmp <- read.delim("gene lists region specific/all/microglia_th_up.txt", header = FALSE)$V1
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
ggsave("plots/microglia analysis/th_up_enrichplt.svg", width = 4, height = 2.5)

#th dn
tmp <- read.delim("gene lists region specific/all/microglia_th_dn.txt", header = FALSE)$V1
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
ggsave("plots/microglia analysis/th_dn_enrichplt.svg", width = 4, height = 2.5)
