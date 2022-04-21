library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(enrichR)


neuron_genes <- read.delim("gene lists final cell type/neuronal_genes.txt", header = FALSE)$V1
cols <- brewer.pal(8, "Dark2")
dbs <- c("BioPlanet_2019", "GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 
up_cols <- brewer.pal(9 ,"Oranges")[4:9]
dn_cols <- brewer.pal(9 ,"Purples")[4:9]

#CA1, dn
tmp <- neuron_genes %>% intersect(read.delim("gene lists/ca1_genes_dn.txt", header = FALSE)$V1)
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
ggsave("plots/neuron analysis/ca1_dn_enrichplt.svg", width = 4, height = 4)

#CA1, up
tmp <- neuron_genes %>% intersect(read.delim("gene lists/ca1_genes_up.txt", header = FALSE)$V1)
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
ggsave("plots/neuron analysis/ca1_up_enrichplt.svg", width = 4, height = 4)

#TH, dn
dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 
tmp <- neuron_genes %>% intersect(read.delim("gene lists/th_genes_dn.txt", header = FALSE)$V1)
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
ggsave("plots/neuron analysis/th_dn_enrichplt.svg", width = 4, height = 2.5)

#TFs
dbs <- c("ChEA_2016") 
tmp <- neuron_genes %>% intersect(read.delim("gene lists/th_genes_dn.txt", header = FALSE)$V1)
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, Overlap, P.value, Combined.Score) %>%
  mutate(Term = gsub(" Ch..*", "", Term),
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
ggsave("plots/neuron analysis/th_TF_dn_enrichplt.svg", width = 3, height = 2.5)

#TH up
dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 
tmp <- neuron_genes %>% intersect(read.delim("gene lists/th_genes_up.txt", header = FALSE)$V1)
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
ggsave("plots/neuron analysis/th_up_enrichplt.svg", width = 4, height = 2.5)

#TFs
dbs <- c("ChEA_2016") 
tmp <- neuron_genes %>% intersect(read.delim("gene lists/th_genes_up.txt", header = FALSE)$V1)
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, Overlap, P.value, Combined.Score) %>%
  mutate(Term = gsub(" Ch..*", "", Term),
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
ggsave("plots/neuron analysis/th_TF_up_enrichplt.svg", width = 3, height = 2.5)

###Include ca1 enriched genes
ca1_enriched_genes <- read.delim("gene lists/ca1_enriched_compared_to_th_genes.txt", header = FALSE)$V1
neuron_genes <- union(ca1_enriched_genes, read.delim("gene lists final cell type/neuronal_genes.txt", header = FALSE)$V1)

#CA1, dn
dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 
tmp <- neuron_genes %>% intersect(read.delim("gene lists/ca1_genes_dn.txt", header = FALSE)$V1)
write.table(tmp, "gene lists region specific/ca1_neuron_and_enriched_dn.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
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
ggsave("plots/neuron analysis/ca1_enriched_dn_enrichplt.svg", width = 4, height = 2.5)

#TFs
dbs <- c("ChEA_2016") 
tmp <- neuron_genes %>% intersect(read.delim("gene lists/ca1_genes_dn.txt", header = FALSE)$V1)
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, Overlap, P.value, Combined.Score) %>%
  mutate(Term = gsub(" Ch..*", "", Term),
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
ggsave("plots/neuron analysis/ca1_enriched_TF_dn_enrichplt.svg", width = 3, height = 2.5)

#CA1, up
dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 
tmp <- neuron_genes %>% intersect(read.delim("gene lists/ca1_genes_up.txt", header = FALSE)$V1)
write.table(tmp, "gene lists region specific/ca1_neuron_and_enriched_up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
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
ggsave("plots/neuron analysis/ca1_enriched_up_enrichplt.svg", width = 4, height = 2.5)

#CA1, up
dbs <- c("ChEA_2016")
tmp <- neuron_genes %>% intersect(read.delim("gene lists/ca1_genes_up.txt", header = FALSE)$V1)
plt_dt <- enrichr(tmp, dbs)
plt_dt <- do.call(rbind, plt_dt)
plt_dt <- plt_dt %>% arrange(P.value) %>% select(Term, Overlap, P.value, Combined.Score) %>%
  mutate(Term = gsub(" Ch..*", "", Term),
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
ggsave("plots/neuron analysis/ca1_enriched_TF_up_enrichplt.svg", width = 3, height = 2.5)
