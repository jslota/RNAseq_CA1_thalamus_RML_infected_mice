
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

#get lists of genes
ca1_enriched <- read.delim("gene lists/ca1_enriched_compared_to_th_genes.txt", header = FALSE)$V1
th_enriched <- read.delim("gene lists/th_enriched_compared_to_ca1_genes.txt", header = FALSE)$V1
neuron <- read.delim("gene lists final cell type/neuronal_genes.txt", header = FALSE)$V1
ca1_altered <- union(read.delim("gene lists/ca1_genes_up.txt", header = FALSE)$V1,
                     read.delim("gene lists/ca1_genes_dn.txt", header = FALSE)$V1)
th_altered <- union(read.delim("gene lists/th_genes_up.txt", header = FALSE)$V1,
                    read.delim("gene lists/th_genes_dn.txt", header = FALSE)$V1)

#CA1 reassigned to neuron
ca1_neuron_markers <- ca1_altered %>% intersect(ca1_enriched) %>% setdiff(neuron)

#TH neuron and TH enriched
th_neuron_markers <- th_altered %>% intersect(neuron) %>% intersect(th_enriched)

#get normalized data
ca1_dt <- as.matrix(read.csv("normalized data/CA1_all_normalized_counts.csv", row.names = 1))
th_dt <- as.matrix(read.csv("normalized data/TH_all_normalized_counts.csv", row.names = 1))

full_dt <- cbind(ca1_dt, th_dt)[c(ca1_neuron_markers, th_neuron_markers),]

hmp_dt <- (full_dt-rowMeans(full_dt))/matrixStats::rowSds(full_dt)

Heatmap(hmp_dt)

#plotcounts
smpls <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1)[colnames(full_dt),c(1,2,4,5)]
summary(rownames(smpls)==colnames(full_dt))

plt_dt <- smpls %>% mutate(region=case_when(Experiment=="RML_CA1"~"CA1",
                                            Experiment=="RML_TH"~"thalamus"),
                           `L1cam`=full_dt["L1cam",],
                           `Cit`=full_dt["Cit",],
                           `Pcp4`=full_dt["Pcp4",],
                           `Rara`=full_dt["Rara",],
                           `Nr4a3`=full_dt["Nr4a3",],
                           `Lzts1`=full_dt["Lzts1",])

#L1cam
ggplot(plt_dt, aes(x=region, y=L1cam, color=region)) +
  geom_boxplot(width=0.15) +
  geom_jitter(width=0.15) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,2)]) +
  ggtitle("L1cam") +
  ylab("log2(read count)") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("plots/supplementary plots/L1cam_plt.svg", width=3.5, height = 3)

#Cit
ggplot(plt_dt, aes(x=region, y=Cit, color=region)) +
  geom_boxplot(width=0.15) +
  geom_jitter(width=0.15) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,2)]) +
  ggtitle("Cit") +
  ylab("log2(read count)") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("plots/supplementary plots/Cit_plt.svg", width=3.5, height = 3)

#Pcp4
ggplot(plt_dt, aes(x=region, y=Pcp4, color=region)) +
  geom_boxplot(width=0.15) +
  geom_jitter(width=0.15) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,2)]) +
  ggtitle("Pcp4") +
  ylab("log2(read count)") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("plots/supplementary plots/Pcp4_plt.svg", width=3.5, height = 3)

#Rara
ggplot(plt_dt, aes(x=region, y=Rara, color=region)) +
  geom_boxplot(width=0.15) +
  geom_jitter(width=0.15) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,2)]) +
  ggtitle("Rara") +
  ylab("log2(read count)") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("plots/supplementary plots/Rara_plt.svg", width=3.5, height = 3)

#Nr4a3
ggplot(plt_dt, aes(x=region, y=Nr4a3, color=region)) +
  geom_boxplot(width=0.15) +
  geom_jitter(width=0.15) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,2)]) +
  ggtitle("Nr4a3") +
  ylab("log2(read count)") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("plots/supplementary plots/Nr4a3_plt.svg", width=3.5, height = 3)

#Lzts1
ggplot(plt_dt, aes(x=region, y=Lzts1, color=region)) +
  geom_boxplot(width=0.15) +
  geom_jitter(width=0.15) +
  scale_color_manual(values = brewer.pal(3, "Dark2")[c(1,2)]) +
  ggtitle("Lzts1") +
  ylab("log2(read count)") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave("plots/supplementary plots/Lzts1_plt.svg", width=3.5, height = 3)

ca1_neuron_markers
th_neuron_markers

tmp <- read.csv("gene lists/DE_res_CA1_vs_TH.csv", row.names = 1)[th_neuron_markers,]

tmp %>% arrange(log2FoldChange)
