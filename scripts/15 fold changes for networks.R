library(dplyr)

fcvals <- read.csv("all_fold_change_values.csv", row.names = 2)[,-1]

#TH astrocyte network
tmp <- read.delim("ppi networks/th_astrocyte_network_genes.txt", header = FALSE)$V1
data.frame(name=tmp, log2FC=fcvals[tmp,]$TH_DE_150_dpi) %>% 
  write.table("ppi networks/th_astrocyte_fold_changes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#CA1 astrocyte network
tmp <- read.delim("ppi networks/ca1_astrocyte_network_genes.txt", header = FALSE)$V1
data.frame(name=tmp, log2FC=fcvals[tmp,]$CA1_DE_150_dpi) %>% 
  write.table("ppi networks/ca1_astrocyte_fold_changes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#TH neuron network
tmp <- read.delim("ppi networks/th_neuron_network_genes.txt", header = FALSE)$V1
data.frame(name=tmp, log2FC=fcvals[tmp,]$TH_DE_150_dpi) %>% 
  write.table("ppi networks/th_neuron_fold_changes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#CA1 neuron network
tmp <- read.delim("ppi networks/ca1_neuron_network_genes.txt", header = FALSE)$V1
data.frame(name=tmp, log2FC=fcvals[tmp,]$CA1_DE_150_dpi) %>% 
  write.table("ppi networks/ca1_neuron_fold_changes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#common astrocyte marker network
tmp <- read.delim("ppi networks/common_astrocyte_markers.txt", header = FALSE)$V1
data.frame(name=tmp, log2FC=fcvals[tmp,]$TH_DE_150_dpi) %>% 
  write.table("ppi networks/common_astrocyte_fold_changes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#common microglia marker network
tmp <- read.delim("ppi networks/common_microglia_markers.txt", header = FALSE)$V1
data.frame(name=tmp, log2FC=fcvals[tmp,]$TH_DE_150_dpi) %>% 
  write.table("ppi networks/common_microglia_fold_changes.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#common vascular marker network
tmp <- read.delim("ppi networks/common_vascular_markers.txt", header = FALSE)$V1
data.frame(name=tmp, log2FC=fcvals[tmp,]$TH_DE_150_dpi) %>% 
  write.table("ppi networks/common_vascular_fold_changes.txt", row.names = FALSE, quote = FALSE, sep = "\t")
