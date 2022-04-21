library(dplyr)

cell_types <- read.delim("cell type databases/final_cell_type_list.txt") %>% 
  rename(human_symbol = hs, mouse_symbol = ms)
gene_ids <- readxl::read_excel("gene lists/idmap.xlsx") %>% 
  rename(human_symbol = query, mouse_symbol = symbol)


#Thalamus
th <- read.csv("results/TH_DE_150_dpi.csv") %>% 
  filter(baseMean>15, abs(log2FoldChange)>0.5, padj<0.05, X %in% cell_types$mouse_symbol) %>% 
  rename(mouse_symbol = X)

th <- th %>%
  left_join(gene_ids) %>% 
  left_join(cell_types) %>% 
  select(mouse_symbol, baseMean, log2FoldChange, lfcSE,
         stat, pvalue, padj, cell_type, human_symbol, name, alias)
write.csv(th, "supplementary_thalamus_150_dpi_prion_altered_transcripts.csv", row.names = FALSE)

#CA1
ca1_neu_mkrs <- union(read.delim("gene lists region specific/ca1_neuron_and_enriched_up.txt", header = FALSE)$V1,
                      read.delim("gene lists region specific/ca1_neuron_and_enriched_dn.txt", header = FALSE)$V1)
cell_types[is.element(cell_types$mouse_symbol, ca1_neu_mkrs),]$cell_type <- "neuron"
ca1 <- read.csv("results/CA1_DE_150_dpi.csv") %>% 
  filter(baseMean>15, abs(log2FoldChange)>0.5, padj<0.05, X %in% cell_types$mouse_symbol) %>% 
  rename(mouse_symbol = X)

ca1 <- ca1 %>%
  left_join(gene_ids) %>% 
  left_join(cell_types) %>% 
  select(mouse_symbol, baseMean, log2FoldChange, lfcSE,
         stat, pvalue, padj, cell_type, human_symbol, name, alias)
write.csv(ca1, "supplementary_CA1_150_dpi_prion_altered_transcripts.csv", row.names = FALSE)
