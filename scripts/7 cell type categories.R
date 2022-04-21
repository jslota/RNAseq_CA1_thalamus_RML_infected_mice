
library(dplyr)

gene_symbol <- na.omit(readxl::read_excel("gene lists/idmap.xlsx")[,c(1,4)])
colnames(gene_symbol) <- c("hs", "ms")

tmp <- union(read.delim("gene lists/all_up.txt", header = FALSE)$V1,
             read.delim("gene lists/all_dn.txt", header = FALSE)$V1)
gene_symbol <- gene_symbol %>% filter(ms %in% tmp)

dbs1 <- c("top_all", "top_mouse", "top_human")
dbs2 <- c("_enrich", "_expression", "_specificity")

gene_symbol <- gene_symbol %>% mutate(cell_type = NA, category = NA)

#add in manually annotated genes
tmp <- read.delim("cell type databases/manually categorized genes.txt")
tmp %>% pull(cell_type) %>% unique()
gene_symbol <- gene_symbol %>% mutate(cell_type = case_when(hs %in% tmp[tmp$cell_type == "mic",]$gene ~ "microglia",
                                                            hs %in% tmp[tmp$cell_type == "ast",]$gene ~ "astrocyte",
                                                            hs %in% tmp[tmp$cell_type == "end",]$gene ~ "endothelial",
                                                            hs %in% tmp[tmp$cell_type == "oli",]$gene ~ "oligodendrocyte",
                                                            hs %in% tmp[tmp$cell_type == "opc",]$gene ~ "opc",
                                                            hs %in% tmp[tmp$cell_type == "neu",]$gene ~ "neuron"),
                                      category = "manual")


#Add in annotated genes from database
tmp <- gene_symbol %>% filter(is.na(cell_type) == TRUE)
for (i in dbs1) {
  for (k in dbs2) {
    cell_db <- readxl::read_excel("cell type databases/cell_type_database.xlsx", 
                                  sheet = paste0(i, k), skip = 2)
    tmp <- tmp %>% mutate(cell_type = case_when(hs %in% cell_db[cell_db$Celltype == "mic",]$gene ~ "microglia",
                                                hs %in% cell_db[cell_db$Celltype == "ast",]$gene ~ "astrocyte",
                                                hs %in% cell_db[cell_db$Celltype == "end",]$gene ~ "endothelial",
                                                hs %in% cell_db[cell_db$Celltype == "oli",]$gene ~ "oligodendrocyte",
                                                hs %in% cell_db[cell_db$Celltype == "opc",]$gene ~ "opc",
                                                hs %in% cell_db[cell_db$Celltype == "neu",]$gene ~ "neuron"),
                          category = paste0(i, k))
    gene_symbol[is.element(gene_symbol$hs, na.omit(tmp)$hs),] <- na.omit(tmp)
    tmp <- tmp %>% filter(is.na(cell_type) == TRUE)
  }
}

uncategorized <- read.delim("cell type databases/manually categorized genes.txt")
setdiff(tmp$hs, uncategorized$gene) %>% writeClipboard()
rm(tmp, uncategorized)

gene_symbol <- gene_symbol %>% arrange(hs)

gene_symbol %>% filter(is.na(cell_type) == TRUE) %>% nrow()

table(gene_symbol$cell_type)

write.table(gene_symbol %>% filter(cell_type == "neuron") %>% pull(ms), "gene lists final cell type/neuronal_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol %>% filter(cell_type == "astrocyte") %>% pull(ms), "gene lists final cell type/astrocyte_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol %>% filter(cell_type == "microglia") %>% pull(ms), "gene lists final cell type/microglia_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol %>% filter(cell_type == "endothelial") %>% pull(ms), "gene lists final cell type/vascular_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol %>% filter(cell_type == "oligodendrocyte") %>% pull(ms), "gene lists final cell type/oligodendrocyte_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol %>% filter(cell_type == "opc") %>% pull(ms), "gene lists final cell type/opc_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol %>% filter(is.na(cell_type) == TRUE) %>% pull(ms), "gene lists final cell type/unaccounted_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol %>% pull(ms), "gene lists final cell type/all_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gene_symbol, "cell type databases/final_cell_type_list.txt", sep = "\t")

