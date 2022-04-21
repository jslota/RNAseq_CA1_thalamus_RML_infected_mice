###Enrichr analysis for each functional category
#2021-09-23
#Jessy Slota

library(enrichR)
library(dplyr)

#get list of all up and down genes
all_up_genes <- read.delim("gene lists/all_up.txt", header = FALSE)$V1
all_dn_genes <- read.delim("gene lists/all_dn.txt", header = FALSE)$V1

#get Human gene symbols for enrichment analysis
gene_symbol <- na.omit(readxl::read_excel("gene lists/idmap.xlsx")[,c(1,4)])
colnames(gene_symbol) <- c("hs", "ms")

#get pathways for enrichment analysis
dbs <- c("BioPlanet_2019", "WikiPathway_2021_Human", "GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021") 

#perform enrichment analysis for each functional category
for (j in Sys.glob("gene lists final cell type/*")) {
  print(gsub(".txt", "", gsub(".*/", "", j)))
  
  tmp <- read.delim(j, header = FALSE)$V1
  
  if (length(intersect(tmp, all_dn_genes)) == 0) {
    up <- intersect(tmp, all_up_genes)
    up <- gene_symbol[is.element(gene_symbol$ms, up),]$ms
    up <- enrichr(up, dbs)
    
    for (i in dbs) {
      up[[i]] <- up[[i]] %>% mutate(database = i, direction = "up") %>%
        tidyr::separate(Overlap, c("list", "set"), sep = "/")
    }
    
    up <- do.call(rbind, up)
    
    write.table(up, paste0("enrichr results/", gsub(".txt", "", gsub(".*/", "", j)), ".txt"), quote = FALSE, sep = "\t")
    
  } else  {
    #up enrichment analysis
    up <- intersect(tmp, all_up_genes)
    up <- gene_symbol[is.element(gene_symbol$ms, up),]$ms
    up <- enrichr(up, dbs)
    
    #dn enrichment analysis
    dn <- intersect(tmp, all_dn_genes)
    dn <- gene_symbol[is.element(gene_symbol$ms, dn),]$ms
    dn <- enrichr(dn, dbs)
    
    for (i in dbs) {
      up[[i]] <- up[[i]] %>% mutate(database = i, direction = "up") %>%
        tidyr::separate(Overlap, c("list", "set"), sep = "/")
      dn[[i]] <- dn[[i]] %>% mutate(database = i, direction = "dn") %>%
        tidyr::separate(Overlap, c("list", "set"), sep = "/")
    }
    
    up <- do.call(rbind, up)
    dn <- do.call(rbind, dn)
    
    res <- rbind(up, dn)
    
    write.table(res, paste0("enrichr results/", gsub(".txt", "", gsub(".*/", "", j)), ".txt"), quote = FALSE, sep = "\t")
  }
  
 }


