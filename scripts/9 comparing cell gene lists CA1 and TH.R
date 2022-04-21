library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(enrichR)
library(treemap)

#compare number of DE genes for each cell type in CA1 and TH
ca1_up <-  read.delim("gene lists/ca1_genes_up.txt", header = FALSE)$V1
ca1_dn <-  read.delim("gene lists/ca1_genes_dn.txt", header = FALSE)$V1
th_up <-  read.delim("gene lists/th_genes_up.txt", header = FALSE)$V1
th_dn <-  read.delim("gene lists/th_genes_dn.txt", header = FALSE)$V1

nm_dt <- list()
overlap <- list()
unique <- list()
jaccard <- list()

for (i in Sys.glob("gene lists final cell type/*")) {
  x <- gsub("_genes.txt", "", gsub("gene lists final cell type/", "", i))
  print(x)
  
  tmp <- read.delim(i, header = FALSE)$V1
  
  nm_dt[[x]] <- data.frame(dir = c("up", "dn"),
                      CA1 = c(0,0),
                      TH = c(0,0))
  nm_dt[[x]][1, "CA1"] <- intersect(tmp, ca1_up) %>% length
  nm_dt[[x]][2, "CA1"] <- intersect(tmp, ca1_dn) %>% length
  nm_dt[[x]][1, "TH"] <- intersect(tmp, th_up) %>% length
  nm_dt[[x]][2, "TH"] <- intersect(tmp, th_dn) %>% length
  nm_dt[[x]] <- nm_dt[[x]] %>% reshape2::melt(id.vars="dir") %>% mutate(cell_type=x)
  
  unique[[x]] <- data.frame(dir = c("up", "dn"),
                           CA1 = c(0,0),
                           TH = c(0,0))
  unique[[x]][1, "CA1"] <- intersect(tmp, ca1_up) %>% setdiff(th_up) %>% length
  unique[[x]][2, "CA1"] <- intersect(tmp, ca1_dn) %>% setdiff(th_dn) %>% length
  unique[[x]][1, "TH"] <- intersect(tmp, th_up) %>% setdiff(ca1_up) %>% length
  unique[[x]][2, "TH"] <- intersect(tmp, th_dn) %>% setdiff(ca1_dn) %>% length
  unique[[x]] <- unique[[x]] %>% reshape2::melt(id.vars="dir") %>% mutate(cell_type=x)

  overlap[[x]] <- data.frame(cat = c("TH_up", "TH_dn"),
                        CA1_up = c(0,0),
                        CA1_dn = c(0,0))
  overlap[[x]][1,2] <- intersect(intersect(tmp, ca1_up), intersect(tmp, th_up)) %>% length
  overlap[[x]][2,2] <- intersect(intersect(tmp, ca1_up), intersect(tmp, th_dn)) %>% length
  overlap[[x]][2,3] <- intersect(intersect(tmp, ca1_dn), intersect(tmp, th_dn)) %>% length
  overlap[[x]][1,3] <- intersect(intersect(tmp, ca1_dn), intersect(tmp, th_up)) %>% length
  overlap[[x]] <- overlap[[x]] %>% reshape2::melt(id.vars="cat") %>% mutate(cell_type=x)
  
  jaccard[[x]] <- data.frame(dir = c("up", "dn"),
                             index = c(0,0))
  jaccard[[x]][1,2] <- length(intersect(intersect(tmp, ca1_up), intersect(tmp, th_up)))/length(union(intersect(tmp, ca1_up), intersect(tmp, th_up)))
  jaccard[[x]][2,2] <- length(intersect(intersect(tmp, ca1_dn), intersect(tmp, th_dn)))/length(union(intersect(tmp, ca1_dn), intersect(tmp, th_dn)))
  jaccard[[x]] <- jaccard[[x]] %>% reshape2::melt(id.vars="dir") %>% mutate(cell_type=x)

}
nm_dt <- do.call(rbind, nm_dt)
overlap <- do.call(rbind, overlap)
unique <- do.call(rbind, unique)
jaccard <- do.call(rbind, jaccard)

#colors for plots
gn <- brewer.pal(9, "Greens")
or <- brewer.pal(9, "Oranges")
pu <- brewer.pal(9, "Purples")

#make heatmap of number DE genes
lvls <- c("microglia", "astrocyte", "vascular", "oligodendrocyte", "opc", "neuronal")

#increased, all
tmp <- nm_dt %>% filter(dir=="up", !cell_type %in% c("all", "unaccounted"))
tmp$cell_type <- factor(tmp$cell_type, levels = lvls)
ggplot(tmp, aes(x=variable, y=cell_type, fill=value, label=value)) +
  geom_tile() +
  geom_text() +
  ggtitle("Up") +
  scale_fill_gradientn(colours = or) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none")
ggsave("plots/comparing cell type ca1 and th/up_genes.svg", height = 4, width = 2)

#increased, unique
tmp <- unique %>% filter(dir=="up", !cell_type %in% c("all", "unaccounted"))
tmp$cell_type <- factor(tmp$cell_type, levels = lvls)
ggplot(tmp, aes(x=variable, y=cell_type, fill=value, label=value)) +
  geom_tile() +
  geom_text() +
  ggtitle("Up") +
  scale_fill_gradientn(colours = or) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none")
ggsave("plots/comparing cell type ca1 and th/unique_up_genes.svg", height = 4, width = 2)

#decreased, all
tmp <- nm_dt %>% filter(dir=="dn", !cell_type %in% c("all", "unaccounted"))
tmp$cell_type <- factor(tmp$cell_type, levels = lvls)
ggplot(tmp, aes(x=variable, y=cell_type, fill=value, label=value)) +
  geom_tile() +
  geom_text() +
  ggtitle("Down") +
  scale_fill_gradientn(colours = pu) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none")
ggsave("plots/comparing cell type ca1 and th/dn_genes.svg", height = 4, width = 2)

#decreased, unique
tmp <- unique %>% filter(dir=="dn", !cell_type %in% c("all", "unaccounted"))
tmp$cell_type <- factor(tmp$cell_type, levels = lvls)
ggplot(tmp, aes(x=variable, y=cell_type, fill=value, label=value)) +
  geom_tile() +
  geom_text() +
  ggtitle("Down") +
  scale_fill_gradientn(colours = pu) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none")
ggsave("plots/comparing cell type ca1 and th/unique_dn_genes.svg", height = 4, width = 2)

#make heatmap of overlap
tmp <- overlap %>% filter(!cell_type %in% c("all", "unaccounted"))
tmp$cell_type <- factor(tmp$cell_type, levels = lvls)

ggplot(tmp, aes(x=cat, y=variable, fill=value, label=value)) +
  geom_tile() +
  geom_text() +
  ggtitle("Overlap") +
  scale_fill_gradientn(colours = gn) +
  facet_wrap(~cell_type, ncol = 2) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none")
ggsave("plots/comparing cell type ca1 and th/overlap_genes.svg", height = 4.5, width = 2.75)

#make heatmap jaccard indices
tmp <- jaccard %>% filter(!cell_type %in% c("all", "unaccounted"))
tmp$cell_type <- factor(tmp$cell_type, levels = lvls)

ggplot(tmp, aes(x=dir, y=cell_type, fill=value, label=round(value, 3))) +
  geom_tile() +
  geom_text() +
  ggtitle("Jaccard index") +
  scale_fill_gradientn(colours = gn) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position = "none")
ggsave("plots/comparing cell type ca1 and th/jaccard_index.svg", height = 4, width = 2)

###Enrichment results for revigo
#All genes
for (i in Sys.glob("enrichr results/*")) {
  x <- gsub(".txt", "", gsub("enrichr results/", "", i))
  print(x)
  tmp <- read.delim(i) %>% filter(database == "GO_Biological_Process_2021")
  rownames(tmp) <- seq(1, nrow(tmp))
  
  up <- tmp %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", tmp$Term))) %>%
    filter(direction == "up" & list > 1 & P.value < 0.05) %>%
    arrange(P.value) %>% 
    select(Term, P.value)
  write.table(up, paste0("lists for revigo/", x, "_up.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  dn <- tmp %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", tmp$Term))) %>%
    filter(direction == "dn" & list > 1 & P.value < 0.05) %>%
    arrange(P.value) %>% 
    select(Term, P.value)
  write.table(dn, paste0("lists for revigo/", x, "_dn.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

#Unique to CA1, TH
for (i in Sys.glob("gene lists final cell type/*")) {
  x <- gsub("_genes.txt", "", gsub("gene lists final cell type/", "", i))
  print(x)
  
  tmp <-  read.delim(i, header = FALSE)$V1
  
  #common up and dn
  res <- tmp %>% intersect(ca1_up) %>% intersect(th_up)
  write.table(res, paste0("gene lists region specific/common/", x, "_up.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- tmp %>% intersect(ca1_dn) %>% intersect(th_dn)
  write.table(res, paste0("gene lists region specific/common/", x, "_dn.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  #unique ca1 up
  res <- tmp %>% intersect(ca1_up) %>% setdiff(th_up)
  write.table(res, paste0("gene lists region specific/unique/", x, "_ca1_up.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/unique_CA1_up_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
  #unique ca1 dn
  res <- tmp %>% intersect(ca1_dn) %>% setdiff(th_dn)
  write.table(res, paste0("gene lists region specific/unique/", x, "_ca1_dn.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/unique_CA1_dn_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
  #unique th up
  res <- tmp %>% intersect(th_up) %>% setdiff(ca1_up)
  write.table(res, paste0("gene lists region specific/unique/", x, "_th_up.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/unique_TH_up_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
  #unique th dn
  res <- tmp %>% intersect(th_dn) %>% setdiff(ca1_dn)
  write.table(res, paste0("gene lists region specific/unique/", x, "_th_dn.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/unique_TH_dn_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
  #all ca1 up
  res <- tmp %>% intersect(ca1_up)
  write.table(res, paste0("gene lists region specific/all/", x, "_ca1_up.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/CA1_up_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
  #all ca1 dn
  res <- tmp %>% intersect(ca1_dn)
  write.table(res, paste0("gene lists region specific/all/", x, "_ca1_dn.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/CA1_dn_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
  #all th up
  res <- tmp %>% intersect(th_up)
  write.table(res, paste0("gene lists region specific/all/", x, "_th_up.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/TH_up_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
  #all th dn
  res <- tmp %>% intersect(th_dn)
  write.table(res, paste0("gene lists region specific/all/", x, "_th_dn.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  res <- enrichr(res, "GO_Biological_Process_2021")
  res <- do.call(rbind, res)
  res <- res %>%
    mutate(Term = gsub("\\).*", "", gsub(".*\\(", "", res$Term))) %>% 
    filter(P.value < 0.05) %>% 
    select(Term, P.value)
  write.table(res, paste0("lists for revigo/TH_dn_", x, ".txt"), row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
  
}


#make tree maps from revigo results
for (i in Sys.glob("revigo results/*")) {
  x <- gsub("revigo results/RevigoTreeMap_", "", gsub(".csv", "", i))
  print(x)
  
  tmp <- read.csv(i, skip = 4)
  tmp$Representative <- gsub("\"", "", gsub(" \"", "", tmp$Representative))
  tmp$Name <- gsub("\"", "", gsub(" \"", "", tmp$Name))
  tmp$Value <- abs(tmp$Value)
  tmp[tmp$Representative == " null",]$Representative <- tmp[tmp$Representative == " null",]$Name
  
  svg(paste0("plots/comparing cell type ca1 and th/", x, "_treemap.svg"), width = 4, height = 3.5)
  treemap(
    tmp,
    index = c("Representative","Name"),
    vSize = "Value",
    type = "categorical",
    vColor = "Representative",
    title = x,
    palette = "Dark2",
    inflate.labels = TRUE,      # set this to TRUE for space-filling group labels - good for posters
    lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
    bg.labels = 125,   # define background color of group labels
    position.legend = "none"
  )
  dev.off()
}

