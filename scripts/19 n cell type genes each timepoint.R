library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)

#Get overlapping n genes between CA1 and TH
ca1_up <-  read.delim("gene lists/ca1_genes_up.txt", header = FALSE)$V1
ca1_dn <-  read.delim("gene lists/ca1_genes_dn.txt", header = FALSE)$V1
th_up <-  read.delim("gene lists/th_genes_up.txt", header = FALSE)$V1
th_dn <-  read.delim("gene lists/th_genes_dn.txt", header = FALSE)$V1

overlap <- list()
for (i in Sys.glob("gene lists final cell type/*")) {
  x <- gsub("_genes.txt", "", gsub("gene lists final cell type/", "", i))
  print(x)
  
  tmp <- read.delim(i, header = FALSE)$V1
  
  
  overlap[[x]] <- data.frame(dir = c("up", "dn"),
                             index = c(0,0))
  overlap[[x]][1,2] <- length(intersect(intersect(tmp, ca1_up), intersect(tmp, th_up)))
  overlap[[x]][2,2] <- length(intersect(intersect(tmp, ca1_dn), intersect(tmp, th_dn)))
  overlap[[x]] <- overlap[[x]] %>% reshape2::melt(id.vars="dir") %>% mutate(cell_type=x)
  
}
overlap <- do.call(rbind, overlap)

overlap <- overlap %>% filter(!cell_type %in% c("all", "unaccounted")) %>% 
  mutate(cell_type=factor(cell_type, levels=c("microglia", "astrocyte", "vascular", "oligodendrocyte", "opc", "neuronal")),
         dir=factor(dir,levels=c("up","dn"))) %>%
  select(cell_type, dir, value) %>% 
  reshape2::dcast(cell_type~dir, value.var="value")
rownames(overlap) <- overlap$cell_type
overlap <- overlap[,-1] %>% as.matrix()

#Get n genes at each timepoint
all_up <- read.delim("gene lists/all_up.txt", header=FALSE)$V1 %>%
  intersect(read.delim("gene lists final cell type/all_genes.txt", header = FALSE)$V1)
all_dn <- read.delim("gene lists/all_dn.txt", header=FALSE)$V1 %>%
  intersect(read.delim("gene lists final cell type/all_genes.txt", header = FALSE)$V1)

files <- c("results/CA1_DE_70_dpi.csv",
           "results/CA1_DE_90_dpi.csv",
           "results/CA1_DE_130_dpi.csv",
           "results/CA1_DE_150_dpi.csv",
           "results/TH_DE_70_dpi.csv",
           "results/TH_DE_90_dpi.csv",
           "results/TH_DE_130_dpi.csv",
           "results/TH_DE_150_dpi.csv")

cell_types <- c("gene lists final cell type/microglia_genes.txt",
                "gene lists final cell type/astrocyte_genes.txt",
                "gene lists final cell type/vascular_genes.txt",
                "gene lists final cell type/oligodendrocyte_genes.txt",
                "gene lists final cell type/opc_genes.txt",
                "gene lists final cell type/neuronal_genes.txt")

out <- list()
for (i in files) {
 print(i)
 tmp <- read.csv(i)
 for (j in cell_types) {
   x <- paste(strsplit(gsub("results/", "", i), "_")[[1]][1], strsplit(gsub("results/", "", i), "_")[[1]][3],
              strsplit(gsub("gene lists final cell type/", "", j), "_")[[1]][1], sep = "_")
   out[[x]]$region <- strsplit(gsub("results/", "", i), "_")[[1]][1]
   out[[x]]$dpi <- strsplit(gsub("results/", "", i), "_")[[1]][3]
   
   ct <- read.delim(j, header = FALSE)$V1
   out[[x]]$cell_type <- strsplit(gsub("gene lists final cell type/", "", j), "_")[[1]][1]
   
   out[[x]]$up <- tmp %>% filter(baseMean>15, log2FoldChange> 0.5, padj<0.05) %>% 
     pull(X) %>%
     intersect(all_up) %>%
     intersect(ct) %>%
     length()
   
   out[[x]]$dn <- tmp %>% filter(baseMean>15, log2FoldChange< -0.5, padj<0.05) %>% 
     pull(X) %>%
     intersect(all_dn) %>%
     intersect(ct) %>%
     length()
 }
}
out <- do.call(rbind, out) %>% 
  as_tibble() %>% 
  mutate(region=as.character(region),
         dpi=as.character(dpi),
         cell_type=factor(cell_type, levels=c("microglia", "astrocyte", "vascular", "oligodendrocyte", "opc", "neuronal")),
         up=as.numeric(up),
         dn=as.numeric(dn))
                                                   
ca1_up <- out %>% filter(region=="CA1") %>% select(dpi, cell_type, up) %>% 
  reshape2::dcast(cell_type~dpi, value.var = "up") %>% 
  select(cell_type, `70`, `90`, `130`, `150`)

ca1_dn <- out %>% filter(region=="CA1") %>% select(dpi, cell_type, dn) %>% 
  reshape2::dcast(cell_type~dpi, value.var = "dn") %>% 
  select(cell_type, `70`, `90`, `130`, `150`)

th_up <- out %>% filter(region=="TH") %>% select(dpi, cell_type, up) %>% 
  reshape2::dcast(cell_type~dpi, value.var = "up") %>% 
  select(cell_type, `70`, `90`, `130`, `150`)

th_dn <- out %>% filter(region=="TH") %>% select(dpi, cell_type, dn) %>% 
  reshape2::dcast(cell_type~dpi, value.var = "dn") %>% 
  select(cell_type, `70`, `90`, `130`, `150`)

rownames(ca1_up) <- ca1_up$cell_type
ca1_up <- ca1_up[,-1] %>% as.matrix()

rownames(ca1_dn) <- ca1_dn$cell_type
ca1_dn <- ca1_dn[,-1] %>% as.matrix()

rownames(th_up) <- th_up$cell_type
th_up <- th_up[,-1] %>% as.matrix()

rownames(th_dn) <- th_dn$cell_type
th_dn <- th_dn[,-1] %>% as.matrix()

cols_up <- circlize::colorRamp2(c(1,50,125,600), brewer.pal(9, "Oranges")[c(2,4,6,9)])
cols_dn <- circlize::colorRamp2(c(1,50,125,600), brewer.pal(9, "Purples")[c(2,4,6,9)])
cols_overlap <- circlize::colorRamp2(c(1,7,25,100), brewer.pal(9, "Greys")[c(2,4,6,9)])

gns <- brewer.pal(9, "Greens")
hdpi <- HeatmapAnnotation(dpi=as.factor(c(70, 90, 130, 150)), show_annotation_name = FALSE,
                          col=list(dpi=c(`70`=gns[2], `90`=gns[4], `130`=gns[7], `150`=gns[9])))
hdir <- HeatmapAnnotation(direction=as.factor(c("up", "dn")), show_annotation_name = FALSE,
                          col=list(direction=c(`up`=brewer.pal(9, "Oranges")[8], `dn`=brewer.pal(9, "Purples")[8])))

pdf("plots/n_genes_cell_types.pdf", width=10, height = 3)
ht_list <- Heatmap(ca1_up, col = cols_up, cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left",
                   name = "n genes increased", column_title = "CA1", top_annotation = hdpi,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if(ca1_up[i, j] > 0)
                       grid.text(sprintf("%.0f", ca1_up[i, j]), x, y, gp = gpar(fontsize = 10))
                   }) +
  Heatmap(th_up, col = cols_up, cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left",
          show_heatmap_legend = FALSE, column_title = "thalamus", top_annotation = hdpi,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(th_up[i, j] > 0)
              grid.text(sprintf("%.0f", th_up[i, j]), x, y, gp = gpar(fontsize = 10))
          }) +
  Heatmap(ca1_dn, col = cols_dn, cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left",
          name = "n genes decreased", column_title = "CA1", top_annotation = hdpi,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(ca1_dn[i, j] > 0)
              grid.text(sprintf("%.0f", ca1_dn[i, j]), x, y, gp = gpar(fontsize = 10))
          }) +
  Heatmap(th_dn, col = cols_dn, cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left",
          show_heatmap_legend = FALSE, column_title = "thalamus", top_annotation = hdpi,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(th_dn[i, j] > 0)
              grid.text(sprintf("%.0f", th_dn[i, j]), x, y, gp = gpar(fontsize = 10))
          }) +
  Heatmap(overlap, col = cols_overlap, cluster_columns = FALSE, cluster_rows = FALSE, row_names_side = "left",
          name = "n genes overlapping", column_title = "Overlapping", top_annotation = hdir,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(overlap[i, j] > 0)
              grid.text(sprintf("%.0f", overlap[i, j]), x, y, gp = gpar(fontsize = 10))
          })
draw(ht_list)
dev.off()

