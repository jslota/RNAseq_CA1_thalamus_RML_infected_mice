library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)

#Thalamus neuron genes
markers <- read.csv("cell type databases/thalamus_neuron_categories.csv")
rownames(markers) <- markers$gene

tmp <- read.delim("gene lists region specific/all/neuronal_th_dn.txt", header = FALSE)$V1
tmp %>% setdiff(markers$gene) %>% length()
tmp %>% writeClipboard()

#fold-changes
fc_dt <- read.csv("all_fold_change_values.csv", row.names = 1) %>% 
  select(X, TH_DE_70_dpi, TH_DE_90_dpi, TH_DE_130_dpi, TH_DE_150_dpi) %>% 
  filter(X %in% markers$gene)

hmp_dt <- fc_dt
rownames(hmp_dt) <- hmp_dt$X
hmp_dt <- hmp_dt[,-1]
hmp_dt <- as.matrix(hmp_dt)
markers <- markers[rownames(hmp_dt),]

gns <- brewer.pal(9, "Greens")
ann_dt <- data.frame(row.names = colnames(hmp_dt),
                     dpi = c(70, 90, 130, 150))
ha <- HeatmapAnnotation(dpi=as.factor(ann_dt$dpi),
                        col=list(dpi=c(`70`=gns[2], `90`=gns[4], `130`=gns[7], `150`=gns[9])))

FC <- circlize::colorRamp2(c(min(hmp_dt), -3, -2.25, -1.5, -0.75, 0, 0.75, 1.5, 2.25, 3, max(hmp_dt)), rev(brewer.pal(11, "PuOr")))

table(markers$group)

tmp <- markers %>% filter(group == "synaptic transmission") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht1 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
                width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "Axon guidance") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht2 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "potassium transport") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht3 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "response to calcium") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht4 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "actin cytoskeleton") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht5 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "neuron death") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht6 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "calcium transport") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht7 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)


pdf("plots/neuron analysis/th_synaptic_transmission.pdf", width = 3.5, height = 7)
draw(ht1)
dev.off()

pdf("plots/neuron analysis/th_axon_guidance.pdf", width = 3.5, height = 6)
draw(ht2)
dev.off()

pdf("plots/neuron analysis/th_potassium_transport.pdf", width = 3.5, height = 5)
draw(ht3)
dev.off()

pdf("plots/neuron analysis/th_respone_to_calcium.pdf", width = 3.5, height = 5)
draw(ht4)
dev.off()

pdf("plots/neuron analysis/th_actin_cytoskeleton.pdf", width = 3.5, height = 5)
draw(ht5)
dev.off()

pdf("plots/neuron analysis/th_neuron_death.pdf", width = 3.5, height = 3.5)
draw(ht6)
dev.off()

pdf("plots/neuron analysis/th_calcium_transport.pdf", width = 3.5, height = 3.5)
draw(ht7)
dev.off()


#CA1 neuron genes
markers <- read.csv("cell type databases/ca1_neuron_categories.csv")
rownames(markers) <- markers$gene

tmp <- union(read.delim("gene lists region specific/ca1_neuron_and_enriched_dn.txt", header = FALSE)$V1,
             read.delim("gene lists region specific/ca1_neuron_and_enriched_up.txt", header = FALSE)$V1)

#fold-changes
fc_dt <- read.csv("all_fold_change_values.csv", row.names = 1) %>% 
  select(X, CA1_DE_70_dpi, CA1_DE_90_dpi, CA1_DE_130_dpi, CA1_DE_150_dpi) %>% 
  filter(X %in% markers$gene)

hmp_dt <- fc_dt
rownames(hmp_dt) <- hmp_dt$X
hmp_dt <- hmp_dt[,-1]
hmp_dt <- as.matrix(hmp_dt)
markers <- markers[rownames(hmp_dt),]

table(markers$group)


tmp <- markers %>% filter(group == "synaptic plasticity") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht1 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "metabolism") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht2 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "stress response") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht3 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "neurogenesis") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht4 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "neuroendrocrine") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht5 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "dendritic spines") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht6 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)

tmp <- markers %>% filter(group == "cilia") %>% pull(gene)
tmp <- hmp_dt[tmp,]
ht7 <- Heatmap(tmp, col = FC, name="Log2(FC)", cluster_columns = FALSE, top_annotation = ha,
               width = ncol(tmp)*unit(5, "mm"), height = nrow(tmp)*unit(5, "mm"), show_column_names = FALSE)


pdf("plots/neuron analysis/ca1_synaptic_plasticity.pdf", width = 3.5, height = 5)
draw(ht1)
dev.off()

pdf("plots/neuron analysis/ca1_metabolism.pdf", width = 3.5, height = 5)
draw(ht2)
dev.off()

pdf("plots/neuron analysis/ca1_stress_response.pdf", width = 3.5, height = 5)
draw(ht3)
dev.off()

pdf("plots/neuron analysis/ca1_neurogenesis.pdf", width = 3.5, height = 5)
draw(ht4)
dev.off()

pdf("plots/neuron analysis/ca1_neuroendocrine.pdf", width = 3.5, height = 5)
draw(ht5)
dev.off()

pdf("plots/neuron analysis/ca1_dendritic_spines.pdf", width = 3.5, height = 5)
draw(ht6)
dev.off()

pdf("plots/neuron analysis/ca1_cilia.pdf", width = 3.5, height = 5)
draw(ht7)
dev.off()