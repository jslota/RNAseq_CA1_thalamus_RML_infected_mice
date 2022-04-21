###Mapping QC
###2021-08-18
###Jessy Slota

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)

#mapping statistics
map_stats <- read.csv("raw data/Mouse_RNAseq_mapping_stats.csv")
rownames(map_stats) <- map_stats$Sample
sample_info <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1)
map_stats <- map_stats %>% filter(Sample %in% rownames(sample_info))


#order samples
s1 <- map_stats[rownames(sample_info[sample_info$Experiment == "RML_CA1",]),]
s1 <- s1[order(s1$Mapped, decreasing = TRUE),]
s2 <- map_stats[rownames(sample_info[sample_info$Experiment == "RML_TH",]),]
s2 <- s2[order(s2$Mapped, decreasing = TRUE),]

samples_orders <- c(rownames(s1), rownames(s2))

#format for plotting
map_stats$Sample <- factor(map_stats$Sample, levels = rev(samples_orders))
map_stats$Experiment <- sample_info$Experiment

map_plot_data <- melt(map_stats, id.vars = c("Sample","Experiment"))

#make plot
ggplot(data = map_plot_data, mapping = aes(x = Sample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  ylab("Raw read counts") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = brewer.pal(5, "Accent")) +
  coord_flip() +
  theme_classic() +
  theme(text = element_text(size = 9),
        legend.title = element_blank(), axis.title.y = element_blank())
ggsave("plots/mapping.svg", width = 6, height = 6)
