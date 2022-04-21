###Normalize RNAseq data
#2021-09-09
#Jessy Slota

library(DESeq2)
###Normalize CA1 and TH data
#ca1
#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_CA1",]
raw_dat <- raw_dat[,rownames(smpl_dat)]

dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Prions)
dds <- DESeq(dds)
norm_dat <- assay(vst(dds))
write.csv(norm_dat, "normalized data/CA1_normalized_counts.csv")

#th
#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_TH",]
raw_dat <- raw_dat[,rownames(smpl_dat)]

#normalize data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Prions)
dds <- DESeq(dds)
norm_dat <- assay(vst(dds))
write.csv(norm_dat, "normalized data/TH_normalized_counts.csv")

