###Allen brain atlas GSEA
#2021-08-24
#Jessy Slota

library(DESeq2)
library(readxl)
library(dplyr)

###Clean up gene sets files
tmp <- ActivePathways::read.GMT("thalamus allen gsea/gene_set_library_up_crisp.gmt")
names(tmp) <- gsub(",", "", gsub(" ", "_", names(tmp)))
names(tmp) <- toupper(names(tmp))
tmp <- tmp[!duplicated(names(tmp))]
for (i in names(tmp)) {
  tmp[[i]]$id <- i
}
ActivePathways::write.GMT(tmp, "thalamus allen gsea/Allen_cleaned_gene_sets.gmt")
rm(tmp)

#load dataset
raw_dat <- read.csv("raw data/Mouse_RNAseq_read_counts.csv", row.names = 1, header = TRUE)
smpl_dat <- read.csv("raw data/Mouse_RNAseq_sample_info.csv", row.names = 1, header = TRUE)

head(smpl_dat)
smpl_dat <- smpl_dat[,c(1:4)]
smpl_dat <- smpl_dat[smpl_dat$Experiment == "RML_TH",]
raw_dat <- raw_dat[,rownames(smpl_dat)]

#Get normalized data
dds <- DESeqDataSetFromMatrix(raw_dat, smpl_dat, ~Prions)
dds <- DESeq(dds)
norm_dat <- vst(dds)

##Get sample specific read counts for GSEA
gsea <- assay(norm_dat)
head(gsea)
summary(rowMeans(gsea) > 1.05*min(gsea))
gsea <- as.matrix(gsea[rowMeans(gsea) > 1.05*min(gsea),])
gsea <- gsea[rowVars(gsea) > mean(rowVars(gsea)),]

gsea <- gsea - rowMeans(gsea)

#write .rnk files for gsea
for (i in colnames(gsea)) {
  write.table(data.frame(gene = rownames(gsea), expr = gsea[,i]), paste0("thalamus allen gsea/", i, ".rnk"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

#Load all GSEA results and merge into one file
res <- list()

for (i in Sys.glob("thalamus allen gsea/*/gsea_report*.tsv")) {
  print(gsub("thalamus allen gsea/*", "", gsub("\\..*", "", i))) #Make sure we are getting correct sample names
  res[[i]] <- read.delim(i) %>% #load files, clean up with dplyr, and add to list
    select(NAME, SIZE, ES, NES, NOM.p.val, FDR.q.val, FWER.p.val) %>%
    mutate(Sample = gsub("thalamus allen gsea/*", "", gsub("\\..*", "", i)))
}

res <- do.call(rbind, res) #make list into data frame
res[res$NES == "---",]$NES <- NA
res$NES <- as.numeric(res$NES)
write.csv(res, "thalamus allen gsea/Full_gsea_res.csv", row.names = FALSE)
