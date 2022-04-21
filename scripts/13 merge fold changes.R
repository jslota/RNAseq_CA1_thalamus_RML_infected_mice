library(dplyr)

genes <- read.delim("gene lists final cell type/all_genes.txt", header = FALSE)$V1
out <- data.frame(X = genes)
cnames <- "X"
for (i in Sys.glob("results/*")) {
  print(i)
  cnames <- c(cnames, gsub(".csv", "", gsub("results/", "", i)))
  tmp <- read.csv(i) %>% select(X, log2FoldChange)
  out <- left_join(out, tmp, by="X")
}
colnames(out) <- cnames
out[is.na(out)] <- 0

write.csv(out, "all_fold_change_values.csv")
