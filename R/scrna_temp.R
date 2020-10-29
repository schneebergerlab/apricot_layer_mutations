library(Seurat)
library(data.table)
library(dplyr)
df <- fread("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/ryu_et_al_plant_physio_2019/GSE123013_5way_merge_raw.tsv",header = TRUE)
df <- data.frame(df)
wt2 <- df %>% select(c('V1', grep("WT2", colnames(df))))
rownames(wt2) <- wt2$V1
wt2$V1 <- NULL
rownames(wt2) <- gsub("_", "-", rownames(wt2))
pbmc <- CreateSeuratObject(counts = as.matrix(wt2), project = "ryu", min.cells = 3, min.features = 200)
