library(openxlsx)
library(dplyr)

load("C:/Projects/Manish/analysis/analysis_19092023/cluster_geneAnnotations_20231220.RData")
#load("C:/Projects/Manish/analysis/analysis/cluster_geneAnnotations.RData")

allData <- do.call(rbind,clstMarker)
rownames(allData) <- 1:nrow(allData)

kkOrthoUpdated <- read.table("C:/Projects/Manish/tables/kk_structure_orthology_marker_genes.txt",
                             sep = "\t",header = TRUE)
leafAtlasUpdated <- read.table("C:/Projects/Manish/tables/leafatlas_structure_orthology_marker_genes.txt",
                               sep = "\t",header = TRUE)
allOrthologs <- read.table("C:/Projects/Manish/tables/structural_orthologs_cur_athal.geneids.tsv",
                           sep = "\t",header = TRUE)

allMatchOrtho <- merge(allData,allOrthologs,
                       by.x = "GeneID",by.y = "currot_id")

fName <- "C:/Projects/Manish/analysis/analysis_19092023/scRNAseqClusterGene_allOrthologs_overlap.tsv"

finalDf <- allMatchOrtho[,setdiff(colnames(allMatchOrtho),c("gene_name_KK","gene_name_LeafAtlas"))]

write.table(x = finalDf, 
            sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,
            file = fName)


