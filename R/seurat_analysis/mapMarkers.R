# this script maps marker genes with the results of scRNA-seq FindMarker
setwd("C:/Projects/Manish/")

allData <- read.table("C:/Projects/Manish/analysis/analysis_19092023/scRNAseqClusterGene_allOrthologs_overlap.tsv",
                      sep = "\t",header = TRUE)

# all annotations

## leafatlas2_epiderm_orthology_marker_genes ->  leafatlas epidermis abxial adxial sequence ortholog
## leafatlas2_epiderm_structure_orthology_marker_genes -> leafatlas epidermis abxial adxial structure ortholog
## leafatlas2_layer_orthology_marker_genes -> leafatlas sequence ortholog (markers found using FindMarkers)
## leafatlas2_layer_structure_orthology_marker_genes -> leafatlas structure ortholog (markers found using FindMarkers)
## leafatlas2_orthology_marker_genes -> leafatlas sequence ortholog (markers found using old publications)
## leafatlas2_structure_orthology_marker_genes -> leafatlas structure ortholog (markers found using old publications)

leafEpOrthoSq <- read.table("tables/leafAtlasUpdated/leafatlas2_epiderm_orthology_marker_genes.txt", 
                            sep = "\t", header = TRUE)
leafEpOrthoSq["Source"] <- "leafatlas2_epiderm_orthology_marker_genes"

leafEpOrthoSt <- read.table("tables/leafAtlasUpdated/leafatlas2_epiderm_structure_orthology_marker_genes.txt", 
                            sep = "\t", header = TRUE)
leafEpOrthoSt["Source"] <- "leafatlas2_epiderm_structure_orthology_marker_genes"

leafLOrthoSq <- read.table("tables/leafAtlasUpdated/leafatlas2_layer_orthology_marker_genes.txt", 
                           sep = "\t", header = TRUE)
leafLOrthoSq["Source"] <- "leafatlas2_layer_orthology_marker_genes"

leafLOrthoSt <- read.table("tables/leafAtlasUpdated/leafatlas2_layer_structure_orthology_marker_genes.txt", 
                           sep = "\t", header = TRUE)
leafLOrthoSt["Source"] <- "leafatlas2_layer_structure_orthology_marker_genes"

leafOrthoSq <- read.table("tables/leafAtlasUpdated/leafatlas2_orthology_marker_genes.txt", 
                          sep = "\t", header = TRUE)
leafOrthoSq["Source"] <- "leafatlas2_orthology_marker_genes"

leafOrthoSt <- read.table("tables/leafAtlasUpdated/leafatlas2_structure_orthology_marker_genes.txt", 
                          sep = "\t", header = TRUE)
leafOrthoSt["Source"] <- "leafatlas2_structure_orthology_marker_genes"

kkSt <- read.table("tables/kk_structure_orthology_marker_genes.txt",
                   sep = "\t", header = TRUE)
kkSt["Source"] <- "kk_structure_orthology_marker_genes"
colnames(kkSt) <- colnames(leafOrthoSt)

kkSq <- read.table("tables/kk_orthology_marker_genes.txt",
                   sep = "\t", header = TRUE)
kkSq["Source"] <- "kk_orthology_marker_genes"
colnames(kkSq) <- colnames(leafOrthoSq)

leaftAtV1St <- read.table("tables/leafatlas_structure_orthology_marker_genes.txt",
                          sep = "\t", header = TRUE)
leaftAtV1St["Source"] <- "leafatlas_structure_orthology_marker_genes"
leaftAtV1Sq <- openxlsx::read.xlsx("tables/leafatlas_orthology_marker_genes.xlsx") 
leaftAtV1Sq["Source"] <- "leafatlas_orthology_marker_genes"


allAnnoLeafAtSt <- rbind(leafEpOrthoSt, leafLOrthoSt, 
                         leafOrthoSt,kkSt,leaftAtV1St)
allAnnoLeafAtSq <- rbind(leafEpOrthoSq, leafLOrthoSq, 
                         leafOrthoSq,kkSq,leaftAtV1Sq)
allAnnoLeafAtSq["reciprocal_match"] <- "Not applicable"
colnames(allAnnoLeafAtSq)[3] <- "currot_id"
allAnnoLeafAtSq <- allAnnoLeafAtSq[,colnames(allAnnoLeafAtSt)]


allAnnoLeafAt <- rbind(allAnnoLeafAtSq, allAnnoLeafAtSt)

finalV1Anno <- merge(allData,allAnnoLeafAt,by.x = "GeneID",by.y = "currot_id")  
openxlsx::write.xlsx(finalV1Anno, file = "C:/Projects/Manish/analysis/analysis_19092023/scRNAseqClusterGene_allLeafAtlasAnno_updated_19092023.xlsx")

