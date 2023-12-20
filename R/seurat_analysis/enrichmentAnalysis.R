library(dplyr)
library(STRINGdb)
library(openxlsx)

load("C:/Projects/Manish/analysis/analysis_19092023/clusterObj_umap_res0.6_20230919.RData")
clstMarkerInfo <- "C:/Projects/Manish/analysis/analysis_19092023/cluster_markers_res0.5_with_GeneNames_20231220.xlsx"

allGeneMapping <- read.table("C:/Projects/Manish/tables/structural_orthologs_cur_athal.geneids.tsv",
                             sep = "\t",header = T)

sheetNms <- getSheetNames(clstMarkerInfo)

wb = createWorkbook()
clstMarker <- list()

wbEnr = createWorkbook()
clstEnr <- list()

string_db <- STRINGdb$new( version="11.5", species=3702,
                           score_threshold=400, input_directory="")

for(i in 1:length(sheetNms)){
  dat <- read.xlsx(xlsxFile = clstMarkerInfo,sheet = sheetNms[i])
  gNames <- merge(dat,allGeneMapping,by.x = "GeneID",by.y = "currot_id")
  geneDiff <- gNames[which(gNames$p_val_adj < 0.05 & gNames$avg_log2FC > 0),]
  gDiff_mapped <- string_db$map( geneDiff,
                                 "athal_id", removeUnmappedRows = TRUE )
  diffSub <- gDiff_mapped[,c("athal_id","p_val_adj","avg_log2FC")]
  colnames(diffSub) <- c("gene","pvalue","logFC")
  enr <- string_db$get_enrichment(diffSub)
  clstMarker[[i]] <- gNames
  clstEnr[[i]] <- diffSub
  addWorksheet(wb, sheetNms[i])
  addWorksheet(wbEnr, paste0("enrichment",sheetNms[i]))
  writeData(wb = wb, sheet = sheetNms[i], x = gNames)
  writeData(wb = wbEnr, sheet = paste0("enrichment",sheetNms[i]), x = enr)
}

wbName <- paste0(getwd(),"/analysis/cluster_markers_res0.5_with_GeneNames_annotations_",
                 gsub("-","",Sys.Date()),".xlsx")
wbEnrName <- paste0(getwd(),"/analysis/cluster_Enrichment_res0.5_with_GeneNames_annotations_",
                    gsub("-","",Sys.Date()),".xlsx")

saveWorkbook(wb,file = wbName)
saveWorkbook(wbEnr, file = wbEnrName)

names(clstMarker) <- paste0("cluster",unique(int.combinedClst$seurat_clusters))

clstMNme <- paste0(getwd(),"/analysis/cluster_geneAnnotations_",
                   gsub("-","",Sys.Date()),".RData")
save(clstMarker, file = clstMNme)