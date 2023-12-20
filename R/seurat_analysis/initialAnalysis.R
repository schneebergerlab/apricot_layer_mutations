# This script performs the downstream analysisof scRNA-seq data
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(openxlsx)
library(dplyr)

set.seed(10)
#### function to filter artefactual cells ####

filterData <- function(sobj){
  countHIGH=quantile(sobj$nCount_RNA, prob=0.95) 
  # set lower limit for nFeature_RNA (number of genes detected per cell)
  featureLOW <- 500
  # set lower limit for nCounts_RNA (number of UMIs per cell)
  countLOW <- 1000
  fObj <- subset(sobj, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)
  return(fObj)
}

# sample names
smplName <- c("mut_11_1","mut_15","wt1","wt19")
# genotype/condition
cndList <- c("MUT11","MUT15",rep("WT",2))

smplList <- sapply(1:length(smplName),function(x){
  pathVal <- paste0("Data/newV3/",smplName[x],"/")
  smpl <- Read10X(data.dir = pathVal)
  smplObj <- CreateSeuratObject(counts = smpl, project = smplName[x], 
                                min.cells = 3, min.features = 250)
  smplObj[["group"]] <- cndList[x]
  smplObj
})

names(smplList) <- smplName

VlnPlot(smplList[[1]], features = c("nFeature_RNA","nCount_RNA"), 
        ncol = 3) 

# # Set up WT object
############ this did not work ############ 
# allWT <- merge(smplList[[3]], 
#                y = c(smplList[[4]]), 
#                add.cell.ids = smplName[3:4])
# 
# VlnPlot(allWT, features = c("nFeature_RNA","nCount_RNA"), 
#         ncol = 2)    
# allWT$group

################################################ 

ctrlWT1 <- filterData(smplList[[3]])
ctrlWT19 <- filterData(smplList[[4]])

# no hardcoded threshold
#ctrlWT1 <- subset(smplList[[3]], subset = nFeature_RNA > 700 & nFeature_RNA < 2500)
#ctrlWT19 <- subset(smplList[[4]], subset = nFeature_RNA > 700 & nFeature_RNA < 2500)

# set up Mut11
mut11 <- filterData(smplList[[1]])

# no hardcoded threshold
#mut11 <- subset(smplList[[1]], subset = nFeature_RNA > 1000 & nFeature_RNA < 2000)

# set up Mut15
mut15 <- filterData(smplList[[2]])

# no hardcoded threshold
#mut15 <- subset(smplList[[2]], subset = nFeature_RNA > 500 & nFeature_RNA < 2000)

VlnPlot(ctrlWT1, features = c("nFeature_RNA","nCount_RNA"), 
          ncol = 2) 
ifnb.list <- list(ctrlWT1, ctrlWT19, mut11,mut15)
save(ifnb.list,file = "seuratObj_filtered_beforeNormalization.RData")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

int.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
int.combined <- IntegrateData(anchorset = int.anchors)

DefaultAssay(int.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
int.combined <- ScaleData(int.combined, verbose = FALSE)
save(int.combined,file = "seuratObj_integrated_scaled.RData")

int.combinedPCA <- RunPCA(int.combined, npcs = 30, verbose = FALSE)
int.combinedUMAP <- RunUMAP(int.combinedPCA, reduction = "pca", dims = 1:30)
int.combinedNbh <- FindNeighbors(int.combinedUMAP, reduction = "pca", dims = 1:30)
int.combinedClst <- FindClusters(int.combinedNbh, resolution = 0.5)
save(int.combinedClst,file = "clusterObj_umap_res0.5.RData")

# plot umap
umapSampleSplt <- DimPlot(int.combinedClst, reduction = "umap", split.by = "orig.ident",
                          label = TRUE)
# add transparency to dots
#umapSampleSplt$layers[[1]]$aes_params$alpha =  .7
ggsave("UMAP_splitBySample.pdf", plot = umapSampleSplt,
       width = 7.64, height = 4.70, units = "in")

umapAll <- DimPlot(int.combinedClst, reduction = "umap",
                   label = TRUE)
ggsave("UMAP_noSplit.pdf", plot = umapAll,
       width = 4.74, height = 4.70, units = "in")

wb = createWorkbook()
clstMarker <- list()
for (i in 1:length(unique(int.combinedClst$seurat_clusters))) {
  clstInfo <- FindMarkers(int.combinedClst,ident.1 = unique(int.combinedClst$seurat_clusters)[i],
                          min.pct = 0.25)
  geneID <- rownames(clstInfo)
  clstInfo["GeneID"] <- geneID
  clstInfo["cluster"] <- paste0("cluster",unique(int.combinedClst$seurat_clusters)[i])
  clstMarker[[i]] <- clstInfo[which(clstInfo$p_val_adj < 0.05),]
  sheet_name = paste('cluster', unique(int.combinedClst$seurat_clusters)[i])
  addWorksheet(wb, sheet_name)
  writeData(wb = wb, sheet = sheet_name, x = clstInfo)
}

wbName <- paste0(getwd(),"/analysis/analysis_19092023/cluster_markers_res0.5_with_GeneNames_",
                 gsub("-","",Sys.Date()),".xlsx")

saveWorkbook(wb,file = wbName,overwrite = T)
names(clstMarker) <- paste0("cluster",unique(int.combinedClst$seurat_clusters))

clstMNme <- paste0(getwd(),"/analysis/analysis_19092023/cluster_geneAnnotations_",
                   gsub("-","",Sys.Date()),".RData")
save(clstMarker, file = clstMNme)
