library(Rsubread)
library(rtracklayer)
library(doParallel)
library(foreach)
library(Seurat)
library(data.table)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(cowplot)

gff <- import.gff3('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.protein_coding.3utr.gff3')
gffgenes <- gff[gff@elementMetadata$type == 'gene']
garb <- seqnames(gffgenes)
chrs <- rep(garb@values, garb@lengths)
spos <- start(gffgenes)
epos <- end(gffgenes)
strnd <- strand(gffgenes)
geneid <- mcols(gffgenes)[,"ID"]
ann <- data.frame('GeneID'=geneid, 'Chr'=chrs, 'Start'=spos, 'End'=epos, 'Strand'=strnd)

SAMPLES=c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/'
samcnt <- list()
cl <- makeCluster(60)
registerDoParallel(cl)
for(sample in SAMPLES){
  bcs <-  sort(read.table(paste(CWD, sample, 'good_bcs_list.txt', sep='/'))$V1)
  bccnt <- foreach(bc=bcs, .packages = c('Rsubread', 'dplyr')) %dopar% {
    featureCounts(files=paste(CWD, sample, paste0(bc,'-1.bam'), sep='/'), annot.ext=ann)$counts
  }
  names(bccnt) <- bcs
  samcnt[[sample]] <- bccnt
}
stopCluster(cl)
# saveRDS(samcnt, paste0(CWD,'samcnt.rds'))

samcnt <- readRDS(paste0(CWD,'samcnt.rds'))
for(sample in SAMPLES){
  samcnt[[sample]] <- data.frame(samcnt[[sample]])
}

scplt <- function(pop, seuobj){
  # Take a processed and clustered seuobj and generate plots for the marker genes
  pdf(paste0(CWD, pop, '.pdf'), height=7, width=7)

  # Plot the clusters
  plot(DimPlot(seuobj, reduction = "umap"))

  # Get all Marker genes
  seuobj.markers <- FindAllMarkers(seuobj, only.pos=TRUE, min.pct = 0.25, logfc.threshold = 1, test.use="MAST")
  garb <- seuobj.markers %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = avg_log2FC)

  # Get plots for all markers
  vlnplts <- VlnPlot(seuobj, features=garb$gene, combine=FALSE)
  ftplts <- FeaturePlot(seuobj, features = garb$gene, keep.scale='all', combine=FALSE)

  # Save plots for individual markers
  for(i in 1:nrow(garb)){
    geneid <- ftplts[[i]]$labels$title
    if(geneid != vlnplts[[i]]$labels$title)   print(paste('Incorrect Genes together', i))
    genedf <- garb %>% filter(gene == geneid)
    plt1 <- ftplts[[i]] / vlnplts[[i]]
    plt1 <- plt1 + plot_annotation(
      title = pop,
      subtitle = paste0("Ave. Log2FC: ", round(genedf$avg_log2FC, 2),"; Cluster cells %: ", genedf$pct.1, "; Outside cells %: ", genedf$pct.2)
    )
    plot(plt1)
  }
  dev.off()
  write.table(garb, paste0(CWD, pop,"_DG.txt"), sep="\t", quote = FALSE, row.names = FALSE)
}

# Run the pipeline for WT_1
pop <- 'WT_1'
df <- samcnt[[pop]]
seuobj <- CreateSeuratObject(counts = df, project = pop, min.cells = 3, min.features = 200)
#VlnPlot(seuobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(seuobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
seuobj <- subset(seuobj, subset = nFeature_RNA > 200 & nFeature_RNA < 3500)
seuobj <- NormalizeData(seuobj)
seuobj <- FindVariableFeatures(seuobj, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(seuobj), 10)
#plot1 <- VariableFeaturePlot(seuobj)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
all.genes <- rownames(seuobj)
seuobj <- ScaleData(seuobj, features = all.genes)
seuobj <- RunPCA(seuobj, features = VariableFeatures(object = seuobj))
#seuobj <- JackStraw(seuobj, num.replicate = 100, dims=50)
#seuobj <- ScoreJackStraw(seuobj, dims = 1:50)
#JackStrawPlot(seuobj, dims = 1:50)
# ndims selected = 30
seuobj <- FindNeighbors(seuobj, dims = 1:30)
seuobj <- FindClusters(seuobj, resolution = 0.5)
seuobj <- RunUMAP(seuobj, dims = 1:30)
seuobj <- RunTSNE(seuobj, dims = 1:30)
wt1seuobj <- seuobj
scplt(pop, seuobj)


# Run the pipeline for WT_19
pop <- 'WT_19'
df2 <- samcnt[[pop]]
seuobj <- CreateSeuratObject(counts = df2, project = pop, min.cells = 3, min.features = 200)
#VlnPlot(seuobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
seuobj <- subset(seuobj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)
seuobj <- NormalizeData(seuobj)
seuobj <- FindVariableFeatures(seuobj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seuobj)
seuobj <- ScaleData(seuobj, features = all.genes)
seuobj <- Seurat::RunPCA(seuobj, features = VariableFeatures(object = seuobj))
#seuobj <- JackStraw(seuobj, num.replicate = 100, dims=50)
#seuobj <- ScoreJackStraw(seuobj, dims = 1:50)
#JackStrawPlot(seuobj, dims = 1:50)
#ElbowPlot(seuobj,ndims = 50)
# ndims selected = 30
seuobj <- FindNeighbors(seuobj, dims = 1:31)
seuobj <- FindClusters(seuobj, resolution = 0.5)
seuobj <- RunUMAP(seuobj, dims = 1:30)
seuobj <- RunTSNE(seuobj, dims = 1:30)
wt19seuobj <- seuobj
scplt(pop, seuobj)

# Run the pipeline for MUT_11_1
pop <- 'MUT_11_1'
df2 <- samcnt[[pop]]
seuobj <- CreateSeuratObject(counts = df2, project = pop, min.cells = 3, min.features = 200)
#VlnPlot(seuobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
seuobj <- subset(seuobj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
seuobj <- NormalizeData(seuobj)
seuobj <- FindVariableFeatures(seuobj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seuobj)
seuobj <- ScaleData(seuobj, features = all.genes)
seuobj <- Seurat::RunPCA(seuobj, features = VariableFeatures(object = seuobj))
#seuobj <- JackStraw(seuobj, num.replicate = 100, dims=50)
#seuobj <- ScoreJackStraw(seuobj, dims = 1:50)
#JackStrawPlot(seuobj, dims = 1:50)
#ElbowPlot(seuobj,ndims = 50)
# ndims selected = 31
seuobj <- FindNeighbors(seuobj, dims = 1:31)
seuobj <- FindClusters(seuobj, resolution = 0.5)
seuobj <- RunUMAP(seuobj, dims = 1:30)
seuobj <- RunTSNE(seuobj, dims = 1:30)
mut11seuobj <- seuobj
scplt(pop, seuobj)

# Run the pipeline for MUT_15
pop <- 'MUT_15'
df2 <- samcnt[[pop]]
seuobj <- CreateSeuratObject(counts = df2, project = pop, min.cells = 3, min.features = 200)
#VlnPlot(seuobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
seuobj <- subset(seuobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
seuobj <- NormalizeData(seuobj)
seuobj <- FindVariableFeatures(seuobj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seuobj)
seuobj <- ScaleData(seuobj, features = all.genes)
seuobj <- Seurat::RunPCA(seuobj, features = VariableFeatures(object = seuobj))
#seuobj <- JackStraw(seuobj, num.replicate = 100, dims=50)
#seuobj <- ScoreJackStraw(seuobj, dims = 1:50)
#JackStrawPlot(seuobj, dims = 1:50)
#ElbowPlot(seuobj,ndims = 50)
# ndims selected = 33
seuobj <- FindNeighbors(seuobj, dims = 1:33)
seuobj <- FindClusters(seuobj, resolution = 0.5)
seuobj <- RunUMAP(seuobj, dims = 1:33)
seuobj <- RunTSNE(seuobj, dims = 1:33)
mut15seuobj <- seuobj
scplt(pop, seuobj)
saveRDS(list("wt1seuobj"=wt1seuobj,
             "wt19seuobj"=wt19seuobj,
             "mut11seuobj"=mut11seuobj,
             "mut15seuobj"=mut15seuobj),
        file=paste0(CWD,"sample_seuobj.rds"))

# Differential genes in when analysing multiple runs simultaneously
sccombplt <- function(pop, seucomb){
  pdf(paste0(CWD, pop, '.pdf'), height=10, width=10)

  # Plot the clusters
  plot(DimPlot(seucomb, reduction = "umap")/DimPlot(seucomb, reduction = "umap", split.by = "orig.ident"))

  genes <- data.frame()
  for(i in sort(unique(seucomb$seurat_clusters))){
    seu.markers <- FindConservedMarkers(seucomb, ident.1 = i, grouping.var = "orig.ident",  logfc.threshold = 1, test.use="MAST", verbose = FALSE)
    varnames <- colnames(seu.markers)[grep("avg_log2FC", colnames(seu.markers))]
    garb <- seu.markers %>%
      filter(max_pval < 0.05) %>%
      filter(!!sym(varnames[1]) > 1, !!sym(varnames[2]) > 1) %>%
      #filter(MUT_11_1_avg_log2FC > 1, MUT_11_1_avg_log2FC>1) %>%
      top_n(n = 3, wt = !!sym(varnames[1])+!!sym(varnames[2]))
    head(garb)
    garb['gene'] <- rownames(garb)
    genes <- bind_rows(genes, garb)
    if(nrow(garb) == 0) next
    for(j in 1:nrow(garb)){
      #print(j)
      geneid <- rownames(garb)[j]
      genedf <- garb %>% filter(gene == geneid)
      plt <- FeaturePlot(seucomb, features = geneid, split.by='orig.ident') / VlnPlot(seucomb, features = geneid, split.by='orig.ident')
      parents <- unique(seucomb$orig.ident)
      subt <- ""
      for(par in parents){
        subt <- paste0(subt, par, " - Ave.Log2FC: ", round(genedf[,paste0(par,"_avg_log2FC")], 2),"; Cluster cells %: ", genedf[,paste0(par, "_pct.1")], "; Outside cells %: ", genedf[,paste0(par, "_pct.2")], "\n")
      }
      plt <- plt + plot_annotation(
        title = paste(pop, "Cluster: ", i),
        subtitle = subt
      )
      plot(plt)
    }
  }
  dev.off()
  write.table(genes, paste0(CWD, pop,"_DG.txt"), sep="\t", quote = FALSE, row.names = FALSE)
}

# Analyse both wild-types together
wt1seu <- CreateSeuratObject(counts = samcnt[['WT_1']], project = 'WT_1', min.cells = 3, min.features = 200)
wt1seu <- subset(wt1seu, subset = nFeature_RNA > 200 & nFeature_RNA < 3500)
wt19seu <- CreateSeuratObject(counts = samcnt[['WT_19']], project = 'WT_19', min.cells = 3, min.features = 200)
wt19seu <- subset(wt19seu, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)
wtseu.list <- list('wt1seu'=wt1seu, 'wt19seu'=wt19seu)

wtseu.list <- lapply(X = wtseu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = wtseu.list, nfeatures = 3000)
wt.anchors <- FindIntegrationAnchors(object.list = wtseu.list, anchor.features = features)
wt.combined <- IntegrateData(anchorset = wt.anchors)
DefaultAssay(wt.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
wt.combined <- ScaleData(wt.combined, verbose = FALSE)
wt.combined <- RunPCA(wt.combined, npcs = 30, verbose = FALSE)
wt.combined <- RunUMAP(wt.combined, reduction = "pca", dims = 1:30)
wt.combined <- FindNeighbors(wt.combined, reduction = "pca", dims = 1:30)
wt.combined <- FindClusters(wt.combined, resolution = 0.5)

DefaultAssay(wt.combined) <- "RNA"
seucomb <- wt.combined
pop = 'WT'
sccombplt(pop, wt.combined)


# Analyse both Mutants together
mut11seu <- CreateSeuratObject(counts = samcnt[['MUT_11_1']], project = 'MUT_11_1', min.cells = 3, min.features = 200)
mut11seu <- subset(mut11seu, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)
mut15seu <- CreateSeuratObject(counts = samcnt[['MUT_15']], project = 'MUT_15', min.cells = 3, min.features = 200)
mut15seu <- subset(mut15seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
mutseu.list <- list('mut11seu'=mut11seu, 'mut15seu'=mut15seu)

mutseu.list <- lapply(X = mutseu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = mutseu.list, nfeatures = 3000)
mut.anchors <- FindIntegrationAnchors(object.list = mutseu.list, anchor.features = features)
mut.combined <- IntegrateData(anchorset = mut.anchors)
DefaultAssay(mut.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
mut.combined <- ScaleData(mut.combined, verbose = FALSE)
mut.combined <- RunPCA(mut.combined, npcs = 30, verbose = FALSE)
mut.combined <- RunUMAP(mut.combined, reduction = "pca", dims = 1:30)
mut.combined <- FindNeighbors(mut.combined, reduction = "pca", dims = 1:30)
mut.combined <- FindClusters(mut.combined, resolution = 0.5)

DefaultAssay(mut.combined) <- "RNA"
seucomb <- mut.combined
pop = 'MUT'
sccombplt(pop, mut.combined)

# Merge all the samples to get marker genes in all populations
seu.list <- list('wt1seu'=wt1seu, 'wt19seu'=wt19seu, 'mut11seu'=mut11seu, 'mut15seu'=mut15seu)
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features)
seu.combined <- IntegrateData(anchorset = seu.anchors)
DefaultAssay(seu.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.combined <- ScaleData(seu.combined, verbose = FALSE)
seu.combined <- RunPCA(seu.combined, npcs = 30, verbose = FALSE)
seu.combined <- RunUMAP(seu.combined, reduction = "pca", dims = 1:30)
seu.combined <- FindNeighbors(seu.combined, reduction = "pca", dims = 1:30)
seu.combined <- FindClusters(seu.combined, resolution = 0.5)

seucomb <- seu.combined
pop='ALL'
pdf(paste0(CWD, pop, '.pdf'), height=10, width=15)
DimPlot(seu.combined, reduction = "umap")/DimPlot(seu.combined, reduction = "umap", split.by = "orig.ident")
genes <- data.frame()
for(i in sort(unique(seucomb$seurat_clusters))){
  seu.markers <- FindConservedMarkers(seucomb, ident.1 = i, grouping.var = "orig.ident",  logfc.threshold = 1, test.use="MAST", verbose = FALSE)
  varnames <- colnames(seu.markers)[grep("avg_log2FC", colnames(seu.markers))]
  garb <- seu.markers %>%
    filter(max_pval < 0.05) %>%
    filter(!!sym(varnames[1]) > 1, !!sym(varnames[2]) > 1, !!sym(varnames[3]) > 1, !!sym(varnames[4]) > 1) %>%
    #filter(MUT_11_1_avg_log2FC > 1, MUT_11_1_avg_log2FC>1) %>%
    top_n(n = 3, wt = !!sym(varnames[1])+!!sym(varnames[2])+!!sym(varnames[3])+!!sym(varnames[4]))
  head(garb)
  garb['gene'] <- rownames(garb)
  genes <- bind_rows(genes, garb)
  if(nrow(garb) == 0) next
  for(j in 1:nrow(garb)){
    #print(j)
    geneid <- rownames(garb)[j]
    genedf <- garb %>% filter(gene == geneid)
    plt <- FeaturePlot(seucomb, features = geneid, split.by='orig.ident') / VlnPlot(seucomb, features = geneid, split.by='orig.ident')
    parents <- unique(seucomb$orig.ident)
    subt <- ""
    for(par in parents){
      subt <- paste0(subt, par, " - Ave.Log2FC: ", round(genedf[,paste0(par,"_avg_log2FC")], 2),"; Cluster cells %: ", genedf[,paste0(par, "_pct.1")], "; Outside cells %: ", genedf[,paste0(par, "_pct.2")], "\n")
    }
    plt <- plt + plot_annotation(
      title = paste(pop, "Cluster: ", i),
      subtitle = subt
    )
    plot(plt)
  }
}
dev.off()
write.table(genes, paste0(CWD, pop,"_DG.txt"), sep="\t", quote = FALSE, row.names = FALSE)

################################################################################
# Merge gene-lists to find anchor genes

################################################################################
# Test the Seurat::merge pipeline
wt.merge.comb <- merge(wt1seu, y = wt19seu, add.cell.ids = c("wt1", "wt19"), project = "wt")
wt.merge.comb <- NormalizeData(wt.merge.comb)
wt.merge.comb <- FindVariableFeatures(wt.merge.comb, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(wt.merge.comb)
wt.merge.comb <- ScaleData(wt.merge.comb, features = all.genes)
wt.merge.comb <- Seurat::RunPCA(wt.merge.comb, features = VariableFeatures(object = wt.merge.comb))
#wt.merge.comb <- JackStraw(wt.merge.comb, num.replicate = 100, dims=50)
#wt.merge.comb <- ScoreJackStraw(wt.merge.comb, dims = 1:50)
#JackStrawPlot(wt.merge.comb, dims = 1:50)
#ElbowPlot(wt.merge.comb,ndims = 50)
# ndims selected = 30
wt.merge.comb <- FindNeighbors(wt.merge.comb, dims = 1:30)
wt.merge.comb <- FindClusters(wt.merge.comb, resolution = 0.5)
wt.merge.comb <- RunUMAP(wt.merge.comb, dims = 1:30)
wt.merge.comb <- RunTSNE(wt.merge.comb, dims = 1:30)
DimPlot(wt.merge.comb, reduction = "umap", group.by = "orig.ident") / DimPlot(wt.merge.comb, reduction = "umap", split.by = "orig.ident")

mut.merge.comb <- merge(mut11seu, y = mut15seu, add.cell.ids = c("mut_11_1", "mut_15"), project = "mut")
mut.merge.comb <- NormalizeData(mut.merge.comb)
mut.merge.comb <- FindVariableFeatures(mut.merge.comb, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mut.merge.comb)
mut.merge.comb <- ScaleData(mut.merge.comb, features = all.genes)
mut.merge.comb <- Seurat::RunPCA(mut.merge.comb, features = VariableFeatures(object = mut.merge.comb))
mut.merge.comb <- FindNeighbors(mut.merge.comb, dims = 1:30)
mut.merge.comb <- FindClusters(mut.merge.comb, resolution = 0.5)
mut.merge.comb <- RunUMAP(mut.merge.comb, dims = 1:30)
mut.merge.comb <- RunTSNE(mut.merge.comb, dims = 1:30)
DimPlot(mut.merge.comb, reduction = "umap", group.by = "orig.ident") / DimPlot(mut.merge.comb, reduction = "umap", split.by = "orig.ident")

# Merging the datasets directly does not help as there is significant variation
# between the 2 WTs as well as between the 2 MUTs. Therefore, only the integration
# pipeline would be used.


################################################################################
# Check the grouping of the cell-type specific genes
sdata <- readRDS("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/sample_seuobj.rds")
wt1seuobj <- sdata[['wt1seuobj']]
wt19seuobj <- sdata[['wt19seuobj']]
mut11seuobj <- sdata[['mut11seuobj']]
mut15seuobj <- sdata[['mut15seuobj']]
remove(sdata)

gns <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/orthology/cur_ortho_genes.txt',
                  header=FALSE,
                  stringsAsFactors = FALSE)
gns$V1 <- gsub(pattern = 'mRNA', replacement = 'Gene', gns$V1)
gns <- unique(gsub(pattern = '^(Gene\\.\\S+)\\..*', replacement = '\\1', gns$V1))

genes <- data.frame()
pdf(paste0(CWD, "sample_ortho_genes_clusters.pdf"), height=8, width=10)
for(seuobj in c(wt1seuobj, wt19seuobj, mut11seuobj, mut15seuobj)){
  pop <- as.character(unique(seuobj$orig.ident))
  seuobj.markers <- FindAllMarkers(seuobj, test.use="MAST")
  #seuobj.markers %>% filter(gene %in% gns, p_val_adj < 0.05)
  garb <- unique(seuobj.markers %>% filter(gene %in% gns, p_val_adj < 0.05) %>% pull(gene))
  #dp <- DimPlot(seuobj, pt.size=2) + ggtitle(pop)
  #plot(dp)
  #orthoplt <- FeaturePlot(seuobj, features = garb, min.cutoff = "q05", max.cutoff = "q95")
  #orthoplt <- orthoplt + plot_annotation(paste0(pop, ": significantly different homolog genes"))
  #plot(orthoplt)
  garb2 <- seuobj.markers %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.05) %>%
    top_n(n = 1, wt=avg_log2FC)
  garb2['pop'] <- pop
  #markplt <- FeaturePlot(seuobj, features = garb2$gene, min.cutoff = "q05", max.cutoff = "q95")
  #markplt <- markplt + plot_annotation(paste0(pop, ": significantly different marker genes"))
  garb2 <- seuobj.markers %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.05) %>%
    top_n(n = 5, wt=avg_log2FC)
  garb2['pop'] = pop
  genes <- bind_rows(genes, garb2)

  #plot(markplt)
}
dev.off()
genes[, c("p_val", "avg_log2FC", "pct.1", "pct.2","p_val_adj")] <- apply(genes[, c("p_val", "avg_log2FC", "pct.1", "pct.2","p_val_adj")], 2, signif, 3)
write.table(genes, paste0(CWD, "marker_genes.tsv"), sep="\t", quote = FALSE, row.names = FALSE)

highfre <- genes %>% group_by(gene) %>% summarise(cnt=n()) %>% filter(cnt>=3)
write.table(highfre, paste0(CWD, "common_marker_genes.tsv"), sep="\t", quote = FALSE, row.names = FALSE)

geneplot_allsamples <- function(gene, objs, objnames){
  pltlist <- lapply(objnames, function(x) {
    plt1 <- FeaturePlot(objs[[x]], features = gene, min.cutoff = "q05", max.cutoff = "q95")
    plt1 <- plt1 + ggtitle(x) + theme_bw() + theme(text = element_text(size = 8))
  })
  pltout <- ggarrange(plotlist = pltlist, ncol= 2, nrow=2)
  pltout <- annotate_figure(pltout, top = gene, fig.lab.size = 50)
  return(pltout)
}
# Plot common marker genes
objnames <- SAMPLES
objs = list("WT_1"=wt1seuobj,
            "WT_19"=wt19seuobj,
            "MUT_11_1"=mut11seuobj,
            "MUT_15"=mut15seuobj)
plts <- lapply(highfre$gene, geneplot_allsamples, objs, objnames)
pdf(paste0(CWD, "common_marker_genes.pdf"), height=8, width=10)
for(plt in plts){
  plot(plt)
}
dev.off()


## get barcode clusters
outdir<-'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/'
snames <- c('wt1'='WT_1', 'wt19'='WT_19', 'mut11'='MUT_11_1', 'mut15'='MUT_15')
for(s in c('wt1', 'wt19', 'mut11', 'mut15')){
  seu <- sdata[[paste0(s, 'seuobj')]]
  bcs <- gsub(pattern = "\\..*", replacement = "", names(Idents(seu)))
  clstrs <- as.integer(Idents(seu))
  df <- data.frame(Idents(seu))
  df['bcs'] <- gsub(pattern = "\\..*", replacement = "", names(Idents(seu)))
  write.table(df, paste0(outdir, snames[s], "_bcs_clstr_id.txt"), quote = FALSE,row.names = FALSE, col.names = FALSE)
}

## get barcode readcounts
for(s in SAMPLES){
  bcs <- gsub("\\..*", '', colnames(samcnt[[s]]))
  rcs <- colSums(samcnt[[s]])
  df <- data.frame(rcs, bcs)
  write.table(df, paste0(outdir, s, "_bc_readcnts.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
}



head(bcs)
head(clstrs)
head(df)
s<-'wt1'



# Plot for poster
plt1 <- FeaturePlot(wt1seuobj, features = "Gene.31380", min.cutoff = "q05", max.cutoff = "q95")
plt1 <- plt1 + ggtitle("Branch 1") + theme_bw() + theme(text = element_text(size = 20))
plt2 <- FeaturePlot(wt19seuobj, features = "Gene.31380", min.cutoff = "q05", max.cutoff = "q95")
plt2 <- plt2 + ggtitle("Branch 2") + theme_bw() + theme(text = element_text(size = 20))
plt3 <- FeaturePlot(mut11seuobj, features = "Gene.31380", min.cutoff = "q05", max.cutoff = "q95")
plt3 <- plt3 + ggtitle("Branch 3") + theme_bw() + theme(text = element_text(size = 20))
plt4 <- FeaturePlot(mut15seuobj, features = "Gene.31380", min.cutoff = "q05", max.cutoff = "q95")
plt4 <- plt4 + ggtitle("Branch 4") + theme_bw() + theme(text = element_text(size = 20))

plt <- ggarrange(plt1, plt2 ,plt3, plt4, nrow = 2, ncol= 2)
ggsave(plot=plt, "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/plot_poster_gene31380.png", width = 12, height=6)





################################################################################
# DE Analysis
################################################################################
# First run Seurat integration analysis (SCT normalisation) transformation to check
# heterogenity in the single cells. If the cellular populations are not heterogenous
# then run normal bulk RNA DESeq analysis.

# Samcnt was based on manual barcode selection, filtering, and count matrix generation
# This was useful for variant calling, but is not required for normal scRNA analysis.
# Therefore, now using the method described by Seurat for reading the scRNA data
# directly from 10X output
# seuobj.list <- lapply(SAMPLES, FUN=function(x){
#   CreateSeuratObject(counts = samcnt[[x]], project = x, min.cells = 3, min.features = 200)
# })
# Load the PBMC dataset
SAMPLES=c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
SAMPLES=c('WT_1')
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/'
seuobj.list <- sapply(SAMPLES, USE.NAMES=TRUE, FUN=function(x){
  df <- Read10X(data.dir = paste0(CWD, 'get_cells/', x, '/',x,'/outs/filtered_feature_bc_matrix'))
  CreateSeuratObject(counts = df, project = x, min.cells = 3, min.features = 200)
  })
seuobj.list <-sapply(SAMPLES, USE.NAMES=TRUE, FUN=function(x){
  subset(seuobj.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 3000 )
})
plist <- lapply(SAMPLES, FUN=function(x){
  VlnPlot(seuobj.list[[x]], features = c("nFeature_RNA", "nCount_RNA"))
})
names(plist) <- SAMPLES

# Define dataset condition
seuobj.list[['WT_1']]$stim <- 'wt'
seuobj.list[['WT_19']]$stim <- 'wt'
seuobj.list[['MUT_11_1']]$stim <- 'mut'
seuobj.list[['MUT_15']]$stim <- 'mut'
# Transform all datasets
seuobj.list <- lapply(X = seuobj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seuobj.list, nfeatures = 3000)
seuobj.list <- PrepSCTIntegration(object.list = seuobj.list, anchor.features = features)
# Integrate
seuobj.anchors <- FindIntegrationAnchors(object.list = seuobj.list, normalization.method = "SCT", anchor.features = features)
seuobj.combined.sct <- IntegrateData(anchorset = seuobj.anchors, normalization.method = "SCT")
seuobj.combined.sct  <- PrepSCTFindMarkers(seuobj.combined.sct)

DefaultAssay(seuobj.combined.sct) <- "integrated"
seuobj.combined.sct <- RunPCA(seuobj.combined.sct,  npcs = 30, verbose = FALSE)
seuobj.combined.sct <- RunUMAP(seuobj.combined.sct, reduction = "pca", dims = 1:30)
seuobj.combined.sct <- FindNeighbors(seuobj.combined.sct, reduction = "pca", dims = 1:30)
seuobj.combined.sct <- FindClusters(seuobj.combined.sct, resolution = 0.5)
p1 <- DimPlot(seuobj.combined.sct, reduction = "umap", split.by = "orig.ident", label = TRUE)
ggsave(paste0(CWD,'integrate_all_samples.png'), p1, width=12, height=6)

DefaultAssay(seuobj.combined.sct) <- "SCT"
theme_set(theme_cowplot())
clst2.cells <- subset(seuobj.combined.sct, idents = "1")
Idents(clst2.cells) <- "stim"
avg.clst2.cells <- as.data.frame(AverageExpression(clst2.cells, assay='SCT', splot='data', verbose=FALSE)$SCT)
avg.clst2.cells$gene <- rownames(avg.clst2.cells)

p1 <- ggplot(avg.clst2.cells, aes(ctrl, mut)) + geom_point() + ggtitle("Cluster 2 Cells")


# Read genes close to somatic mutations
mut15.genes <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_15.all.layer_somatic_variants.overllaping_100000bp.genes.txt')
mut11.1.genes <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1.all.layer_somatic_variants.overllaping_100000bp.genes.txt')
wt1.genes <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_1.all.layer_somatic_variants.overllaping_100000bp.genes.txt')
wt19.genes <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_19.all.layer_somatic_variants.overllaping_100000bp.genes.txt')
mut.genes <- union(mut15.genes$V1, mut11.1.genes$V1)
mut.genes <- intersect(mut.genes, rownames(seuobj.combined.sct))

# Select all DE genes between the WT and MUT clusters
# Use of scale.data was suggested her https://github.com/satijalab/seurat/issues/5321,
# but other places suggested using slot=data. So using slot=data.
seuobj.combined.sct$celltype.stim <- paste(Idents(seuobj.combined.sct), seuobj.combined.sct$stim, sep = "_")
seuobj.combined.sct$celltype <- Idents(seuobj.combined.sct)
Idents(seuobj.combined.sct) <- "celltype.stim"
for(i in 0:16){
  clst1.response <- FindMarkers(seuobj.combined.sct, assay='SCT', slot='data', ident.1 = paste0(i, "_mut"), ident.2 = paste0(i, "_wt"),verbose = FALSE, features=mut.genes)
  print(i)
  garb <- clst1.response %>%
    filter(abs(avg_log2FC)>0.5) %>%
    filter(p_val_adj < 0.05)
  print(garb)
}
FeaturePlot(seuobj.combined.sct, features = c('Gene.4534', 'Gene.2239', 'Gene.32994', 'Gene.50057'), split.by = "orig.ident", min.cutoff = 1, max.cutoff = 3, cols = c("grey", "red"), pt.size=0.5, order=TRUE)

Idents(seuobj.combined.sct) <- "celltype"
Idents(seuobj.combined.sct) <- paste(Idents(seuobj.combined.sct), seuobj.combined.sct$orig.ident, sep = "_")
markerdf <- data.frame()
for(s1 in c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')){
  for(s2 in c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')){
    if(s1==s2) next()
    print(c(s1, s2))
    for(i in 1:9){
      clst1.response <- FindMarkers(seuobj.combined.sct, assay='SCT', slot='data', ident.1 = paste(i, s1, sep='_'), ident.2 = paste(i, s2, sep='_'), verbose = FALSE)
      if (nrow(clst1.response) == 0){next()}
      garb <- clst1.response %>%
        filter(abs(avg_log2FC)>0.5) %>%
        filter(p_val_adj < 0.05)
      if(nrow(garb) > 0){
        # print(i)
        # print(garb)
        garb['gene'] <- rownames(garb)
        garb['s1'] <- s1
        garb['s2'] <- s2
        garb['i'] <- i
        markerdf <- bind_rows(markerdf, garb)
      }
    }
  }
}


FeaturePlot(seuobj.combined.sct, features = c('Gene.4605', 'Gene.5501'), split.by='orig.ident', min.cutoff = 1, max.cutoff = 3, cols = c("grey", "red"), pt.size=0.5, order=TRUE)


seuobj.combined.sct$type <- paste(Idents(seuobj.combined.sct), 'ctrl', sep = "_")
seuobj.combined.sct$type[seuobj.combined.sct$orig.ident == 'MUT_15'] <- paste(Idents(seuobj.combined.sct)[seuobj.combined.sct$orig.ident == 'MUT_15'], 'mut', sep = "_")
Idents(seuobj.combined.sct) <- "type"
for(i in 0:16){
  clst1.response <- FindMarkers(seuobj.combined.sct, assay='SCT', slot='data', ident.1 = paste0(i, "_mut"), ident.2 = paste0(i, "_ctrl"),verbose = FALSE, features=mut.genes, max.cells.per.ident=500)
  print(i)
  if (nrow(clst1.response) == 0){next()}
  garb <- clst1.response %>%
    filter(abs(avg_log2FC)>0.5) %>%
    filter(p_val_adj < 0.05)
  print(garb)
}


garb <- subset(seuobj.combined.sct, idents = c('0_ctrl', '0_mut'))
VlnPlot(garb, features = 'nFeature_RNA', pt.size = 0.1, combine = TRUE)


clst1.response <- FindMarkers(seuobj.combined.sct, assay='SCT', slot='data', ident.1 = "1_mut", ident.2 = "1_wt",verbose = FALSE, features=mut.genes)
clst1.response %>%
  filter(abs(avg_log2FC)>0.5) %>%
  filter(p_val_adj < 0.05)
# p_val avg_log2FC pct.1 pct.2    p_val_adj
# Gene.2239 8.605091e-103 -0.6330453 0.060 0.388 1.474999e-98
# Gene.4534  6.686506e-86 -0.5982580 0.102 0.425 1.146134e-81

markers.genes <- FindAllMarkers(seuobj.combined.sct, assay='SCT', slot='data', verbose = FALSE, features=mut.genes)
markers.genes['gene'] <- rownames(markers.genes)
markers.genes.filt <- markers.genes %>%
  filter(abs(avg_log2FC)>0.5) %>%
  filter(p_val_adj < 0.05)


VlnPlot(seuobj.combined.sct, features = 'nFeature_RNA', split.by = "stim", pt.size = 0.1, combine = TRUE)
VlnPlot(seuobj.combined.sct, features = 'nFeature_RNA', group.by = "orig.ident", pt.size = 0, idents=c('2_mut', '2_ctrl'), combine = FALSE)



# The above figure shows that MUT_15 is outlier and has some weird expression profile
# when compared to other samples. So now, I use that as the treatment and the rest
# as the control and then find all the genes that have different expression profile in
# MUT_15
seuobj.list <- lapply(SAMPLES, FUN=function(x){
  CreateSeuratObject(counts = samcnt[[x]], project = x, min.cells = 3, min.features = 200)
})
names(seuobj.list) <- SAMPLES
# Define dataset condition
seuobj.list[['WT_1']]$stim <- 'ctrl'
seuobj.list[['WT_19']]$stim <- 'ctrl'
seuobj.list[['MUT_11_1']]$stim <- 'ctrl'
seuobj.list[['MUT_15']]$stim <- 'mut'
# Transform all datasets
seuobj.list <- lapply(X = seuobj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seuobj.list, nfeatures = 3000)
seuobj.list <- PrepSCTIntegration(object.list = seuobj.list, anchor.features = features)
# Integrate
seuobj.anchors <- FindIntegrationAnchors(object.list = seuobj.list, normalization.method = "SCT", anchor.features = features)
seuobj.combined.sct <- IntegrateData(anchorset = seuobj.anchors, normalization.method = "SCT")

################################################################################
# DE Analysis
################################################################################
# Load object created by Anshupa and get BCs and clusters
load('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/clusterObj_umap_res0.5.RData')
orig.ident = int.combinedClst@meta.data$orig.ident
bcs <- names(Idents(int.combinedClst))
clustid <- as.numeric(Idents(int.combinedClst))
wt1.bcs <- gsub('-.*', '', bcs[grep('_1', bcs)])
wt1.clst <- clustid[grep('_1', bcs)]
wt19.bcs <- gsub('-.*', '', bcs[grep('_2', bcs)])
wt19.clst <- clustid[grep('_2', bcs)]
mut11.1.bcs <- gsub('-.*', '', bcs[grep('_3', bcs)])
mut11.1.clst <- clustid[grep('_3', bcs)]
mut15.bcs <- gsub('-.*', '', bcs[grep('_4', bcs)])
mut15.clst <- clustid[grep('_4', bcs)]
write.table(data.frame(wt1.clst, wt1.bcs), '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/wt1_bcs_clstr_id.txt', quote=FALSE, row.names=FALSE)
write.table(data.frame(wt19.clst, wt19.bcs), '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/wt19_bcs_clstr_id.txt', quote=FALSE, row.names=FALSE)
write.table(data.frame(mut11.1.clst, mut11.1.bcs), '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/mut_11_1_bcs_clstr_id.txt', quote=FALSE, row.names=FALSE)
write.table(data.frame(mut15.clst, mut15.bcs), '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/mut_15_bcs_clstr_id.txt', quote=FALSE, row.names=FALSE)

samples <- c('wt1', 'wt19', 'mut_11_1', 'mut_15')
cnts <- int.combinedClst@assays$RNA@counts
for(s in samples){
  ids <- which(orig.ident == s)
  scnts <- cnts[,ids]
  outdf <- data.frame(colSums(scnts))
  outdf['bc'] <- gsub('-.*', '', rownames(outdf))
  write.table(outdf, paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/', s, '_bc_readcnts.txt'), quote=FALSE, row.names=FALSE, col.names=FALSE)
}