# Title     : Iso-Seq analysis
# Objective : Compilation of all R tasks used to analyse Iso-Seq data
# Created by: goel
# Created on: 12/9/21


library(Rsubread)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ggpubr)

################################################################################
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'

gff <- import.gff3('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.protein_coding.3utr.gff3')
gffgenes <- gff[gff@elementMetadata$type == 'gene']
garb <- seqnames(gffgenes)
chrs <- rep(garb@values, garb@lengths)
spos <- start(gffgenes)
epos <- end(gffgenes)
strnd <- strand(gffgenes)
geneid <- mcols(gffgenes)[,"ID"]
ann <- data.frame('GeneID'=geneid, 'Chr'=chrs, 'Start'=spos, 'End'=epos, 'Strand'=strnd)

clslist <- list('WT_1'=seq(8), 'WT_19'=seq(8), 'MUT_11_1'=seq(6), 'MUT_15'=seq(6))
readcounts <- lapply(c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'), function(s){
  rc <- lapply(paste0('clstrs_', clslist[[s]]), function(x){
    featureCounts(files=paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/',s,'/',x,'.iso_seq.bam'), annot.ext=ann)$counts
  })
  df <- data.frame(rc)
  rowSums(df)
})

str(readcounts)
df <- as.data.frame(readcounts)
colnames(df) <- c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
df['rc'] <- rowSums(df)
df <- df %>% filter(rc>10)
df[,1:4][df[,1:4] > 300] <- 300
pltlist <- lapply(c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'), function(s){
  ggplot(df) +
    geom_histogram(aes_string(s), bins=50, alpha=0.5, lwd=0.3, colour='black', fill='grey') +
    theme_pubr() +
    ggtitle(s) +
    xlab("Read count") +
    ylab("Number of genes")
})
ggarrange(plotlist=pltlist, nrow=2, ncol=2)
ggsave(paste0(CWD,"gene_read_cunts.pdf"))
ggsave(paste0(CWD,"gene_read_cunts.png"))