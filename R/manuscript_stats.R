# Title     : Collection of different analysis used to generate stats for manuscripts
# Objective : Collection of different analysis used to generate stats for manuscripts
# Created by: goel
# Created on: 12/15/21
library(ggplot2)
library(ggpubr)
library(extrafont)
library(remotes)
library(GenomicRanges)
library(Biostrings)
library(dplyr)
loadfonts()
CURGENSIZE=read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.chrlen.txt',
                      header=FALSE)

# Stats from the distribution of branch-specific mutations
df <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.txt',
                  header = FALSE)
colnames(df) <- c('chr', 'start', 'ref', 'alt', 'altrc', 'altaf', 'sample')
df['type'] = c(rep('SNV', 39), rep('INDEL', 66))
df['highAF'] = df$altaf>0.3
dfgr <- makeGRangesFromDataFrame(df, seqnames.field = 'chr', start.field = 'start', end.field = 'start')
## Plot location of SMs on genome
MAX=max(CURGENSIZE$V2)
plt1 <- ggplot(df) +
  theme_classic(base_family='Arial') +
  theme(axis.text.x = element_text(angle=90)) +
  geom_jitter(aes(x=chr, y=start/MAX, shape=type, color=highAF), size=2) +
  geom_bar(data=CURGENSIZE%>%filter(grepl('CUR', V1)), aes(x=V1, y=V2/MAX), color='blue', alpha=0.1, stat = "identity") +
  xlab("Chromosome") +
  ylab("position") +
  ggtitle("Position of SNPs on chromosome")
ggsave(plot = plt1, filename="/local/goel/Dropbox/projects/apricot_leaf/manuscript/figures/sm_pos.pdf", height=4, width=4, dpi=300)
ggsave(plt1, filename= "/local/goel/Dropbox/projects/apricot_leaf/manuscript/figures/sm_pos.png", height=4, width=4, dpi=300)

## Distribution of SMs in genes
GENECOORDS = read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.bed', header=FALSE)
colnames(GENECOORDS) <- c('chr', 'start', 'end', 'id')
GENECOORDS = makeGRangesFromDataFrame(GENECOORDS)
EXONCOORDS = read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.exons.bed', header=FALSE)
colnames(EXONCOORDS) <- c('chr', 'start', 'end', 'id')
EXONCOORDS <- makeGRangesFromDataFrame(EXONCOORDS)

### Number of SMs in Genes and Exons
df['ingenes'] <- FALSE
df['inexon'] <- FALSE
df[unique(queryHits(findOverlaps(dfgr, GENECOORDS))), 'ingenes'] <- TRUE
df[unique(queryHits(findOverlaps(dfgr, EXONCOORDS))), 'inexon'] <- TRUE

## Get whether the SNVs are transitions or travsversions
df['subs'] <- 'NA'
change <- paste0(df$ref, df$alt)
df[!grepl('-|\\+', change), 'subs'] <- 'Tv'
df[change %in% c('AG', 'CT', 'GA', 'TC'), 'subs'] <- 'Ts'
table(df$subs)

## Check whether the G/C to A/T transitions are in dyrimidine region
seq <- readDNAStringSet("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/snps_with_neighbour_sequence.fa")
seq = data.frame(seq)
seq['pos'] <- rownames(seq)
rownames(seq) <- NULL

posvec <- paste(df$chr, df$start-2, df$start+1, sep='_')
df['seq'] <- NA
df[1:39, 'seq'] <- seq$seq[match(posvec, seq$pos)[1:39]]
### Check if the SNV is in dipyrimidine region
isdip <- apply(df, MARGIN = 1, function(row){
  if(row['subs'] != 'Ts'){
    return(NA)
  }
  else if (! row['ref'] %in% c('C', 'G')){
  return(NA)
  }
  else if (row['ref'] == 'C'){
    print(paste('C', row['seq'], row['seq'][3]))
    if (substring(row['seq'], 1, 1) %in% c('C', 'T')) return(TRUE) else return(FALSE)
  }
  else if (row['ref'] == 'G'){
    print(paste('G', row['seq'], row['seq'][3]))
    if (substring(row['seq'], 3, 3) %in% c('G', 'A')) return(TRUE) else return(FALSE)
  }
})
df['isdipy'] <- isdip

## Which somatic mutations are in syntenic/SR regions
synr = read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syn_range.txt', header=FALSE)
colnames(synr) <- c('chr', 'start', 'end')
synr = makeGRangesFromDataFrame(synr)
srr = read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/sr_range.txt', header=FALSE)
colnames(srr) <- c('chr', 'start', 'end')
srr = makeGRangesFromDataFrame(srr)
df['struct'] <- NA
df[unique(queryHits(findOverlaps(dfgr, synr))), 'struct'] <- 'SYN'
df[unique(queryHits(findOverlaps(dfgr, srr))), 'struct'] <- 'SR'

write.table(df, file = '/local/goel/Dropbox/projects/apricot_leaf/manuscript/mutations.stats.csv', sep='\t', quote=FALSE, row.names=FALSE)


length(unique(queryHits(findOverlaps(dfgr, GENECOORDS))))
length(unique(queryHits(findOverlaps(dfgr, EXONCOORDS))))
table(df[unique(queryHits(findOverlaps(dfgr, GENECOORDS))), 'V8'])
table(df[unique(queryHits(findOverlaps(dfgr, EXONCOORDS))), 'V8'])
head(df)
