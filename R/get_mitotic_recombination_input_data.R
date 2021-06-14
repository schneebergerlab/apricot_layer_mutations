# Using RTIGER get mitotic recombination in single-cells
library(RTIGER)
library(dplyr)
library(Gviz)
library(ggplot2)
library(ggpubr)
library(data.table)
library(Sushi)

# setupJulia()
sourceJulia()

CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/'
setwd(CWD)
SAMPLES=c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'


fins = list.files(CWD, pattern = '_b30_q10.bt2.txt', recursive = TRUE, full.names = TRUE)
fins.selected = grep('read_count', fins, invert = TRUE, value=TRUE)
samples = basename(fins.selected)
sids = c()
for(f in samples){
  s = strsplit(f, '\\.')[[1]][1]
  print(s)
  sids = c(sids, s)
}

# Get sequenced mutation count per sample

snp_cnt <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp.selected.txt') %>% group_by(V1) %>% summarise(cnt=n())
good_sids = c()
for(sample in SAMPLES){
  rc = fread(paste0(INDIR, sample, '/barcodes_read_count_deduped')) %>%
    filter(V2>=200000)
  good_sids = c(good_sids, paste0(sample, '_', rc$V1,'_b30_q10'))
}

isgood <- sids %in% good_sids
fins.selected <- fins.selected[isgood]
sids <- sids[isgood]

celldf <- data.frame()
for(f in fins.selected){
  df <- fread(f) %>% group_by(V1) %>% summarise(cnt=n())
  df$cnt <- df$cnt/snp_cnt$cnt
  df$sample <- strsplit(basename(f), '\\.')[[1]][1]
  m <- mean(df$cnt)
  df$cnt <- df$cnt - m
  df$`x` <- 1:8
  celldf <- bind_rows(celldf, df)
}
celldf$off <- "Normal"
celldf$off[celldf$cnt > 3*sd(celldf$cnt)] <- "High"
celldf$off[celldf$cnt < -3*sd(celldf$cnt)] <- "Low"

# abs(celldf$cnt) > 3*sd(celldf$cnt)
offsample <- celldf %>% filter(off %in% c('High', 'Low')) %>% pull(sample)
offcelldf <- celldf %>% filter(sample %in% offsample)
ggplot(offcelldf) +
  geom_line(aes(x=x, y=cnt, group=sample))

chr_len = c('CUR1G'= 46975282,
            'CUR2G'= 33806098,
            'CUR3G'= 26861604,
            'CUR4G'= 26096899,
            'CUR5G'= 18585576,
            'CUR6G'= 27136638,
            'CUR7G'= 25539660,
            'CUR8G'= 23045982)



pdf(paste0(CWD, 'aneuploidy_candidates.pdf'), width=9, height=6)
for(s in unique(offsample)){
  pop <- substr(s, 1, nchar(s)-25)
  bc <- substr(s,nchar(s)-23,nchar(s)-8)
  sample_snpcnt <- offcelldf %>% filter(sample == s)
  col <- sample_snpcnt$off
  names(col) <- sample_snpcnt$V1
  
  df <- fread(paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/',pop,'/barcodes/',bc,'/',bc,'.bedgraph'))  %>% filter(!grepl('pilon', V1))
  df_small <- sample_n(df, 0.2*nrow(df))
  # df$V4[df$V4>10] <- 10
  df_small$V4[df_small$V4>20] <- 20
  df_small$col = col[df_small$V1]
  df_small <- df_small[order(df_small$V2),]
  df_small <- df_small[order(df_small$V1),]
  
  plt <- ggplot(df_small) +
    theme_classic() +
    facet_grid(V1~.) +
    theme(text = element_text(size=8),
          legend.position = c(0.8, 0.4))+ 
    geom_line(aes(V2, V4, color=df_small$col)) +
    scale_color_manual(values = c("High"='red', 'Normal'='grey', 'Low'='blue'), name='Syntenic SNP\nfrequency') +
    ylab("Coverage") +
    xlab('Position') +
    ggtitle(paste(pop, bc))
  
  df0 <- df %>% filter(V4==0) %>% group_by(V1)  %>% mutate(len=V3-V2+1) %>% summarise(total_len = sum(len))
  df0$total_len <- 1 - df0$total_len/chr_len[df0$V1]
  plt1 <- ggplot(df0) +
    theme_bw() +
    theme(text = element_text(size=8),
          axis.title.y = element_blank(),
          axis.text.y = element_text(angle = 90)) +
    geom_bar(aes(x=V1, y=total_len), stat="identity", alpha=0.5)+
    scale_x_discrete(limits=rev) +
    ylim(0,1)+
    ylab("% sequenced") +
    coord_flip()
    
  plot(ggarrange(plt, plt1,ncol = 2,widths = c(0.7, 0.3)))
  
  # c1 <- df %>% filter(V1=='CUR1G')
  # 
  # c1$V4[c1$V4>10] <- 10
  # x <- as.vector(rbind(df$V2, df$V3))
  # y <- rep(df$V4, rep(2, nrow(df)))
  # plot(x, y, type='l',ylim=c(0,10))
  # plot(c1$V2, c1$V4, type='l', ylim=c(0,10))
}
dev.off()

################################################################################
# Run RTIGER on individual chromosomes
# For each cell-chromosome pair, check if enough SNP markers are sequenced.
# Remove chromosomes, with few sequenced markers.
################################################################################

chr_len = c('CUR1G'= 46975282,
            'CUR2G'= 33806098,
            'CUR3G'= 26861604,
            'CUR4G'= 26096899,
            'CUR5G'= 18585576,
            'CUR6G'= 27136638,
            'CUR7G'= 25539660,
            'CUR8G'= 23045982)
chrs <- names(chr_len)
setwd('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/')

## For Mapping Quality: 10
fins = unlist(sapply(SAMPLES, function(sample){
  list.files(paste0(CWD, sample), pattern = '_b30_q10.bt2.txt', recursive = TRUE, full.names = TRUE)
}, USE.NAMES = FALSE))
fins = grep('read_count', fins, invert = TRUE, value=TRUE)
for(chr in chrs) dir.create(paste0(chr, '/input'), recursive = TRUE)
a <- lapply(fins, function(f){
  df = fread(f)
  if(nrow(df)==0) return()
  df2 = df %>% group_by(V1) %>% summarize(cnt=n())
  for(chr in chrs){
    if(chr %in% df$V1){
      if(df2[df2$V1==chr,]$cnt >= 5000){
        write.table(df %>% filter(V1==chr), file = paste0(chr, '/input/', chr, '_', basename(f)), quote = FALSE, col.names = FALSE, row.names=FALSE, sep = '\t')
      }  
    }
  }
})


## For Mapping Quality: 40
library(doParallel)
library(foreach)
cl <- makeCluster(40)
registerDoParallel(cl)

fins = unlist(sapply(SAMPLES, function(sample){
  list.files(paste0(CWD, sample), pattern = '_b30_q40.bt2.txt', recursive = TRUE, full.names = TRUE)
}, USE.NAMES = FALSE))
fins = grep('read_count', fins, invert = TRUE, value=TRUE)
for(chr in chrs) dir.create(paste0(chr, '/input_q40'), recursive = TRUE)

foreach(f=fins, .packages = c('data.table', 'dplyr')) %dopar% {
  df = fread(f)
  if(nrow(df)==0) return()
  df2 = df %>% group_by(V1) %>% summarize(cnt=n())
  for(chr in chrs){
    if(chr %in% df$V1){
      if(df2[df2$V1==chr,]$cnt >= 5000){
        write.table(df %>% filter(V1==chr), file = paste0(chr, '/input_q40/', chr, '_', basename(f)), quote = FALSE, col.names = FALSE, row.names=FALSE, sep = '\t')
      }  
    }
  }
}
stopCluster(cl)

a <- lapply(fins, function(f){
  df = fread(f)
  if(nrow(df)==0) return()
  df2 = df %>% group_by(V1) %>% summarize(cnt=n())
  for(chr in chrs){
    if(chr %in% df$V1){
      if(df2[df2$V1==chr,]$cnt >= 5000){
        write.table(df %>% filter(V1==chr), file = paste0(chr, '/input_q40/', chr, '_', basename(f)), quote = FALSE, col.names = FALSE, row.names=FALSE, sep = '\t')
      }  
    }
  }
})


# Create Inverted data set for Chromosome CUR6G to test RTIGER
setwd(paste0(CWD, 'rtiger_out/CUR6G/'))
dir.create('reversed_input_q40', showWarnings = FALSE)
fins <- list.files('input_q40')
lapply(fins, function(f){
  df <- read.table(paste0('input_q40/', f))
  df$V2 <- chr_len['CUR6G'] - df$V2
  df <- df[order(df$V2),]
  write.table(df, paste0('reversed_input_q40/', f),quote = FALSE,sep = '\t',row.names = FALSE, col.names = FALSE)
})
