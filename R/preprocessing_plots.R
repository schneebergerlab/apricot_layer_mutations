library(ggplot2)
library(dplyr)
library(ggpubr)
library(scales)

bincount <- function(x, nbin){
  x <- sort(x)
  s <- abs(max(x) - min(x))/nbin
  cnts <- sapply(seq(min(x), max(x), s), function(i){
    length(x[(x>=i) & (x<(i+s))])
  })
}
dirs <- c('cur_wtA','cur_wtB','cur_mutC','cur_mutD')
indir <- '/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/'

cols <- c('cadetblue1', 'coral', 'green', 'deeppink')
names(cols) <- dirs
plts <- list()

for(dir in dirs){
  df <- read.table(paste0('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/', dir, '/barcodes_read_count'))
  # df$V2 <- df %>% transmute(V2 = ifelse(V2>10000, 10000, V2)) %>% pull(V2)
  for(i in c(10, 100 , 500, 1000)){
    dffilt <- df %>% filter(V2>=i)
    nbin = 100
    plts[[paste0(dir, i)]] <- ggplot(dffilt)+
      theme_pubr(base_size = 8, base_family = 'Arial') +
      theme(plot.margin = margin(2,2,2,2,'mm'))+
      geom_histogram(aes(x=V2), bins=nbin, alpha=0.1, color=cols[dir], fill=cols[dir], lwd=0.5) +
      # scale_y_log10()+
      # annotate("label", x=max(dffilt$V2)/2,y=10^as.integer(log10(max(table(dffilt$V2)))), label=paste0('Min #reads:', i, '\n#Cell:', nrow(dffilt)), size=8*5/14)+
      annotate("label", x=40000, y=max(bincount(dffilt$V2, nbin)), label=paste0('Min #reads:', i, '\n#Cell:', nrow(dffilt)), size=8*5/14)+
      xlab("Read Count")+
      ylab("Barcode Count")+
      xlim(-1000,60000)+
      ggtitle(dir)
  }
}
plt <- ggarrange(plotlist = plts, nrow = 4, ncol = 4)
ggsave(plot = plt, filename = paste0(indir, 'barcode_read_count_distribution.pdf'), width = 250, height = 250, units = 'mm', device = cairo_pdf)


################################################################################
chrs <- paste0('CUR', 1:8, 'G')

samplebincnt <- lapply(dirs, function(dir){
  barlist <- read.table(paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/',dir,'/barcodes_list'), stringsAsFactors = FALSE)$V1
  barlist <- basename(barlist)
  
  bcbincnt <- lapply(barlist, function(bc){
    print(bc)
    bg <- tryCatch(expr = read.delim(paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/',dir,'/barcodes/',bc,'/',bc,'.bedgraph'), header = FALSE, stringsAsFactors = FALSE) %>% filter(V1 %in% chrs) %>% filter(V4 != 0), error= function(e) 0)
    if(is.double(bg)){
      return(sapply(chrs, function(x){
        0
      }))  
    } else if(is.data.frame(bg)){
      sapply(chrs, function(x){
        df <- bg %>% filter(V1==x)
        length(unique(as.integer(c(df$V2/10000, df$V3/10000))))
      })
    }else{
      print(paste0('ERROR in ', bc))
      print(bg)
    }
  })
  names(bcbincnt) <- barlist
  bcbincnt
})
names(samplebincnt) <- dirs

par(mfrow=c(4,2))
for(dir in dirs){
  bcbincntsum <- data.frame(sapply(samplebincnt[[dir]], sum))
  colnames(bcbincntsum) <- 'bincnt'
  bcbincntsum['bc'] <- rownames(bcbincntsum)
  bcbincntsum['rc'] <- read.table(paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/',dir,'/barcodes_read_count'), stringsAsFactors = FALSE)$V2
  plot(bcbincntsum$rc, bcbincntsum$bincnt, xlab = "read count", ylab = "# 10kb bins", main = dir, pch=20)
  ct <- seq(100, 25000, 100)
  cumbinsnt <- sapply(ct, function(x){
    sum(bcbincntsum$bincnt>x)
  })
  plot(ct, cumbinsnt, xlab = "# 10kb bins", ylab = "# barcodes", main=dir, pch=20)
}
