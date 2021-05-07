library(dplyr)
library(ggplot2)
library(pheatmap)
library(igraph)
df <- read.table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mut_rc.txt',
                 header = TRUE,
                 stringsAsFactors = FALSE)
df <- df[grepl('CUR', df$chr),]
df$end = as.integer(df$end)

samples = c('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
sample_df <- list()
for(sample in samples){
  df2 = read.table(paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/', sample, '/', sample,'_good_candidates.txt'),
                   header = FALSE,
                   stringsAsFactors = FALSE)
  df2[,'V2'] <- NULL
  colnames(df2) <- c('chr', 'end', 'ref', 'alt', 'depth', 'freq')
  sdf = inner_join(df, df2, by=c('chr', 'end', 'ref', 'alt', 'depth', 'freq'))
  print(dim(sdf))
  sdf = sdf[, c('chr', 'start', 'end', 'ref', 'alt', 'depth', 'freq', colnames(sdf)[grep(paste0(sample, '\\.'), colnames(sdf))])]
  print(dim(sdf))
  sample_df[[sample]] <- sdf
}

par(mfrow=c(2,2))
for(sample in samples){
  cnts = colSums(sample_df[[sample]][,seq(9, ncol(sample_df[[sample]]), 2)])
  empty = sum(cnts==0)
  hist(cnts,
       main = paste(sample, 'mutations per cell\nNumber of cells with 0 muts: ', empty),
       xlab = 'Number of Mutations',
       ylab='Number of cells')
}

for(sample in samples){
  cnts = colSums(sample_df[[sample]][sample_df[[sample]]$depth>=10,seq(9, ncol(sample_df[[sample]]), 2)])
  cnts = cnts[which(cnts!=0)]
  high = sum(cnts>=2)
  hist(cnts,
       main = paste(sample, 'High_conf mutations per cell\nNumber of cells >=2 muts: ', high),
       xlab = 'Number of Mutations',
       ylab ='Number of cells')
}

# par(mfrow=c(2,2))
for(sample in samples){
  s1df <- sample_df[[sample]][sample_df[[sample]]$depth>=10,seq(9, ncol(sample_df[[sample]]), 2)]
  s1df <- t(s1df)
  s1df <- s1df[rowSums(s1df) > 0, ]
  pheatmap(s1df, show_colnames=FALSE,
           filename=paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/',
                           sample,'/', sample,
                           '_mut_heatmap.pdf'),
           main=paste0(sample, '\nNumber of mutations: ', ncol(s1df), '\n Number of cells: ', nrow(s1df)))
  s1.pca <- prcomp(s1df)
  p1var <- as.numeric(format((s1.pca$sdev^2/sum(s1.pca$sdev^2))[1], digits=3))*100
  p2var <- as.numeric(format((s1.pca$sdev^2/sum(s1.pca$sdev^2))[2], digits=3))*100
  p3var <- as.numeric(format((s1.pca$sdev^2/sum(s1.pca$sdev^2))[3], digits=3))*100
  
  s1.km <- kmeans(s1.pca$x, centers=3, iter.max = 10000)
  pdf(paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/',
             sample,'/', sample,
             '_mut_pca_clust.pdf'))
  plot(s1.pca$x[,1], s1.pca$x[,2],
       col=s1.km$cluster, 
       pch=16,
       cex=1,
       alpha=0.2,
       xlab=paste0("PC1 (", p1var,"%)"),
       ylab=paste0("PC2 (", p2var,"%)"),
       main = paste0(sample, " PC1 / PC2 - plot\nClusters generated using all PCs"))
       # main = paste0(sample, " PC1 / PC2 - plot\nCluster using first three PC\nVariance considered: ", p1var+p2var+p3var))
  dev.off()
}

pdf(paste0('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/',
           'mut_cooccurence.pdf'), width = 10, height = 10)
par(mfrow= c(2,2))
par(mar=c(0,0,3,0))
par(mgp = c(0,0,0))
for(sample in samples){
  
  s1df <- sample_df[[sample]][sample_df[[sample]]$depth>=10,]
  bad_pos <- c()
  for(i in 1:nrow(s1df)){
    for(j in 1:nrow(s1df)){
      if(i==j) next
      if(s1df[i,'chr'] == s1df[j,'chr']){
        if(abs(s1df[i,'end'] -  s1df[j,'end']) <= 1000){
          bad_pos <- c(bad_pos, i)
          break
        }
      }
    }
  }
  if(length(bad_pos) >0 ){
    s1df <- s1df[-bad_pos,]
  }  
  meta_df <- s1df[, 1:7]
  s1df <- s1df[, seq(9, ncol(sample_df[[sample]]), 2)]
  s1df <- s1df[, colSums(s1df) > 0]
  
  mutFreq <- sapply(1:nrow(s1df), function(i){
    sum(s1df[i,] > 0)
  })
  edge <- c()
  weight <- c()
  for(i in 1:(nrow(s1df)-1)){
    for(j in (i+1):nrow(s1df)){
      # print(c(i,j))
      s <- sum(s1df[i,] & s1df[j,])
      if(s>0){
        edge <- c(edge, i, j)
        weight <- c(weight, s)
      }
    }
  }
  g1 <- graph(edge, directed=FALSE, n = nrow(s1df))
  E(g1)$weight <- weight
  
  for(i in c(0, 3, 5, 10)){
    plt_g1 <- g1
    plt_g1 <- delete.edges(plt_g1, which(E(plt_g1)$weight < i))
    coords <- layout_with_kk(plt_g1)
    plot(plt_g1, edge.width=E(plt_g1)$weight,
         edge.color=E(plt_g1)$weight,
         vertex.size=(mutFreq)^0.666,
         layout=coords,
         main=paste0(sample, "\nEdges with weight < ", i, " removed"))
  }
  
  for(i in c(0.2, 0.3, 0.4, 0.5)){
    high_af <- which(meta_df$freq > i)
    plt_g1 <- g1
    plt_g1 <- delete.vertices(plt_g1, high_af)
    coords <- layout_with_kk(plt_g1)
    plot(plt_g1, edge.width=E(plt_g1)$weight,
         edge.color=E(plt_g1)$weight,
         vertex.size=(mutFreq[-high_af])^0.666,
         layout=coords,
         main=paste0(sample, "\nMutations with AF > ", i, " removed"))
  }
}
dev.off()

