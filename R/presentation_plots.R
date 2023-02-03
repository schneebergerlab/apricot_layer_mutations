library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(scales)
library(pheatmap)
library(grid)

df <- read.delim2("/local/goel/Dropbox/PostDoc_LMU/presentations/department_seminar_02_2022/dpcr_result_layer_specific.txt",
                 header=TRUE)
colnames(df) <- c('sample', 'LI', 'LII', 'LIII')
df['ID'] <- 1:nrow(df)
df['sample'] <- gsub('_R1', '', df$sample)
df <- melt(df, id.vars = c('sample', 'ID'))
df['value'] <- as.numeric(df$value)
plt <- ggplot(df) +
        theme( panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA)) +
        theme_pubclean() +
        geom_line(aes(x=variable, y = value, group=ID, color=sample), lwd=1) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_x_discrete(expand=c(0,0.1)) +
        ylab("Allele Frequency") +
        xlab("Sample layer")
ggsave(plt, filename = '/local/goel/Dropbox/PostDoc_LMU/presentations/department_seminar_02_2022/layer_specific.png',
       width = 12,
       height = 12,
       units = 'cm',
       dpi=300,
       bg='transparent')


df <- read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/tests/eukaryotes.txt',
                 header=TRUE)
df <- df %>% filter(Status == "Chromosome")
df['year'] <- as.integer(gsub("/.*", "", df$Release.Date))
df <- df %>% filter(year>=2010 & year<2022)
pltdf <- data.frame('year'=names(table(df$year)), "Genomes"=table(df$year))
pltdf <- data.frame(table(df$year))
colnames(pltdf) <- c('year', 'Genomes')
plt <- ggplot(pltdf) +
        theme_pubclean() +
        theme(axis.text.x = element_text(angle=60)) +
        geom_line(aes(year, Genomes, group=1), linetype = "dotted") +
        geom_smooth(aes(year, Genomes, group=1), method = 'loess', se=FALSE) +
        scale_x_discrete(breaks = pretty_breaks(5))
ggsave(plt, filename = '/local/goel/Dropbox/PostDoc_LMU/presentations/department_seminar_02_2022/genome_assembly_count.png',
       width = 8,
       height = 8,
       units = 'cm',
       dpi=1200,
       bg='transparent')

## Draw heatmap for poster @ PAG2023
rcvalues <- read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.read_counts.txt')
# rcvalues <- t(read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.read_counts.txt'))
afvalues <- read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.allele_freq.txt')
annorow <- read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.annocol.txt')
annocol <- read.delim('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.annorow.txt')
rownames(rcvalues) <- paste(rcvalues$X, rcvalues$X.1, rcvalues$X.2, sep='_')
rownames(afvalues) <- paste(afvalues$X, afvalues$X.1, afvalues$X.2, sep='_')
rownames(annorow) <- colnames(rcvalues[4:32])
rownames(annocol) <- rownames(rcvalues)
ann_colors = list(
  var_type = c(Indel='#B9D61A', SNP='#361AD5'),
  branch=c('mut_11_1'='#3182bd','mut_11_2'='#6baed6','mut_15'='#9ecae1','mut_4'='#c6dbef',
           'wt_1'='#e6550d', 'wt_7'='#fd8d3c', 'wt_18'='#fdae6b', 'wt_19'='#fdd0a2'),
  tissue=c('L1'='#006d2c', 'L2'='#31a354', 'L3'='#74c476', 'leaf'='#e6550d')
)
color = colorRampPalette(c("white", "seagreen"))(50)
## These plots were updated on 09.01.2023 after sending the poster for PAG
pheatmap(t(rcvalues[4:32]), cluster_rows = FALSE, cluster_cols=FALSE, annotation_row = annorow[2:3], annotation_col = annocol[3], color=color, annotation_colors=ann_colors, show_rownames = FALSE, show_colnames=FALSE, gaps_row=c(0,0,0,8, 15, 22), gaps_col=c(79), fontsize_row=40, filename='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.read_counts.pdf', width=12, height=4, border_color=NA, main='Read count at somatic mutation positions')
pheatmap(t(afvalues[4:32]), cluster_rows = FALSE, cluster_cols=FALSE, annotation_row = annorow[3], annotation_col = annocol[4], color=color, annotation_colors=ann_colors, show_rownames = FALSE, show_colnames=FALSE, gaps_row=c(0,0,0,8, 15, 22), gaps_col=c(79), fontsize_row=60, annotation_names_row = FALSE, annotation_names_col = FALSE, filename ='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.allele_freq.png', width=8, height=4, border_color=NA, main='Allele frequency at layer-specific somatic mutation positions')


