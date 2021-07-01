#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()
parser$add_argument("indir", help='path to directory containing read counts at SNP marker positions')
# parser$add_argument("pattern", help='Pattern to select files')
parser$add_argument('chr', help="Chromosome being analysed", type='character')
parser$add_argument("outdir", help='path to output directory')
parser$add_argument("-R", help='Rigidity value', default=500, type="integer")


args <- parser$parse_args()

## Set path to folder containing the input file
file_paths = list.files(path = args$indir, full.names = TRUE)
# file_paths = list.files(path = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/CUR1G/input/', full.names = TRUE)[1:3]
sampleIDs <- basename(file_paths)

# Create the expDesign object
expDesign = data.frame(files=file_paths, name=sampleIDs)

chr_len = c('CUR1G'= 46975282,
            'CUR2G'= 33806098,
            'CUR3G'= 26861604,
            'CUR4G'= 26096899,
            'CUR5G'= 18585576,
            'CUR6G'= 27136638,
            'CUR7G'= 25539660,
            'CUR8G'= 23045982)

library(Gviz)
library(RTIGER)
sourceJulia()
options(ucscChromosomeNames=FALSE)

# print(args$R)

mydat <- RTIGER::RTIGER(expDesign = expDesign,
               rigidity = args$R,
               outputdir = args$outdir,
               seqlengths = chr_len[args$chr],
               save.results = TRUE)

saveRDS(mydat, paste0(args$outdir, '/RTIGER_R', args$R, '.RDS'))

# Try using this RTIGER call. Would need to set seqlengths to match the chromosome that you are analysing.
# myDat2 <- RTIGER::RTIGER(expDesign = expDesign,
#        rigidity = 500,
#        outputdir = 'rtiger_co',
#        seqlengths = chr_len['CUR1G'],
#        save.results = TRUE)
