#!/bin/bash

<<Comment
bedbt2=$1
refcur=$2
inpath=$3
bc=$(basename $4)
mkdir $bc
cd $bc
bam-readcount -b 30 -q 10 -w 0 -l $bedbt2 -f $refcur ${inpath}/${bc}/${bc}_dedup.bt2.sorted.bam \
| awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2]}' > ${bc}_readcount.bt2.txt
Comment


bed=$1
ref=$2
bam=$3
out=$4
MAPQ=$5
bam-readcount -b 30 -q $MAPQ -w 0 -l $bed -f $ref $bam \
| awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2]}' > $out



#
#
#bam-readcount -b 30 -q 10 -w 0 -l ../../strict_syn_snp.selected.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/MUT_11_1/barcodes/AAACCTGAGCCCAACC/AAACCTGAGCCCAACC.DUPmarked.deduped.bam \
#| awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2]}' > AAACCTGAGCCCAACC_read_counts_b30_q10.bt2.txt
#
