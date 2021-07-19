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
#bam-readcount -b 30 -q $MAPQ -w 0 -l $bed -f $ref $bam \
#| awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2]}' > $out

bam-readcount -b 30 -q $MAPQ -w 0 -l $bed -f $ref $bam \
| awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > $out

#
#awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); n6=split($11,f, ":"); n7=split($12,g, ":"); n8=split($13,h, ":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}' > $out
#
#
#
#
#bam-readcount -b 30 -q 10 -w 0 -l TMP.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta MUT_11_1/MUT_11_1.sorted.bt2.bam |
#awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}'