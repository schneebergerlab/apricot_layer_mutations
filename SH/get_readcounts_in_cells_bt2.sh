#!/bin/bash

bedbt2=$1
refcur=$2
inpath=$3
bc=$(basename $4)
mkdir $bc
cd $bc
bam-readcount -b 30 -q 10 -w 0 -l $bedbt2 -f $refcur ${inpath}/${bc}/${bc}_dedup.bt2.sorted.bam \
| awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2]}' > ${bc}_readcount.bt2.txt