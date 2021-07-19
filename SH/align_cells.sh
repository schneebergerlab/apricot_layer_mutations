#sample=$1
#curidxbt2=$2
#echo $1 | xargs -n 1 -P 1 -I {} /bin/bash -c "
xargs -a barcodes_list -n 1 -P 40 -I {} /bin/bash -c "
    cd {}
    bc=\$(basename {})
    # Align reads using bowtie
#    bowtie2 --end-to-end \
#      --very-sensitive \
#      --threads 1 \
#      -x \$2 \
#      -1 \${bc}_R1.fastq.gz \
#      -2 \${bc}_R2.fastq.gz \
#      --rg-id \$bc \
#    | samtools sort -@1 -O BAM - \
#    > \${bc}.sorted.bam
#    samtools index -@1 \${bc}.sorted.bam
#
#    ## Mark duplicated reads
#    samtools sort -n -O BAM \${bc}.sorted.bam \
#    | samtools fixmate -c -m -O BAM - - \
#    | samtools sort -O BAM - \
#    | samtools markdup -S -s -O BAM - - \
#    | samtools sort -O BAM - > \${bc}.DUPmarked.bam
#    samtools index \${bc}.DUPmarked.bam
#
#    # Get non-dup reads in fastq.gz from the markdup bam files
#    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/remove_duplicates_from_marked_bam_file.py \${bc}.DUPmarked.bam -p \${bc}_dedup
#
#    gzip -f \${bc}_dedup_R1.fastq
#    gzip -f \${bc}_dedup_R2.fastq
#
#    # Get non-dup bam file from the markdup bam file
#    samtools view -G 1024 -O BAM \${bc}.DUPmarked.bam > \${bc}.DUPmarked.deduped.bam
#    samtools index \${bc}.DUPmarked.deduped.bam
#
#    # Get BAM coverage
#    bamCoverage --numberOfProcessors 1 -b \${bc}.DUPmarked.deduped.bam -of bedgraph -o \${bc}.bedgraph
#    hometools bamcov \${bc}.DUPmarked.deduped.bam mean_read_cov_q30_Q10.txt -d \"samtools depth -d 0 -Q 10 -q 30 \"
#    hometools bamcov \${bc}.DUPmarked.deduped.bam mean_read_cov_q30_Q40.txt -d \"samtools depth -d 0 -Q 40 -q 30 \"

    bowtie2 --end-to-end \
      --very-sensitive \
      --threads 1 \
      -x \$3 \
      -1 \${bc}_dedup_R1.fastq.gz \
      -2 \${bc}_dedup_R2.fastq.gz \
      --rg-id \$bc \
    | samtools sort -@1 -O BAM - \
    > \${bc}.DUPmarked.deduped.ORA.bam
    samtools index -@1 \${bc}.DUPmarked.deduped.ORA.bam

    hometools bamcov \${bc}.DUPmarked.deduped.ORA.bam mean_read_cov_q30_Q10.ORA.txt -d \"samtools depth -d 0 -Q 10 -q 30 \"
    hometools bamcov \${bc}.DUPmarked.deduped.ORA.bam mean_read_cov_q30_Q40.ORA.txt -d \"samtools depth -d 0 -Q 40 -q 30 \"


" -- {} $2 $3