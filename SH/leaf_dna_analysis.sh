Analyse the leaf dna from the four branches sequenced in the second run
NOTE: FOR ALL ANALYSIS CURROT assembly would be used as the primary assembly

################################################################################
############# Pre-run: Get repeat regions in the currot assembly ###############
################################################################################
### Not required as done in the scdna_analysis.sh

################################################################################
############################ Pre-process reads #################################
################################################################################
hometools='/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py'

### Get sequencing quality plots
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot_branch_specific/
INDIR=/srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_5350/X204SC21122726-Z01-F001/raw_data/
declare -A sample_libid
sample_libid['A']='wt7'
sample_libid['B']='wt18'
sample_libid['C']='mut4'
sample_libid['D']='mut11_2'

cd $CWD
for lib in A B C D; do
    {
        cd $CWD
        mkdir ${sample_libid[$lib]}; cd ${sample_libid[$lib]}
#        ln -s ${INDIR}/*_${lib}/*fq.gz .
        fastqc -t 10 *fq.gz
    } &
done

### Merge, trim and align reads to CUR genomes
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot_branch_specific/
cd $CWD
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
for s in wt7 wt18 mut4 mut11_2; do
        cd $CWD
        bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo ${s}.log -eo ${s}.err "
            cd $s
            cat ${INDIR}/${s}/A5350*_1.fq.gz > merged_read_1.fq.gz
            cat A5350*_2.fq.gz > merged_read_2.fq.gz
            skewer -r 0.1 -d 0.05 -k 8 -q 20 -l 75 -m pe -t 15 -x \
            /srv/netscratch/dep_mercier/grp_schneeberger/projects/hyper_co/data/reads/shqMergedAdapters_Primers_representative_rc.fa \
            merged_read_1.fq.gz \
            merged_read_2.fq.gz \
            -z -o ${s}_ql

            bowtie2 --end-to-end \
                --very-sensitive \
                --threads 20 \
                -x $curidxbt2 \
                -1 ${s}_ql-trimmed-pair1.fastq.gz \
                -2 ${s}_ql-trimmed-pair2.fastq.gz \
                --rg-id ${s}\'\\t\'PL:ILLUMINA\'\\t\'SM:1 \
            | samtools view -b - \
            > ${s}.bam

            ## Mark duplicated reads
            samtools sort -@ 10 -n -O BAM ${s}.bam \
            | samtools fixmate -@ 10 -c -m -O BAM - - \
            | samtools sort -@ 10 -O BAM - \
            | samtools markdup -@ 10 -S -s -O BAM - - \
            | samtools sort -@ 10 -O BAM - > ${s}.DUPmarked.bam
            samtools index ${s}.DUPmarked.bam

            samtools view -G 1024 -O BAM ${s}.DUPmarked.bam > ${s}.deduped.bam
            samtools index ${s}.deduped.bam
        "
done

### Align reads to telomere sequence
telfa=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/telomeric_sequence.fa
for s in wt7 wt18 mut4 mut11_2; do
        cd $CWD; cd $s
        bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo ${s}_tel.log -eo ${s}_tel.err "
            bowtie2 --end-to-end \
                --very-sensitive \
                --threads 20 \
                -x $telfa \
                -1 ${s}_ql-trimmed-pair1.fastq.gz \
                -2 ${s}_ql-trimmed-pair2.fastq.gz \
                --rg-id ${s}\'\\t\'PL:ILLUMINA\'\\t\'SM:1 \
            | samtools view -h -@20 -F4 - \
            | samtools sort -@20 -O BAM - \
            > ${s}.telo.bam
        "
done


### Get average of read depth in the genome and in the telomere
for s in wt7 wt18 mut4 mut11_2; do
    cd ${CWD}/${s}
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=15000]" -M 20000 -oo ${s}_cov.log -eo ${s}_cov.err "
       $hometools bamcov -g ${s}.deduped.bam ${s}_vs_cur.deduped.cov
    "
done
# Telomere cov ratio:
# wt7: 54.2; wt18: 58.8; mut4: 67; mut11_2: 54.35334785924399

####################################################################
############ Get read-counts at variant positions
####################################################################
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
CHRBED=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.chrlen.bed
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
for s in wt7 wt18 mut4 mut11_2; do
    cd ${CWD}/$s
    bsub -q ioheavy -n 8  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc.log -eo ${s}_bamrc.err -m 'hpc001 hpc002 hpc003 hpc005 hpc004' "
        $hometools pbamrc -n 8 -b 30 -q 10 -w 0 -S -I -f $refcur -l $CHRBED ${s}.deduped.bam bam_read_counts_b30_q10.bt2.txt

      # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
      /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q10.bt2.txt
  "
done

for s in wt7 wt18 mut4 mut11_2; do
    cd ${CWD}/$s
    bsub -q ioheavy -n 8  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc_q0.log -eo ${s}_bamrc_q0.err -m 'hpc001 hpc002 hpc003 hpc005 hpc004' "
        $hometools pbamrc -n 8 -b 0 -q 0 -w 0 -S -I -f $refcur -l $CHRBED ${s}.deduped.bam bam_read_counts_b0_q0.bt2.txt
  "
done

## Get read mapping depth histogram
for s in wt7 wt18 mut4 mut11_2; do
    cd ${CWD}/$s
#    cut -f 4 bam_read_counts_b30_q10.bt2.txt | sort -n | uniq -c | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.pdf -x Mapping Depth -y Frequency -t ${s}_bt2 -xlim 0 400 &
    cut -f 4 bam_read_counts_b30_q10.bt2.txt | sort -n | uniq -c | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.png -x Mapping Depth -y Frequency -t ${s} -xlim 0 400 -n 400 &
done

### Variant identification
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/leaf_dna_analysis.py


####################################################################
############ Step 4: Gene conversions
####################################################################

# <editor-fold desc="Initial attempts for gene conversion identification that did not work">
### As initial check, get the alt allele-frequency distribution at the
### heterozygous positions identified by syri
#syriout=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out
#syrisnp=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt
#grep 'SNP' $syriout | cut -f1,2,3,4,5 | grep -v 'N' > $syrisnp
#CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
#for sample in wt7 wt18 mut4 mut11_2 ; do
#    cd ${CWD}/${sample}
#    bsub -q multicore40 -n 40  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc_q0.log -eo ${s}_bamrc_q0.err -m 'hpc001 hpc002 hpc003 hpc005 hpc004' "
#        $hometools pbamrc -n 8 -b 30 -q 10 -w 0 -I -n 40 -f $refcur -l $syrisnp ${sample}.deduped.bam ${sample}.syri_snps.bamrc
#    "
#done
#
#for sample in MUT_11_1 MUT_15 WT_1 WT_19 ; do
#    cd ${CWD}/${sample}
#    bsub -q multicore40 -n 40  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc_q0.log -eo ${s}_bamrc_q0.err -m 'hpc001 hpc002 hpc003 hpc005 hpc004' "
#        $hometools pbamrc -n 8 -b 30 -q 10 -w 0 -I -n 40 -f $refcur -l $syrisnp ${sample}.sorted.bt2.bam ${sample}.syri_snps.bamrc
#    "
#done
#SNPS=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/snps_close.bed
## Get reads overlapping SNPs
#for sample in wt7 wt18 mut4 mut11_2 ; do
#    cd ${CWD}/${sample}
#    N=20
#    bsub -q multicore40 -n 20  -R "span[hosts=1] rusage[mem=40000]" -M 50000 -oo ${sample}_snps_reads.log -eo ${sample}_snps_reads.err "
#        samtools view -b -L $SNPS -f 2 -@ $N ${sample}.deduped.bam \
#        | samtools sort -n -O SAM -@ $N - \
#        | grep -vP '^@' \
#        | cut -f 1 \
#        | uniq >  snps_read_name.txt
#
#        samtools view -H ${sample}.deduped.bam > snps_reads.sam
#        samtools sort -O SAM -@ $N -n ${sample}.deduped.bam \
#        | grep -Ff snps_read_name.txt >> snps_reads.sam
#        samtools sort -@ $N -O BAM snps_reads.sam > snps_reads.bam
#        samtools index -@ $N snps_reads.bam
#        rm snps_reads.sam
#    "
#done
#for sample in MUT_11_1 MUT_15 WT_1 WT_19 ; do
#    cd ${CWD}/${sample}
#    N=20
#    bsub -q multicore40 -n 20 -R "span[hosts=1] rusage[mem=40000]" -M 50000 -oo ${sample}_snps_reads.log -eo ${sample}_snps_reads.err "
#        samtools view -b -L $SNPS -f 2 -@ $N ${sample}.sorted.bt2.bam \
#        | samtools sort -n -O SAM -@ $N - \
#        | grep -vP '^@' \
#        | cut -f 1 \
#        | uniq >  snps_read_name.txt
#
#        samtools view -H ${sample}.sorted.bt2.bam > snps_reads.sam
#        samtools sort -O SAM -@ $N -n ${sample}.sorted.bt2.bam \
#        | grep -Ff snps_read_name.txt >> snps_reads.sam
#        samtools sort -@ $N -O BAM snps_reads.sam > snps_reads.bam
#        samtools index -@ $N snps_reads.bam
#        rm snps_reads.sam
#    "
#done
#
## Get pileup data at SNP positions
#awk '{print $1":"$3"-"$3}' /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/snps_close.bed > snps_close.txt
#mkdir snps_split
#for i in {1..8}; do
#    grep CUR${i}G snps_close.txt > snps_split/snps_CUR${i}G.txt
#done
#
#for sample in wt7 wt18 mut4 mut11_2 MUT_11_1 MUT_15 WT_1 WT_19; do
#    cd ${CWD}/${sample}
#    N=20
#    for i in {1..8}; do
#        snps="../snps_split/snps_CUR${i}G.txt"
#        bsub -q multicore40 -n $N -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${sample}_snps_pileup.log -eo ${sample}_snps_pileup.err "
#            xargs -a $snps -P $N -I {} samtools mpileup -f $refcur -q 40 -E -Q 26 snps_reads.bam -r {} -O --output-QNAME > snps_CUR${i}G.pileup 2> garb
#            sort -k2,2n -o snps_CUR${i}G.pileup snps_CUR${i}G.pileup
#        "
#    done
#done
#
## Get pileup data at ShV positions (added when analysing layer-specific gene conversions)
#cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling
#awk '{if($4~/^-/){print $1"\t"$2-1"\t"$2"\t\t"$4} else {print $0}}' shv_close.bed > shv_close.for_mpileup.bed
#awk '{print $1":"$3"-"$3}' shv_close.bed > shv_close.txt
#awk '{if($4~/^-/){print $1":"$2"-"$2} else {print $1":"$3"-"$3}}' shv_close.bed > shv_close.txt
#mkdir -p shv_split
#for i in {1..8}; do
#    grep CUR${i}G shv_close.txt > shv_split/shv_CUR${i}G.txt
#done
#
### TODO: RUN THE BELOW PILEUP COMMAND FOR LEAF DATA
##for sample in wt7 wt18 mut4 mut11_2 MUT_11_1 MUT_15 WT_1 WT_19; do
##    cd ${CWD}/${sample}
##    N=20
##    for i in {1..8}; do
##        snps="../snps_split/snps_CUR${i}G.txt"
##        bsub -q multicore40 -n $N -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${sample}_snps_pileup.log -eo ${sample}_snps_pileup.err "
##            xargs -a $snps -P $N -I {} samtools mpileup -f $refcur -q 40 -E -Q 26 snps_reads.bam -r {} -O --output-QNAME > snps_CUR${i}G.pileup 2> garb
##            sort -k2,2n -o snps_CUR${i}G.pileup snps_CUR${i}G.pileup
##        "
##    done
##done
#
## Get pileup data at selected gene conversion positions
#for sample in wt7 wt18 mut4 mut11_2 MUT_11_1 MUT_15 WT_1 WT_19; do
#    cd ${CWD}/${sample}
#    snps="../gene_conversions_positions.txt"
#    {
##    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${sample}_geneconv_pileup.log -eo ${sample}_geneconv_pileup.err "
#        xargs -a $snps -P 10 -I {} samtools mpileup -f $refcur -q 40 -E -Q 26 snps_reads.bam -r {} -O --output-QNAME > geneconv.pileup 2> garb
#        sort -k1,1 -k2,2n -o geneconv.pileup geneconv.pileup
##    "
#    } &
#done
#
## Get bam corresponding to geneconv.reads.txt and pileup
#for sample in wt7 wt18 mut4 mut11_2 MUT_11_1 MUT_15 WT_1 WT_19; do
#    cd ${CWD}/${sample}
#    {
#        samtools view -H snps_reads.bam > geneconv.reads.sam
#        samtools view snps_reads.bam | grep -Ff geneconv.reads.txt >> geneconv.reads.sam
#        samtools view -b -T $refcur geneconv.reads.sam > geneconv.reads.bam
#        samtools index geneconv.reads.bam
#        # Not using mapping quality filter because that removes noise that is needed to filter mismapping reads
#        samtools mpileup -f $refcur -q 40 -E -Q 0 -O --output-QNAME geneconv.reads.bam > geneconv.reads.pileup
#    } &
#done
#
## Get reads overlapping gene-conversion from geneconv.reads.sam and align them to ORA genome
#refora=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/ora.genome.v1.fasta
#for sample in wt7 wt18 mut4 mut11_2 MUT_11_1 MUT_15 WT_1 WT_19; do
#    cd ${CWD}/${sample}
#    bsub -q multicore20 -n 8 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${sample}_geneconv_ora_pileup.log -eo ${sample}_geneconv_ora_pileup.err "
#        samtools sort -@ 8 -n geneconv.reads.bam \
#        | samtools fastq -@ 8 -1 geneconv.reads.R1.fq.gz -2 geneconv.reads.R2.fq.gz -n -
#        bowtie2 --end-to-end \
#            --very-sensitive \
#            --threads 8 \
#            -x $refora \
#            -1 geneconv.reads.R1.fq.gz \
#            -2 geneconv.reads.R2.fq.gz \
#            --rg-id ${sample} \
#            --rg PL:ILLUMINA \
#            --rg SM:1 \
#        | samtools view -h -@8 -F4 - \
#        | samtools sort -@8 -O BAM - \
#        > geneconv.reads.ora.sorted.bam
#        samtools index geneconv.reads.ora.sorted.bam
#        samtools mpileup -f $refora -q 40 -E -Q 0 -O --output-QNAME geneconv.reads.ora.sorted.bam > geneconv.reads.ora.pileup
#    "
#done
# </editor-fold>

# <editor-fold desc="Using the same methodology for gene conversion identification in leaves as used in layers">

# Adapting python.layer_specific_dna_analysis.layer_specific_gene_conversion_all_samples for leaf_dna_analysis

# Get read counts at syri snp and indel positions
syrisnp=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt
syriind=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.indels.txt
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
for sample in wt7 wt18 mut4 mut11_2 ; do
    cd ${cwd}/${sample}
    shvs=../shv_close.for_mpileup.bed
    bsub -q multicore20 -n 10 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${sample}_bamrc_q10.log -eo ${sample}_bamrc_q10.err "
        # Bam-readcount at syri snp and indel positions
        /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -b 30 -q 10 -w 0 -I -n 10 -f $refcur -l $syrisnp ${sample}.deduped.bam ${sample}.syri_snps.bamrc
        /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -b 30 -q 10 -w 0 -I -n 10 -f $refcur -l $syriind ${sample}.deduped.bam ${sample}.syri_indels.bamrc
        # Bam-readcount at syri snp and indel positions (shvs)
        samtools mpileup -f $refcur -q 10 -E -Q 30 ${sample}.deduped.bam -l $shvs -O --output-QNAME > shvs.pileup
        sort -k1,1 -k2,2n -o shvs.pileup shvs.pileup
    "
done

for sample in MUT_11_1 MUT_15 WT_1 WT_19 ; do
    cd ${cwd}/${sample}
    shvs=../shv_close.for_mpileup.bed
    bsub -q multicore40 -n 10 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${sample}_bamrc_q10.log -eo ${sample}_bamrc_q10.err "
        # Bam-readcount at syri snp and indel positions
        /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -b 30 -q 10 -w 0 -I -n 10 -f $refcur -l $syrisnp ${sample}.sorted.bt2.bam ${sample}.syri_snps.bamrc
        /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -b 30 -q 10 -w 0 -I -n 10 -f $refcur -l $syriind ${sample}.sorted.bt2.bam ${sample}.syri_indels.bamrc
        # Bam-readcount at syri snp and indel positions (shvs)
        samtools mpileup -f $refcur -q 10 -E -Q 30 ${sample}.sorted.bt2.bam -l $shvs -O --output-QNAME > shvs.pileup
        sort -k1,1 -k2,2n -o shvs.pileup shvs.pileup
    "
done

# Get pileup data at candidate gene convergence positions in all samples
# Being used in python.leaf_dna_analysis.gene_conversion_all_samples
bedfile=${cwd}/all_sample_candidate_high_gene_convergence.bed
for sample in MUT_11_1 MUT_15 WT_1 WT_19 ; do
    cd ${cwd}/${sample}
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${sample}_pileup.log -eo ${sample}_pileup.err  "
        samtools mpileup -A -f $refcur -q 0 -E -Q 0 ${sample}.sorted.bt2.bam -l $bedfile -O --output-QNAME > all_sample_candidate_high_gene_convergence.q0.Q0.pileup
        sort -k1,1 -k2,2n -o all_sample_candidate_high_gene_convergence.q0.Q0.pileup all_sample_candidate_high_gene_convergence.q0.Q0.pileup
        samtools mpileup -A -f $refcur -q 10 -E -Q 13 ${sample}.sorted.bt2.bam -l $bedfile -O --output-QNAME > all_sample_candidate_high_gene_convergence.q10.Q13.pileup
        sort -k1,1 -k2,2n -o all_sample_candidate_high_gene_convergence.q10.Q13.pileup all_sample_candidate_high_gene_convergence.q10.Q13.pileup
    "
done
for sample in wt7 wt18 mut4 mut11_2 ; do
    cd ${cwd}/${sample}
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${sample}_pileup.log -eo ${sample}_pileup.err  "
        samtools mpileup -A -f $refcur -q 0 -E -Q 0 ${sample}.deduped.bam -l $bedfile -O --output-QNAME > all_sample_candidate_high_gene_convergence.q0.Q0.pileup
        sort -k1,1 -k2,2n -o all_sample_candidate_high_gene_convergence.q0.Q0.pileup all_sample_candidate_high_gene_convergence.q0.Q0.pileup
        samtools mpileup -A -f $refcur -q 10 -E -Q 13 ${sample}.deduped.bam -l $bedfile -O --output-QNAME > all_sample_candidate_high_gene_convergence.q10.Q13.pileup
        sort -k1,1 -k2,2n -o all_sample_candidate_high_gene_convergence.q10.Q13.pileup all_sample_candidate_high_gene_convergence.q10.Q13.pileup
    "
done



# Merge candidate gene conversion list for manual curation
cd $cwd
head -n1 wt_1/candidate_gene_conversion.txt \
| cut -f 3,7,10,11,13 --complement \
| awk '{print $0"\tsample\tselected\tploidy\tremark"}' > candidate_gene_conversion.all_sample.tsv
for s in ${samples[@]}; do
    tail +2 ${s}/candidate_gene_conversion.txt \
    | cut -f 3,7,10,11,13 --complement \
    | awk -v s=$s '{print $0"\t"s}' >> candidate_gene_conversion.all_sample.tsv
done
sort -k1,1 -k2,2n -o candidate_gene_conversion.all_sample.tsv candidate_gene_conversion.all_sample.tsv

# </editor-fold>

####################################################################
############ Step 5: Candidate gene check
####################################################################
# Get orthogroups for the flowering time genes
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/pheno_gene/
cd $CWD
grep -Ff /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/misc/arabidopsis_flower_gene_list.txt /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/orthology/cur_athal/OrthoFinder/Results_Sep29/Orthogroups/Orthogroups.txt > flower_gene.ortho.txt
grep -o  'mRNA[^ ]*' flower_gene.ortho.txt \
| awk -F '.' '{print $1"."$2}' \
| sed 's/mRNA/Gene/g' \
| sort \
| uniq \
> rp_flogen.txt
grep -w -Ff rp_flogen.txt /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/function/gene_annotations.txt > rp_flogen.annotations.txt

# Merge the Arabidopsis and Cur data
# Code in get_candidate_for_mutant_phenotype.py