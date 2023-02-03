NOTE: FOR ALL ANALYSIS CURROT assembly would be used as the primary assembly

################################################################################
############# Pre-run: Get repeat regions in the currot assembly ###############
################################################################################
### Not required as done in the scdna_analysis.sh

################################################################################
############################ Pre-process reads #################################
################################################################################
hometools='/srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/python /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py'

### Get sequencing quality plots
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/fruit_layer_specific/
samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
cd $CWD
for s in ${samples[@]}; do
	for l in l1 l2 l3; do
		cd ${CWD}/${s}/${s}_${l}
		fastqc -t 10 *fq.gz &
	done
done

### Merge, trim and align reads to CUR genomes

CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/fruit_layer_specific/
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
adpx=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
adpy=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG
samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
samples=( wt_7 wt_18 )

for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd $CWD
        mkdir -p ${s}/${s}_${l}
        cd ${s}/${s}_${l}
        bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo ${s}_${l}.log -eo ${s}_${l}.err "
            cat ${indir}/${s}/${s}_${l}/*_1.fq.gz > merged_read_1.fq.gz
            cat ${indir}/${s}/${s}_${l}/*_2.fq.gz > merged_read_2.fq.gz
            fastqc -t 20 merged_read_1.fq.gz merged_read_2.fq.gz
            skewer -r 0.1 -d 0.05 -k 8 -q 20 -l 75 -m pe -t 20 -x $adpx -y $adpy \
                merged_read_1.fq.gz \
                merged_read_2.fq.gz \
                -z -o ${l}_ql

            bowtie2 --end-to-end \
                --very-sensitive \
                --threads 20 \
                -x $curidxbt2 \
                -1 ${l}_ql-trimmed-pair1.fastq.gz \
                -2 ${l}_ql-trimmed-pair2.fastq.gz \
                --rg-id ${s}_${l}\'\\t\'PL:ILLUMINA\'\\t\'SM:1 \
            | samtools view -b - \
            > ${l}.bam

            ## Mark duplicated reads
            samtools sort -@ 10 -n -O BAM ${l}.bam \
            | samtools fixmate -@ 10 -c -m -O BAM - - \
            | samtools sort -@ 10 -O BAM - \
            | samtools markdup -@ 10 -S -s -O BAM - - \
            | samtools sort -@ 10 -O BAM - > ${l}.DUPmarked.bam
            samtools index ${l}.DUPmarked.bam

            samtools view -G 1024 -O BAM ${l}.DUPmarked.bam > ${l}.deduped.bam
            samtools index ${l}.deduped.bam
        "    
    done
done
    

# Copied the trimmed-read and alignments to the following folder and would now use that as CWD
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
#samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
samples=( wt_7 wt_18 )
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
        bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=15000]" -M 20000 -oo ${l}_cov.log -eo ${l}_cov.err "
            $hometools bamcov -g ${l}.deduped.bam ${l}_vs_cur.deduped.cov
        "
    done
done


####################################################################
############ Get read-counts at variant positions
####################################################################
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
CHRBED=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.chrlen.bed
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
        bsub -q ioheavy -n 8 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${l}_bamrc.log -eo ${l}_bamrc.err  "
            $hometools pbamrc -n 8 -b 30 -q 10 -w 0 -S -I -f $refcur -l $CHRBED ${l}.deduped.bam bam_read_counts_b30_q10.bt2.txt

          # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
          /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q10.bt2.txt
        "
  done
done

## Get read mapping depth histogram, allele count at leaf mutation positions
MUTSNEW=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.unique.regions
samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
#samples=( wt_7 wt_18 )
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
#        cut -f 4 bam_read_counts_b30_q10.bt2.txt \
#        | sort -n \
#        | uniq -c \
#        | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.pdf -x Mapping Depth -y Frequency -t ${s}_${l} -xlim 0 400 &
#        ## Get allele-frequency at the updated leaf mutation positions
#        hometools pbamrc -b 30 -q 10 -w 0 -f $refcur -I -n 1 -l $MUTSNEW ${l}.deduped.bam read_count_at_leaf_mut_pos.all_branches.updated.txt &
        hometools pbamrc -b 0 -q 0 -w 0 -f $refcur -I -n 1 -l $MUTSNEW ${l}.deduped.bam read_count_at_leaf_mut_pos.all_branches.b0_q0.updated.txt &
    done
done




#
#
#grep 'MUT_11_1' /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions > ${CWD}/mut_11_1.mutations.regions
#
##MUTS=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions
#MUTSNEW=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.unique.regions
#for s in l1 l2 l3; do
#    cd ${CWD}/mut_11_1_${s}
#
#    # Get allele-frequency at the mutation positions in MUT_11_1
##    hometools pbamrc -b 30 -q 10 -w 0 -f $refcur -I -n 1 -l ../mut_11_1.mutations.regions ${s}.deduped.bam read_count_at_leaf_mut_pos.txt &
#
##    Get allele-frequency at the updated mutation positions
#    hometools pbamrc -b 30 -q 10 -w 0 -f $refcur -I -n 1 -l $MUTSNEW ${s}.deduped.bam read_count_at_leaf_mut_pos.all_branches.updated.txt &
#
#done
#



####################################################################
############ De-Novo mutation identification
####################################################################
# Get Layer specific somatic mutations:
layer_specific_dna_analysis.py -> layer_specific_sm_calling_all_samples()
# Get Layer 3 specific somatic mutations:
layer_specific_dna_analysis.py -> layer_3_variant_calling()
# Get Layer specific somatic mutations conserved in multiple samples:
layer_specific_dna_analysis.py -> layer_conserved_variants()
# Get layer specific somatic mutations using only fold-change cutoff and without leaf background noise removal:
layer_specific_dna_analysis.py -> layer_specific_fc_check()


####################################################################
#### Statistics
####################################################################
# Get SNP positions
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples
cat */high_conf_layer_specific_somatic_mutations.selected.tsv  | grep Y | sort -k1,1 -k2,2 | grep -v -E '[+-]' |cut -f 1,2,3,4 | uniq > SNP_positions.txt

# For somatic mutations identified in leaves
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/plot_mutation_changes_spectra.py \
    -f manual_selected_somatic_mutations.csv \
    -rc 3 \
    -qc 4 \
    -t Nucleotide substitute \
    -W 8 \
    -ymax 40 \
    -o somatic_snps_substitution.png &

# Get allele frequency at layer specific variant sites in the leaf data
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples
tail +2 all_layer_somatic_variants.txt  | cut -f 1,2 | uniq | sort -k1,1 -k2,2n | awk '{print $1,$2,$2}' > all_layer_somatic_variants.positions
tail +2 all_layer_somatic_variants.txt  | cut -f 1,2 | uniq | sort -k1,1 -k2,2n | awk '{print $1,$2-1,$2}' > all_layer_somatic_variants.bed


for sample in 'WT_1' 'WT_19' 'MUT_15' 'MUT_11_1' ; do
    hometools pbamrc -n 4 -b 0 -q 0 -w 0 -I -f $refcur -l all_layer_somatic_variants.positions ../${sample}/${sample}.sorted.bt2.bam ${sample}.all_layer_somatic_variants.read_count.txt &
    echo $sample
done
for sample in 'wt7' 'wt18' 'mut4' 'mut11_2' ; do
    hometools pbamrc -n 4 -b 0 -q 0 -w 0 -I -f $refcur -l all_layer_somatic_variants.positions ../${sample}/${sample}.deduped.bam ${sample}.all_layer_somatic_variants.read_count.txt &
    echo $sample
done

samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        hometools pbamrc -b 0 -q 0 -w 0 -f $refcur -I -n 1 -l all_layer_somatic_variants.positions ${s}/${s}_${l}/${l}.deduped.bam ${s}.${l}.all_layer_somatic_variants.read_count.txt &
    done
done



####################################################################
############ Gene-conversion identification
####################################################################
# Get read counts at syri snp positions
syrisnp=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
        bsub -q multicore40 -n 10  -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${l}_snp_bamrc.log -eo ${l}_snp_bamrc.err "
            /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -b 30 -q 10 -w 0 -I -n 10 -f $refcur -l $syrisnp ${l}.deduped.bam ${l}.syri_snps.bamrc
        "
    done
done

# Get pileup data at syri snp positions
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
        N=20
        for i in {1..8}; do
            snps="../../../snps_split/snps_CUR${i}G.txt"    # Split from SNPs that are close (within 1kb) of other SNPs
            bsub -q multicore40 -n $N -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${s}_${l}_${i}_snps_pileup.log -eo ${s}_${l}_${i}_snps_pileup.err "
                xargs -a $snps -P $N -I {} samtools mpileup -f $refcur -q 40 -E -Q 26 ${l}.deduped.bam -r {} -O --output-QNAME > snps_CUR${i}G.pileup 2> garb
                sort -k2,2n -o snps_CUR${i}G.pileup snps_CUR${i}G.pileup
            "
        done
    done
done
