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

#           GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
            /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q10.bt2.txt


        "
  done
done

# Get read counts at candidate high coverage somatic mutations in all samples
# Being used in layer_specific_dna_analysis => sm_after_masking_layers()
regfile=${CWD}/all_samples_candidate_high_cov_sms.txt
hometools=/srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
        bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${l}_bamrc.log -eo ${l}_bamrc.err  "
            $hometools pbamrc -n 1 -b 0 -q 0 -w 0 -I -f $refcur -l ${regfile} ${l}.deduped.bam bam_read_counts_b0_q0.bt2.txt
        "
  done
done

# Get pileup data at candidate high coverage somatic mutations (second set) in all samples
# Being used in layer_specific_dna_analysis => sm_after_masking_layers()
bedfile=${CWD}/all_samples_candidate_high_cov_sms2.bed
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
        bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${l}_pileup.log -eo ${l}_pileup.err  "
            samtools mpileup -A -f $refcur -q 0 -E -Q 0 ${l}.deduped.bam -l $bedfile -O --output-QNAME > all_samples_candidate_high_cov_sms2.q0.Q0.pileup
            sort -k1,1 -k2,2n -o all_samples_candidate_high_cov_sms2.q0.Q0.pileup all_samples_candidate_high_cov_sms2.q0.Q0.pileup
            samtools mpileup -A -f $refcur -q 10 -E -Q 13 ${l}.deduped.bam -l $bedfile -O --output-QNAME > all_samples_candidate_high_cov_sms2.q10.Q13.pileup
            sort -k1,1 -k2,2n -o all_samples_candidate_high_cov_sms2.q10.Q13.pileup all_samples_candidate_high_cov_sms2.q10.Q13.pileup
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


#grep 'MUT_11_1' /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions > ${CWD}/mut_11_1.mutations.regions
#
##MUTS=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions
#MUTSNEW=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.unique.regions


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
# Merge the different somatic mutations list
layer_specific_dna_analysis.py -> merge_variant_calls()


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
    hometools pbamrc -n 4 -b 0 -q 0 -w 0 -I -f $refcur -l all_layer_somatic_variants.positions ../${sample}/${sample}.sorted.bt2.bam ${sample}.all_layer_somatic_variants.read_count.txt
    echo $sample
done
for sample in 'wt7' 'wt18' 'mut4' 'mut11_2' ; do
    hometools pbamrc -n 4 -b 0 -q 0 -w 0 -I -f $refcur -l all_layer_somatic_variants.positions ../${sample}/${sample}.deduped.bam ${sample}.all_layer_somatic_variants.read_count.txt
    echo $sample
done

samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        hometools pbamrc -n 10 -b 0 -q 0 -w 0 -f $refcur -I -l all_layer_somatic_variants.positions ${s}/${s}_${l}/${l}.deduped.bam ${s}.${l}.all_layer_somatic_variants.read_count.txt
    done
done

# Get allele frequency of all_sm_in_all_sample in leaf, L1, and L2 samples
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples_readcounts/
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
muts=${indir}/all_sm_in_all_samples.regions
# Leaf scdna
cd $cwd
for s in 'WT_1' 'WT_19' 'MUT_15' 'MUT_11_1' ; do
    echo $s
    hometools pbamrc -n 10 -b 0 -q 0 -w 0 -I -f $refcur -l $muts ../${s}/${s}.sorted.bt2.bam ${s}.all_sm_in_all_samples.read_count.txt
done
# Leaf normal
for s in 'wt7' 'wt18' 'mut4' 'mut11_2' ; do
    echo $s
    hometools pbamrc -n 10 -b 0 -q 0 -w 0 -I -f $refcur -l $muts ../${s}/${s}.deduped.bam ${s}.all_sm_in_all_samples.read_count.txt
done
# Layers
for s in wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15; do
    for l in l1 l2; do
        hometools pbamrc -n 10 -b 0 -q 0 -w 0 -f $refcur -I -l $muts ../layer_samples/${s}/${s}_${l}/${l}.deduped.bam ${s}.${l}.all_sm_in_all_samples.read_count.txt
    done
done

####################################################################
############ Gene-conversion identification
####################################################################
# Initial testing with only syri snps found few gene conversions in L1
# Extending the gene conersion identification pipeline to syri indels as well
# Get read counts at syri snp and indel positions
syrisnp=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt
syriind=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.indels.txt
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${CWD}/${s}/${s}_${l}
        bsub -q multicore40 -n 10  -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${l}_shv_bamrc.log -eo ${l}_shv_bamrc.err "
#            /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -b 30 -q 10 -w 0 -I -n 10 -f $refcur -l $syrisnp ${l}.deduped.bam ${l}.syri_snps.bamrc
            /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -b 30 -q 10 -w 0 -I -n 10 -f $refcur -l $syriind ${l}.deduped.bam ${l}.syri_indels.bamrc
        "
    done
done


# Get pileup data at syri snp and indel positions
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
for s in ${samples[@]}; do
    for l in l1 l2 l3; do
        cd ${cwd}/${s}/${s}_${l}
        N=1
        shvs=../../../shv_close.for_mpileup.bed
        bsub -q ioheavy -n $N -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${s}_${l}_shv_pileup.log -eo ${s}_${l}_shv_pileup.err "
#        for i in {1..8}; do
#            snps="../../../snps_split/snps_CUR${i}G.txt"    # Split from SNPs that are close (within 1kb) of other SNPs
#            shvs="../../../shv_split/shv_CUR${i}G.txt"    # Split from SNPs that are close (within 1kb) of other SNPs
                # Split from SNPs that are close (within 1kb) of other SNPs
#                xargs -a $snps -P $N -I {} samtools mpileup -f $refcur -q 40 -E -Q 26 ${l}.deduped.bam -r {} -O --output-QNAME > snps_CUR${i}G.pileup 2> garb
#                sort -k2,2n -o snps_CUR${i}G.pileup snps_CUR${i}G.pileup
#                xargs -a $shvs -P $N -I {} samtools mpileup -f $refcur -q 10 -E -Q 30 ${l}.deduped.bam -r {} -O --output-QNAME > shvs_CUR${i}G.pileup 2> garb
                samtools mpileup -f $refcur -q 10 -E -Q 30 ${l}.deduped.bam -l $shvs -O --output-QNAME > shvs.pileup
#                sort -k2,2n -o shvs_CUR${i}G.pileup shvs_CUR${i}G.pileup
                sort -k1,1 -k2,2n -o shvs.pileup shvs.pileup
        "
    done
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

################################################################################
############ SVs identification
################################################################################

# Call SVs using Manta, GRIDSS
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'
samples=( wt_1 wt_7 wt_18 wt_19 mut_11_1 mut_11_2 mut_15 )

## Manta : I do joint-calling for all L1 samples, joint-calling for all L2 samples,
## and paired L1 vs L2 and L2 vs L1 calling for each sample.
manta_config='/srv/netscratch/dep_mercier/grp_schneeberger/software/manta-1.6.0.centos6_x86_64/bin/configManta.py'
python='/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python'
### Joint-calling L1 and L2 samples
cd $cwd
mkdir -p l1_manta
fins=''
for s in ${samples[@]}; do
    fins="${fins} --bam ${s}/${s}_l1/l1.deduped.bam"
done
bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo l1_manta.log -eo l1_manta.err "
    $python \
        $manta_config \
        ${fins} --referenceFasta $refcur --runDir l1_manta
        cd l1_manta
    $python  ./runWorkflow.py -m local -j 40 -g 50
"
cd $cwd
mkdir -p l2_manta
fins=''
for s in ${samples[@]}; do
    fins="${fins} --bam ${s}/${s}_l2/l2.deduped.bam"
done
bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo l2_manta.log -eo l2_manta.err "
    $python \
        $manta_config \
        ${fins} --referenceFasta $refcur --runDir l2_manta
        cd l2_manta
    $python  ./runWorkflow.py -m local -j 40 -g 50
"

### L1 vs L2 and L2 vs L1 paired SV calling for each sample
for s in ${samples[@]}; do
    #### Call SVs in L1 while considering L2 as normal
    cd $cwd; cd $s
    mkdir -p l1_l2_manta
    bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo l1_l2_manta.log -eo l1_l2_manta.err "
        $python \
            $manta_config \
            --normalBam ${s}_l2/l2.deduped.bam \
            --tumorBam ${s}_l1/l1.deduped.bam \
            --referenceFasta $refcur --runDir l1_l2_manta
            cd l1_l2_manta
        $python  ./runWorkflow.py -m local -j 40 -g 50
    "
    #### Call SVs in L2 while considering L1 as normal
    cd $cwd; cd $s
    mkdir -p l2_l1_manta
    bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo l2_l1_manta.log -eo l2_l1_manta.err "
        $python \
            $manta_config \
            --normalBam ${s}_l1/l1.deduped.bam \
            --tumorBam ${s}_l2/l2.deduped.bam \
            --referenceFasta $refcur --runDir l2_l1_manta
            cd l2_l1_manta
        $python  ./runWorkflow.py -m local -j 40 -g 50
    "
done
### Layer-specific calls identified using: layer_specific_dna_analysis.py -> get_layer_specific_svs()

################################################################################
############ TE insertion identification
################################################################################

# Find TE insertion/deletion sites using TEPID
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
yhindex=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.X15_01_65525S
## Get TE location and sequence
cd $cwd
grep -v "#" /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.ann.gff3 \
 | egrep -v 'Low|Simple|RNA|other|Satellite'  \
 | cut -f 1,4,5,7,9 \
 | sed 's/;/\t/g ; s/ID=//g ; s/Subfamily=//g ; s/Family=//g ; s/Class=//g' > TMP.bed
paste <(cut -f 1-6 TMP.bed) <(cut -f7- TMP.bed | tr '\t' '/') \
 > cur.TE.bed

## Map reads to the reference genome
for s in ${samples[@]}; do
    for l in l1 l2; do
        cd ${cwd}/${s}/${s}_${l}
        mkdir -p tepid; cd tepid
        bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
            echo Running tepid-map
            /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/tepid/bin/tepid-map \
                -x $refcur \
                -y $yhindex \
                -p 40 \
                -s 400 \
                -n ${s}_${l} \
                -1 ../${l}_ql-trimmed-pair1.fastq.gz \
                -2 ../${l}_ql-trimmed-pair2.fastq.gz \
                -z
        "
    done
done

## Call TE indels
for s in ${samples[@]}; do
    for l in l1 l2; do
        cd ${cwd}/${s}/${s}_${l}
        cd tepid
        bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=30000]" -M 40000 -oo ${s}_${l}.log -eo ${s}_${l}.err "
        echo Getting discordant reads
        samtools view -b -@ 40 -F 1294 ${s}_${l}.bam \
        | samtools sort -n -@ 40 -O BAM -o ${s}_${l}.discordants.bam -
        echo Running tepid discover
        /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/tepid/bin/tepid-discover \
            --strict \
            -k \
            -D ${s}_${l}.discordants.bam \
            -p 40 \
            -n ${s}_${l} \
            -c ${s}_${l}.bam \
            -s ${s}_${l}.split.bam \
            -t ../../../cur.TE.bed
        "
    done
done


################################################################################
############ MUT specific SVs identification
################################################################################
# Call somatic SVs in mutant branch considering the wt branch as the background (normal)
# For this, first I do joint-variant calling for all mut samples and separately for all
# wt samples. Then for each mut sample, do somatic mutation calling using all WT samples
# as the background.
# Phenotype causing mutations should then be present in the jointly-called and as well
# as the somatic mutation call. And ideally, should be present in multiple mutant sample.

# Call SVs using Manta
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'
wts=( wt_1 wt_7 wt_18 wt_19 )
muts=( mut_11_1 mut_11_2 mut_15 )

# Manta joint-calling for all wt samples, joint-calling for all mut samples
manta_config='/srv/netscratch/dep_mercier/grp_schneeberger/software/manta-1.6.0.centos6_x86_64/bin/configManta.py'
python='/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python'
## Joint-calling wt samples
cd $cwd
mkdir -p wt_manta; cd wt_manta
fins=''
for s in ${wts[@]}; do
    for l in l1 l2; do
        fins="${fins} --bam ../${s}/${s}_${l}/${s}_${l}.deduped.RG.withSM.bam"
    done
done
bsub -q multicore40 -n 10 -R "span[hosts=1] rusage[mem=10000]" -M 20000 -oo wt_manta.log -eo wt_manta.err "
    $python \
        $manta_config \
        ${fins} --referenceFasta $refcur --runDir .
    $python  ./runWorkflow.py -m local -j 10 -g 10
"
## Joint-calling mut samples
cd $cwd
mkdir -p mut_manta; cd mut_manta
fins=''
for s in ${muts[@]}; do
    for l in l1 l2; do
        fins="${fins} --bam ../${s}/${s}_${l}/${s}_${l}.deduped.RG.withSM.bam"
    done
done
bsub -q multicore40 -n 10 -R "span[hosts=1] rusage[mem=10000]" -M 20000 -oo mut_manta.log -eo mut_manta.err "
    $python \
        $manta_config \
        ${fins} --referenceFasta $refcur --runDir .
    $python  ./runWorkflow.py -m local -j 10 -g 10
"

# Get mut somatic mutations with the wts as normal. Manta allows only one normal
# sample per run. So, I use wt19 L1 and L2 as the normals because this is the
# closest wt branch to the mutant branches.

## For each mut sample (L1 and L2), run manta using the wt19 samples as normals
for s in ${muts[@]}; do
    for l in l1 l2; do
    cd ${cwd}/${s}/${s}_${l}/
    mkdir -p ${l}_manta; cd ${l}_manta
    bsub -q multicore40 -n 10 -R "span[hosts=1] rusage[mem=10000]" -M 20000 -oo ${l}_manta.log -eo ${l}_manta.err "
        $python \
            $manta_config \
            --normalBam ../../../wt_19/wt_19_${l}/wt_19_${l}.deduped.RG.withSM.bam\
            --tumorBam ../${s}_${l}.deduped.RG.withSM.bam \
            --referenceFasta $refcur --runDir .
        $python  ./runWorkflow.py -m local -j 10 -g 10
    "
    done
done

### Layer-specific calls identified using: layer_specific_dna_analysis.py -> get_mut_branch_specific_svs()

################################### OLD CODE ###################################
## SV identification using gridss

## GRIDSS : I do joint-calling for all L1 samples, joint-calling for all L2 samples,
## and paired L1 vs L2 and L2 vs L1 calling for each sample (similar pipeline as Manta)
gridssjar=/srv/netscratch/dep_mercier/grp_schneeberger/software/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar
### GRIDSS STEP 1: setupreference --> This creates BWA mem index. Not required as all files are already there
### GRIDSS STEP 2: Joint-calling L1 and L2 samples
#### Create bam files with correct RG
for s in ${samples[@]}; do
    for l in l1 l2; do
        cd ${cwd}/${s}/${s}_${l}
#        bsub -q ioheavy -n 4 -R "span[hosts=1] rusage[mem=32000]" -M 40000 -oo ${s}_${l}_RG.log -eo ${s}_${l}_RG.err "
            #TODO: The SM tag is used by the downstream methods. It would be better to change it to <branch>_<layer> rather than just 1
            samtools addreplacerg -@ 4 -r ID:${s}_${l}'\t'PL:ILLUMINA'\t'SM:${s}_${l} --no-PG -o ${s}_${l}.deduped.RG.withSM.bam ${l}.deduped.bam &
#        "
    done
done

### Joint-calling L1 and L2 samples
#### L1 call
cd $cwd
mkdir -p l1_gridss; cd l1_gridss
#### Get file lists for L1 jobs
bam_fins=''
labels=''
for s in ${samples[@]}; do
    bam_fins="${bam_fins} ../${s}/${s}_l1/${s}_l1.deduped.RG.bam"
    labels="${labels} ${s}_l1"
done
labels=$(echo $labels | tr ' ' ',')
#### Run gridss for L1 joint calling
bsub -q ioheavy -n 8 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo l1_gridss.log -eo l1_gridss.err "
    export PATH=/opt/share/software/packages/samtools-1.16.1/bin/:$PATH
    gridss \
        -s preprocess,assemble,call \
        -j $gridssjar \
        -r $refcur \
        -o l1.gridss.vcf \
        --labels $labels \
        -t 8 \
        $bam_fins
"
#### L2 call
cd $cwd
mkdir -p l2_gridss; cd l2_gridss
#### Get file lists for L2 jobs
bam_fins=''
labels=''
for s in ${samples[@]}; do
    bam_fins="${bam_fins} ../${s}/${s}_l2/${s}_l2.deduped.RG.bam"
    labels="${labels} ${s}_l2"
done
labels=$(echo $labels | tr ' ' ',')
#### Run gridss for L2 joint calling
bsub -q ioheavy -n 8 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo l2_gridss.log -eo l2_gridss.err "
    export PATH=/opt/share/software/packages/samtools-1.16.1/bin/:$PATH
    gridss \
        -s preprocess,assemble,call \
        -j $gridssjar \
        -r $refcur \
        -o l2.gridss.vcf \
        --labels $labels \
        -t 8 \
        $bam_fins
"

### Pairwise SV calling L1 vs L2
#### Run gridss on L1 and L2 separately. This will be required to create panel of normals
for s in ${samples[@]}; do
    for l in l1 l2; do
        cd ${cwd}/${s}/${s}_${l}
        mkdir -p gridss; cd gridss
        bsub -q ioheavy -n 8 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${l}_gridss.log -eo ${l}_gridss.err "
            export PATH=/opt/share/software/packages/samtools-1.16.1/bin/:$PATH
            gridss \
                -s preprocess,assemble,call \
                -j $gridssjar \
                -r $refcur \
                -o ${l}.gridss.vcf \
                -t 8 \
                ../${s}_${l}.deduped.RG.bam
        "
    done
done
#### Create PON for L1
cd $cwd
mkdir -p l1_gridss_pon; cd l1_gridss_pon
input=$(ls -1 ../*/*_l1/gridss/*.vcf | awk ' { print "INPUT=" $0 }')
mkdir -p pondir
java -Xmx8g \
    -cp $gridssjar \
    gridss.GeneratePonBedpe \
    $input \
    O=pondir/gridss_pon_breakpoint.bedpe \
    SBO=pondir/gridss_pon_single_breakend.bed \
    REFERENCE_SEQUENCE=$refcur \
    NORMAL_ORDINAL=0
#### Create PON for L2
cd $cwd
mkdir -p l2_gridss_pon; cd l2_gridss_pon
input=$(ls -1 ../*/*_l2/gridss/*.vcf | awk ' { print "INPUT=" $0 }')
mkdir -p pondir
java -Xmx8g \
    -cp $gridssjar \
    gridss.GeneratePonBedpe \
    $input \
    O=pondir/gridss_pon_breakpoint.bedpe \
    SBO=pondir/gridss_pon_single_breakend.bed \
    REFERENCE_SEQUENCE=$refcur \
    NORMAL_ORDINAL=0


#### Call pairwise somatic mutations between layers
for s in ${samples[@]}; do
    #### Call SVs in L1 while considering L2 as normal
    cd $cwd; cd $s
    mkdir -p l1_l2_gridss; cd l1_l2_gridss
    bsub -q ioheavy -n 8 -R "span[hosts=1] rusage[mem=10000]" -M 20000 -oo l1_l2_gridss.log -eo l1_l2_gridss.err "
        export PATH=/opt/share/software/packages/samtools-1.16.1/bin/:$PATH
        gridss \
            -s preprocess,assemble,call \
            -j $gridssjar \
            -r $refcur \
            -o l1_l2.gridss.vcf \
            --labels ${s}_l2,${s}_l1 \
            -t 8 \
            ../${s}_l2/${s}_l2.deduped.RG.bam \
            ../${s}_l1/${s}_l1.deduped.RG.bam \

        /srv/netscratch/dep_mercier/grp_schneeberger/software/gridss/gridss_somatic_filter_mg \
          --pondir ./pondir/ \
          --input wt_7_l2.vcf \
          --output high_confidence_somatic.vcf.gz \
          --fulloutput high_and_low_confidence_somatic.vcf.gz \
          --scriptdir /srv/netscratch/dep_mercier/grp_schneeberger/software/gridss/ \
          --ref BSgenome.Parmeniaca.LMU.cur \
          -n 1 \
          -t 2

    "
    #### Call SVs in L2 while considering L1 as normal
    cd $cwd; cd $s
    mkdir -p l2_l1_gridss; cd l2_l1_gridss
    bsub -q ioheavy -n 8 -R "span[hosts=1] rusage[mem=10000]" -M 20000 -oo l2_l1_gridss.log -eo l2_l1_gridss.err "
        export PATH=/opt/share/software/packages/samtools-1.16.1/bin/:$PATH
        gridss \
            -s preprocess,assemble,call \
            -j $gridssjar \
            -r $refcur \
            -o l2_l1.gridss.vcf \
            --labels ${s}_l1,${s}_l2 \
            -t 8 \
            ../${s}_l1/${s}_l1.deduped.RG.bam \
            ../${s}_l2/${s}_l2.deduped.RG.bam \
    "
done

## TE insertions identification using TEfinder
# Find TE insertion/deletion sites using TEfinder
## Find TE insertion/deletion sites using TEfinder
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
RPout=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.out.gff
#INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'

cd $cwd
grep -v -iE '(Motif\:[ATGC]+\-rich)|(Motif\:\([ATGC]+\)n)' $RPout > TEs.gtf
awk -F '\t' '{print $9}' TEs.gtf | awk -F '"' '{print $2}' | sort | uniq > List_of_TEs.txt

for s in ${samples[@]}; do
    for l in l1 l2; do
        cd ${cwd}/${s}/${s}_${l}
        mkdir -p tefinder; cd tefinder
        bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=40000]" -m "hpc001 hpc003 hpc004 hpc005 hpc006 hpc007" -M 100000 -oo ${s}.log -eo ${s}.err "
            # TEfinder_mg has minor modifications to ensure proper parallelisation
            rm -r ./tmp
            TEfinder_mg2 \
                -picard $picard \
                -k 30 \
                -intermed yes \
                -maxHeapMem 20000 \
                -threads 40 \
                -workingdir ./tmp \
                -alignment ../${s}_${l}.deduped.RG.bam \
                -fa $refcur\
                -gtf ../../../TEs.gtf \
                -te ../../../List_of_TEs.txt
        "
    done
done
