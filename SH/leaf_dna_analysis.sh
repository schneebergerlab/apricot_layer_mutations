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

      # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
      /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b0_q0.bt2.txt
  "
done


## Get read mapping depth histogram
for s in wt7 wt18 mut4 mut11_2; do
    cd ${CWD}/$s
    cut -f 4 bam_read_counts_b30_q10.bt2.txt | sort -n | uniq -c | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.pdf -x Mapping Depth -y Frequency -t ${s}_bt2 -xlim 0 400 &
done




####################################################################
############ Step 4: Get mutated positions (for SNPs)
####################################################################
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")

## Step 4a: Get candidate variant positions using heuristic cut-offs
cd $cwd
## Use python3.7 environment

python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_candidates_variant_positions.py \
    WT_1/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt \
    WT_19/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt \
    MUT_11_1/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt \
    MUT_15/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt \
    -s WT_1_bt2 WT_19_bt2 MUT_11_1_bt2 MUT_15_bt2 \
    -n 5 -m 50 -M 250 &

## Step 4b: Select candidates supported by multiple cells
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")

for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
    barcodes_list=${indir}${sample}/barcodes_list
    bedbt2=$(rf ${sample}_bt2_SNPs_candidate.regions)
    inpath=${indir}${sample}/barcodes/
    mkdir cells_readcount
    cd cells_readcount

    bsub -q multicore40 -n 25 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo rc2.log -eo rc2.err "
        xargs -a $barcodes_list \
        -P 25 \
        -I {} \
        bash -c '
            bc=\$(basename \${1})
            mkdir \$bc; cd \$bc
            echo \$2 \$3 \${1}/\${bc}.DUPmarked.deduped.bam \${bc}_readcount.bt2.txt
            /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/get_readcounts_in_cells_bt2.sh \$2 \$3 \${1}/\${bc}.DUPmarked.deduped.bam \${bc}_readcount.bt2.txt
            cd ..
        ' -- {} $bedbt2 $refcur
    "
done

## Use python3.7 environment
for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
    echo $sample
    rf cells_readcount/*/*_readcount.bt2.txt >bt2_rcfiles.txt

    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/select_candidate_supported_by_cells3.py \
        ${sample}_bt2_filtered_SNPs_candidate.sorted.bed \
        bt2_rcfiles.txt \
        -o multi_cell_bt2 \
        -n 5 &
done

## Step 4c: Select good candidates based on ALT allele depth and frequency comparisons in other samples
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_read_counts_of_candidates_in_other_samples.py
cd $cwd
cat */*good_candidates.txt | awk '{if($6>=10) print}' >high_cov_mutants.txt
cat */*good_candidates.txt | awk '{if($6>=20) print}' >high_cov_mutants_AF20.txt
sort -k1,1 -k2,2n -k3,3n high_cov_mutants.txt > high_cov_mutants_sorted.txt
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/plot_mutation_changes_spectra.py \
    -f high_cov_mutants.txt \
    -rc 4 \
    -qc 5 \
    -samples all_samples_BT2 \
    -t Mutation changes in $sample \
    -W 8 \
    -ymax 40 \
    -o mutation_changes_${sample}.png &

python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/plot_mutation_changes_spectra.py \
    -f high_cov_mutants_AF20.txt \
    -rc 4 \
    -qc 5 \
    -samples all_samples_BT2 \
    -t Mutation changes in $sample \
    -W 8 \
    -ymax 40 \
    -o mutation_changes_all_sample_AF20.png &


# Filter candidates with low sequencing depth based on number of neighboring SNPs, mapping quality, and BAQ
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_noisy_candidates.py

# Get mutations spectra for the final-list of SNPs
manual_selected_somatic_mutations.csv
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/plot_mutation_changes_spectra.py \
    -f manual_selected_somatic_mutations.csv \
    -rc 3 \
    -qc 4 \
    -t Nucleotide substitute \
    -W 8 \
    -ymax 40 \
    -o somatic_snps_substitution.png &
