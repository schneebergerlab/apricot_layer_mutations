NOTE: FOR ALL ANALYSIS CURROT assembly would be used as the primary assembly

################################################################################
############# Pre-run: Get repeat regions in the currot assembly ###############
################################################################################
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/
bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo repeat.log -eo repeat.err "
  mkdir repeat_region
  cd repeat_region
  BuildDatabase -name currot ../cur.genome.v1.fasta
  RepeatModeler -database currot -pa 10 -LTRStruct

  # Find repeatitive regions using RepeatMasker
  cd ..
  RepeatMasker \
    -pa 40 \
    -s \
    -norna \
    -lib ./repeat_region/currot-families.fa \
    -frag 100000 \
    -dir ./repeat_region/ \
    -poly \
    -u \
    cur.genome.v1.fasta
"
cd repeat_region
tail +4 cur.genome.v1.fasta.out | awk '{print $5"\t"$6-1"\t"$7}' >cur.genome.v1.fasta.out.bed
## Add buffer of 20 BP around indel positions
cat cur.genome.v1.fasta.out.bed |
    awk '{printf $1"\t"} {if($2<20) {printf 0"\t"} else {printf $2-20"\t"}} {print $3+20}' \
        >cur.genome.v1.fasta.out.padded_20.bed
cd ..
ln -s repeat_region/cur.genome.v1.fasta.out.padded_20.bed cur.genome.v1.repeats.bed

# Divide the genome in regions of length 1million
hometools genome_ranges cur.genome.v1.fasta

# Use genmap to get regions which have low mappability (high uniqueness) in the genomes
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/
genmap index -F cur.genome.v1.fasta -I cur.genome.v1.fasta.genmap.idx
genmap map -I cur.genome.v1.fasta.genmap.idx/ -O cur.genome.v1.fasta.genmap.E0_K51.map -E 0 -K 51 -w -bg -T 60
genmap map -I cur.genome.v1.fasta.genmap.idx/ -O cur.genome.v1.fasta.genmap.E1_K51.map -E 1 -K 51 -w -bg -T 60
cat cur.genome.v1.fasta.genmap.E0_K51.map.bedgraph | awk '{if($4==1) print $0}' >cur.genome.v1.fasta.genmap.E0_K51.map.unique.bedgraph

# List of contigs in the currot genome
grep '>' cur.genome.v1.fasta | sed 's/>//g' >cur.genome.v1.fasta.contig_list

# Create genome dict for GATK
java -jar /srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar CreateSequenceDictionary R=cur.genome.v1.fasta O=cur.genome.v1.dict

# Create blast database
makeblastdb -dbtype nucl -in cur.genome.v1.fasta -input_type fasta -title cur.genome.v1.fasta
makeblastdb -dbtype nucl -in ora.genome.v1.fasta -input_type fasta -title ora.genome.v1.fasta


################################################################################
###################### RUN CELLRANGER TO CORRECT BARCODE #######################
################################################################################
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/
cd $indir
cellranger-dna mkref cur.genome.v1.fasta currot_contig_defs.json

## Step 2: Use the indexed genomes and run the cellranger-dna CNV identification pipeline. One of the internediate steps for this pipeline is to barcode correction. We can use the final bam file to get the corrected barcodes.
## Currot genome would be used as the reference genome for all analysis (preliminary testing showed that there are not major differences between currot vs orangered)

run_barcode_correction() {
    bsub -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=100000]" -M 100000 -oo $1.log -eo $1.err "
    export _JAVA_OPTIONS=-Xgcthreads4
        cellranger-dna cnv --id=$1 --reference=$2 --fastq=$3 --sample=$4 --localcores=40 --localmem=98
        "
}
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/refdata-cur.genome.v1/'
indir="/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/"
python="/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.7/bin/python"
filter="/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_short_reads.py"

samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
<<Useful
##  All reads smaller than 50bps were filtered out as cellranger-dna cannot work with smaller reads and there were naked barcode sequenced as reads. (Question: why did naked barcode ended up getting sequenced?)
for sample in ${samples[@]}; do
  cd ${indir}/${sample}/
  R1s=($(ls ${sample}*R1*))
  R2s=($(ls ${sample}*R2*))
  cnt=$(echo ${R2s[@]} | wc -w)
  for ((i=0;i<cnt;i++)); do
    bsub -q normal -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo ${i}.log -eo ${i}.err "
      $python $filter ${R1s[${i}]} ${R2s[$i]} 50
    "
  done
done
Useful

## Cellranger pipeline on all samples
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/cellranger_out'
for sample in ${samples[@]}; do
    cd $cwd
    R2s=($(ls ${indir}/${sample}/renamed_L50_${sample}*R2* | xargs -n 1 basename))
    pres=''
    for R2 in ${R2s[@]}; do
        pres=$(echo ${pres}${R2%_S*},)
    done
    pres=$(echo ${pres} | sed 's/,$//g')
    run_barcode_correction $sample $refcur ${indir}$sample/ $pres
done

## Step 3: Separate reads from the corrected barcode bam file to separate folders corresponding to each barcode, get read count, and align them to the reference
T10xbam2fq=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/binaries/T10xbam2fq
asCellseparator=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/binaries/asCellseparator
curidx=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.mm2_Xsr.idx
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
oraidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/ora.genome.v1.fasta'

MIN_RC=10000

get_cell_reads() {
    bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 20000 -oo get_cell_reads.log -eo get_cell_reads.err "
            samtools sort -@ 19 -n $1 -o ${2}_RNsorted_bam.bam
            samtools view -@ 19 ${2}_RNsorted_bam.bam | $T10xbam2fq - ${2}
            mkdir barcodes
            $asCellseparator 16 ${2}_fqfrom10xBam_bxCorrected_R1.fq.gz ${2}_fqfrom10xBam_bxCorrected_R2.fq.gz $MIN_RC 10000000 barcodes/
            "
}

merge_fastqs() {
    bsub -q normal -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo merge_fastqs.log -eo merge_fastqs.err "
        readlink -f barcodes/* > barcodes_list
        while read r; do
            cd \$r
            bc=\$(basename \$r)
            cat part*R1*.fq.gz > \${bc}_R1.fastq.gz
            cat part*R2*.fq.gz > \${bc}_R2.fastq.gz
        done < barcodes_list

        cd ${2}/${1}
        rm barcodes_read_count
        while read r; do
            bc=\$(basename \$r)
            cnt=\$(zgrep -c '@' \$r/\${bc}_R1.fastq.gz)
            echo \$bc \$cnt >> barcodes_read_count
        done < barcodes_list
        "
}

samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/cellranger_out/'

for sample in ${samples[@]}; do
#for sample in WT_1; do
    cd $cwd/
    mkdir $sample
    cd $sample
    #    get_cell_reads ${indir}/${sample}/outs/possorted_bam.bam $sample
    #    merge_fastqs $sample $cwd
#    rm align_cells.log align_cells.err
#    nohup /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/align_cells.sh $sample $curidxbt2 $oraidxbt2 &
    bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=40000]" -M 50000 -oo align_cells.log -eo align_cells.err "
        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/align_cells.sh $sample $curidxbt2 $oraidxbt2
  "
#  break
done

for sample in ${samples[@]}; do
    {
        cd ${cwd}/${sample}
        rm barcodes_read_count_deduped
        while read r; do
            bc=$(basename $r)
            cnt=$(zgrep -c '@' ${r}/${bc}_dedup_R1.fastq.gz)
            echo $bc $cnt >>barcodes_read_count_deduped
        done <barcodes_list
    } &
done

####################################################################
############ Step  3: Aligned pooled data reads and get read-counts at variant positions
####################################################################

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
curidx='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.mm2_Xsr.idx'
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
CHRBED=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.chrlen.bed

cd $cwd

## Instead of getting reads and realigning them, I can just merge the read alignments of individual cells
## Step 3a: Merge cell-wise BAM Files
for sample in ${samples[@]}; do
    cd $cwd
    mkdir $sample
    cd $sample
    rf ${indir}${sample}/barcodes/*/*.DUPmarked.deduped.bam >cell_bam_path.txt
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo align_reads_bt2.log -eo align_reads_bt2.err "
      
#      samtools merge -O BAM -@40 -f -b cell_bam_path.txt ${sample}.sorted.bt2.bam
#      samtools index -@40 ${sample}.sorted.bt2.bam
    "
done

# Manually update RG if required for GATK
for sample in ${samples[@]}; do
    cd ${cwd}/${sample}
    samtools view -h -@40 ${sample}.sorted.bt2.bam | awk '{if($1=="@RG"){print $0"\tPL:ILLUMINA\tSM:1"} else {print $0}}' | samtools view -b -@40 > ${sample}.sorted.RG.bt2.bam &
    mv ${sample}.sorted.RG.bt2.bam ${sample}.sorted.bt2.bam
    samtools index -@ 50 ${sample}.sorted.bt2.bam &
done


for sample in ${samples[@]}; do
    cd $cwd
    mkdir $sample
    cd $sample
#    rf ${indir}${sample}/barcodes/*/*.DUPmarked.deduped.bam > cell_bam_path.txt
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo align_reads_bt2.log -eo align_reads_bt2.err "
      # bam-readcount -b 30 -q 10 -w 0 -f $refcur ${sample}.sorted.bt2.bam | awk '{printf \$1\" \"\$2\" \"\$3\" \"\$4; for(i=6;i<=10;i++) {n1=split(\$i,a,\":\"); printf \" \"a[2]};  for(i=11;i<=NF;i++) {n1=split(\$i,a,\":\"); printf \" \"a[1]\" \"a[2]}; printf \"\n\"}' > bam_read_counts_b30_q10.bt2.txt

      # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
      /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q10.bt2.txt
  "
done

for s in ${samples[@]}; do
    cd ${cwd}/$s
    bsub -q ioheavy -n 8  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc_q0.log -eo ${s}_bamrc_q0.err "
        $hometools pbamrc -n 8 -b 0 -q 0 -w 0 -S -I -f $refcur -l $CHRBED ${s}.sorted.bt2.bam bam_read_counts_b0_q0.bt2.txt

      # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
      /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b0_q0.bt2.txt
  "
done

## Step 3d: Get read mapping depth histogram
for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
#    cut -d' ' -f 4 bam_read_counts_b30_q10.bt2.txt | sort -n | uniq -c | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.pdf -x Mapping Depth -y Frequency -t ${sample}_bt2 -xlim 0 400 &
    cut -d' ' -f 4 bam_read_counts_b30_q10.bt2.txt | sort -n | uniq -c | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.png -x Mapping Depth -y Frequency -t ${sample} -xlim 0 400 -n 400 &
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

####################################################################
############ Step 5: Get mutated positions (for Indels)
####################################################################
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
cd $cwd

# Step 5a: Get candidate indels
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/somatic_indel_identification.py \
    candidates \
    WT_1/bam_read_counts_b30_q10.bt2.txt \
    WT_19/bam_read_counts_b30_q10.bt2.txt \
    MUT_11_1/bam_read_counts_b30_q10.bt2.txt \
    MUT_15/bam_read_counts_b30_q10.bt2.txt \
    -s WT_1_bt2 WT_19_bt2 MUT_11_1_bt2 MUT_15_bt2 \
    -n 5 -m 50 -M 250 --cores 4 &

## Step 5b: Select candidates indels supported by multiple cells
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")

for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
    barcodes_list=${indir}${sample}/barcodes_list
    bedbt2=$(rf ${sample}_bt2_indel_candidate.regions)
    inpath=${indir}${sample}/barcodes/
    mkdir cells_readcount_indels
    cd cells_readcount_indels

    bsub -q bigmem -n 40 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo rc2.log -eo rc2.err "
        xargs -a $barcodes_list \
        -P 40 \
        -I {} \
        bash -c '
            bc=\$(basename \${1})
            mkdir \$bc; cd \$bc
            /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/get_readcounts_in_cells_bt2_indels.sh \$2 \$3 \${1}/\${bc}.DUPmarked.deduped.bam \${bc}_readcount.bt2.txt 10
            cd ..
        ' -- {} $bedbt2 $refcur
    "
done

for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
    echo $sample
    rf cells_readcount_indels/*/*_readcount.bt2.txt > bt2_indels_rcfiles.txt
    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/somatic_indel_identification.py \
    multicell \
    ${sample}_bt2_indel_candidate.sorted.filtered.bed \
    bt2_indels_rcfiles.txt \
    -o multi_cell &
done

## Step 5c: Select good candidates based on ALT allele depth and frequency comparisons in other samples
cd $cwd
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/somatic_indel_identification.py \
    filterbg \
    --cores 4

awk '{if($6>=10) print $0"\tWT_1"}' WT_1/*_indel_good_candidates.txt > high_conf_indels_TMP.csv
awk '{if($6>=10) print $0"\tWT_19"}' WT_19/*_indel_good_candidates.txt >> high_conf_indels_TMP.csv
awk '{if($6>=10) print $0"\tMUT_11_1"}' MUT_11_1/*_indel_good_candidates.txt >> high_conf_indels_TMP.csv
awk '{if($6>=10) print $0"\tMUT_15"}' MUT_15/*_indel_good_candidates.txt >> high_conf_indels_TMP.csv

## Step 5d: Manually select indels from high candidates supported by >=10 reads. For indels supported by <10 reads, do additional filetering similar to SNP filtering.
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling

# Filter indels supported by <10 reads. This script might not work in one click
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_noisy_indel_candidates.py

# Manually curated list of indels
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_conf_indels_manually_selected.csv
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_good_candidate_indels.selected_manually.bed


# Get read-count of somatic mutations in scRNA-data
cat mutations.txt | awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > mutations.regions
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells/
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling
    hometools pbamrc -q 10 -b 30 -l mutations.regions -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -w 0 -n 10 ${INDIR}/${s}/${s}/outs/possorted_genome_bam.bam ${s}_rna_mutations.bam_rc.txt &
done
rm snps_expressions_in_rna.min_snp_af0.1.txt snps_expressions_in_rna.max_snp_af0.1.txt
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    echo $s >> snps_expressions_in_rna.min_snp_af0.1.txt
    echo $s >> snps_expressions_in_rna.max_snp_af0.1.txt
    head -39 mutations.regions | awk '{if($7>0.1){print $0}}' | grep $s$ | cut -f2 | grep -f /dev/stdin  ${s}_rna_mutations.bam_rc.txt >> snps_expressions_in_rna.min_snp_af0.1.txt
    head -39 mutations.regions | awk '{if($7<0.1){print $0}}' | grep $s$ | cut -f2 | grep -f /dev/stdin  ${s}_rna_mutations.bam_rc.txt >> snps_expressions_in_rna.max_snp_af0.1.txt
done


# Combine SNP and Indel somatic mutations into one output file

cat manual_selected_somatic_mutations.csv > somatic_snvs.txt
cut -f 1,3,4,5,6,7,8 manual_selected_low_cov_somatic_mutations.txt >> somatic_snvs.txt
cut -f 1,3,4,5,6,7,8 high_conf_indels_manually_selected.csv > somatic_indels.txt
cut -f 1,3,4,5,6,7,8 all_good_candidate_indels.selected_manually.bed >> somatic_indels.txt
cat somatic_snvs.txt somatic_indels.txt > mutations.txt

### Get neighboring bases at the SNV positions
head -39 mutations.bed \
| awk '{print $1"\t"$2-1"\t"$3+1}' \
| hometools exseq --f /dev/stdin /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -o snps_with_neighbour_sequence.fa
################################################################################
############ STEP 6: Mitotic recombinations in individual cells  ###############
################################################################################
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")

cd $cwd
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/strict_syntenic_markers.py /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out strict_syn_snp.txt
awk '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}' strict_syn_snp.txt > strict_syn_snp.bed

samtools merge -@80 -O BAM merged_samples.bam ../WT_1/WT_1.sorted.bt2.bam ../WT_19/WT_19.sorted.bt2.bam ../MUT_11_1/MUT_11_1.sorted.bt2.bam ../MUT_15/MUT_15.sorted.bt2.bam
samtools view -C -T $refcur merged_samples.bam merged_samples_cur_ref.cram

bam-readcount -b 30 -w 0 -q 40 -l strict_syn_snp.txt -f $refcur merged_samples.bam \
| awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2]}' > strict_syn_snp_readcount.txt

/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/readcount_at_syn_snp_pos.py strict_syn_snp.txt strict_syn_snp_readcount.txt strict_syn_snp_allele_readcount.txt

cut -f 7 strict_syn_snp_allele_readcount.txt | sort -n | uniq -c | hometools plthist -xlim 100 800 -o strict_syn_snp_allele_readcount_total_depth.pdf
# Select positions with total read depth between 450-650
awk '{if($7>=450 && $7<=650) print $0}' strict_syn_snp_allele_readcount.txt > strict_syn_snp_allele_readcount.depth450-650.txt

cut -f 8 strict_syn_snp_allele_readcount.depth450-650.txt | sort -n | uniq -c | hometools plthist -xlim 0 1 -o strict_syn_snp_allele_readcount_maf.pdf
# Select positions with MAF between 0.4-0.6
awk '{if($8>=0.4 && $8<=0.6) print $0}' strict_syn_snp_allele_readcount.depth450-650.txt > strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt
awk '{print $1,$2,$2,$3,$5,$4,$6,$7,$8}' strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt > strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.bed

for sample in ${samples[@]}; do
    cd $cwd
    mkdir $sample
    cd $sample
#    bed=../../strict_syn_snp.selected.txt
    bed=../../strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.bed
    bsub -q bigmem -n 40 -R "span[hosts=1] rusage[mem=20000]" -M 20000 -oo rc.log -eo rc.err "
        xargs -a ${indir}/${sample}/barcodes_list \
            -P 40 \
            -I {} \
            bash -c '
                bc=\$(basename \${1})
                mkdir \$bc; cd \$bc
                /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/get_readcounts_in_cells_bt2.sh \$2 \$3 \${1}/\${bc}.DUPmarked.deduped.bam \${bc}_read_counts_b30_q40.depth450-650.af0.4-0.6.bt2.txt 40
                /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/extract_allele_count_from_bam_readcount.py \${bc}_read_counts_b30_q40.depth450-650.af0.4-0.6.bt2.txt \$2 ${sample}_\${bc}_b30_q40.depth450-650.af0.4-0.6.bt2.txt
                cd ..
            ' -- {} $bed $refcur
    "
done

# Run syri with ORA as reference
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/with_ora_ref
minimap2 -ax asm5 -t 50 --eqx ../ora.genome.v1.filtered.fasta ../cur.genome.v1.filtered.fasta \
| samtools sort -O BAM -@ 50 - \
> out.bam
samtools index -@ 50 out.bam
nohup python3 /srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri \
syri -c out.bam \
  -r ../ora.genome.v1.filtered.fasta \
  -q ../cur.genome.v1.filtered.fasta \
  -k -F B --nc 8 &

cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/
refora='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/ora.genome.v1.fasta'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")

cd $cwd
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/strict_syntenic_markers.py /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/with_ora_ref/syri.out ora_ref_strict_syn_snp.txt
awk '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}' ora_ref_strict_syn_snp.txt > ora_ref_strict_syn_snp.bed

rf ${indir}/*/barcodes/*/*.DUPmarked.deduped.ORA.bam > cell_bam_path_ora_mapped.txt
sed -i 's/^-I //g' cell_bam_path_ora_mapped.txt

samtools merge -@40 -O BAM -b cell_bam_path_ora_mapped.txt merged_samples_ora_ref.bam
samtools view -@60 -C -T $refora merged_samples_ora_ref.bam > merged_samples_ora_ref.cram
rm merged_samples_ora_ref.bam
## Can use the hometools pbamrc for parallelisation of bam-readcount
bam-readcount -b 30 -w 0 -q 40 -l ora_ref_strict_syn_snp.txt -f $refora merged_samples_ora_ref.bam \
| awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2]}' > ora_ref_strict_syn_snp_readcount.txt

/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/readcount_at_syn_snp_pos.py ora_ref_strict_syn_snp.txt ora_ref_strict_syn_snp_readcount.txt ora_ref_strict_syn_snp_allele_readcount.txt

cut -f 7 ora_ref_strict_syn_snp_allele_readcount.txt | sort -n | uniq -c | hometools plthist -xlim 100 800 -o ora_ref_strict_syn_snp_allele_readcount_total_depth.pdf

# Select positions with total read depth between 450-650
awk '{if($7>=450 && $7<=650) print $0}' ora_ref_strict_syn_snp_allele_readcount.txt > ora_ref_strict_syn_snp_allele_readcount.depth450-650.txt

cut -f 8 ora_ref_strict_syn_snp_allele_readcount.depth450-650.txt | sort -n | uniq -c | hometools plthist -xlim 0 1 -o ora_ref_strict_syn_snp_allele_readcount_maf.pdf
# Select positions with MAF between 0.4-0.6
awk '{if($8>=0.4 && $8<=0.6) print $0}' ora_ref_strict_syn_snp_allele_readcount.depth450-650.txt > ora_ref_strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt
awk '{print $1,$2,$2,$3,$5,$4,$6,$7,$8}' ora_ref_strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt > ora_ref_strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.bed

for sample in ${samples[@]}; do
    cd $cwd
    mkdir $sample
    cd $sample
    bed=../../ora_ref_strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.bed
    bsub -q bigmem -n 40 -R "span[hosts=1] rusage[mem=20000]" -M 20000 -oo rc_ora.log -eo rc_ora.err "
        xargs -a ${indir}/${sample}/barcodes_list \
            -P 40 \
            -I {} \
            bash -c '
                bc=\$(basename \${1})
                mkdir \$bc; cd \$bc
                /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/get_readcounts_in_cells_bt2.sh \$2 \$3 \${1}/\${bc}.DUPmarked.deduped.ORA.bam \${bc}_read_counts_ORA_b30_q40.depth450-650.af0.4-0.6.bt2.txt 40
                /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/extract_allele_count_from_bam_readcount.py \${bc}_read_counts_ORA_b30_q40.depth450-650.af0.4-0.6.bt2.txt \$2 ${sample}_\${bc}_ORA_b30_q40.depth450-650.af0.4-0.6.bt2.txt
                cd ..
            ' -- {} $bed $refora
    "
done


# Get aneuploidy change plots
Rscript /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/R/get_mitotic_recombination_input_data.R

# Select chromosome with low coverage variance
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/select_chromosomes_for_mitotic_recombination.py

# Run RTIGER
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/'
chrs=('CUR1G' 'CUR2G' 'CUR3G' 'CUR4G' 'CUR5G' 'CUR6G' 'CUR7G' 'CUR8G')
for chr in ${chrs[@]}; do
    cd ${cwd}/${chr}
    mkdir rtiger_co_q40_lowvar
    bsub -q bigmem -R "span[hosts=1] rusage[mem=15000]" -M 15000 -oo chr.log -eo chr.err "
        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/R/run_rtiger_on_chromosome_files.R input_q40_lowvar $chr cur rtiger_co_q40_lowvar -R 750 &
    "
done

chrs=( ORA1G ORA2G ORA3G ORA4G ORA5G ORA6G ORA7G ORA8G )
for chr in ${chrs[@]}; do
    cd ${cwd}/${chr}
    mkdir rtiger_co_q40_lowvar
    bsub -q bigmem -R "span[hosts=1] rusage[mem=15000]" -M 15000 -oo chr.log -eo chr.err "
        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/R/run_rtiger_on_chromosome_files.R input_q40_lowvar $chr ora rtiger_co_q40_lowvar -R 750 &
    "
done

## Run code here:
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_rtiger_out_and_get_stats.py

# Testing reversed CUR6G
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/CUR6G'
cd $cwd
mkdir reveresed_rtiger_co_q40
bsub -q normal -R "span[hosts=1] rusage[mem=25000]" -M 25000 -oo chr.log -eo chr.err "
    /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/R/run_rtiger_on_chromosome_files.R reversed_input_q40 CUR6G reversed_rtiger_co_q40 -R 2000 &
    "


################################################################################
########### STEP 7: Find candidate mutation for phenotypic change ##############
################################################################################

# Use the output of 4b (SNP candidates supported by multiple cells) and select SNPs that are present in both MUTs and absent in both WTs
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_candidate_for_mutant_phenotype.py


# Try finding candidate mutations by merging the data from both branches
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
# Merge the mutant samples
cd $CWD
samtools merge -r -l 9 -@ 50 -O BAM mutant.bam ../MUT_11_1/MUT_11_1.sorted.bt2.bam ../MUT_15/MUT_15.sorted.bt2.bam
samtools index -@ 60 mutant.bam
chr=$( samtools view -H  mutant.bam | grep @SQ | cut -f 2 | sed 's/SN://g' )
rm tmp.bed
for c in ${chr[@]}; do
    echo -e "${c}\t1\t100000000" >> tmp.bed
done
for c in ${chr[@]}; do
    bsub -q multicore20 -n 2 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${c}.log -eo ${c}.err "
#        samtools view -b mutant.bam $c > ${c}.bam;
#        samtools index ${c}.bam
        bam-readcount -b 30 -q 10 -w 0 -f $refcur -l tmp.bed ${c}.bam | awk '{printf \$1\" \"\$2\" \"\$3\" \"\$4; for(i=6;i<=10;i++) {n1=split(\$i,a,\":\"); printf \" \"a[2]};  for(i=11;i<=NF;i++) {n1=split(\$i,a,\":\"); printf \" \"a[1]\" \"a[2]}; printf \"\n\"}' > ${c}.rc.txt
        python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py -o ${c}.rc.filtered.txt ${c}.rc.txt
    "
done
rm mutant.filtered.txtd
for c in ${chr[@]}; do
    cat ${c}.rc.txt >> mutant.rc.txt
    cat ${c}.rc.filtered.txt >> mutant.filtered.txt
done

cut -d' ' -f 4 mutant.filtered.txt | sort -n | uniq -c | hometools plthist -o mutant.read_depth.pdf -x Depth -y Frequency -xlim 0 600
# Homo depth: 140-460
python -c "
import sys
sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/')
from get_candidate_for_mutant_phenotype_funcs import sample_snps
sample_snps(\"mutant.filtered.txt\", \"mutant.filtered\", 10, 140, 460)
"


# Allele count with samples without any filtering
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
cd $CWD
chr=$( samtools view -H  mutant.bam | grep @SQ | cut -f 2 | sed 's/SN://g' )
## Merged mutant bam
for c in ${chr[@]}; do
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${c}.log -eo ${c}.err "
#        samtools view -b mutant.bam $c > ${c}.bam;
#        samtools index ${c}.bam
        bam-readcount -b 0 -q 0 -w 0 -f $refcur -l tmp.bed ${c}.bam | awk '{printf \$1\" \"\$2\" \"\$3\" \"\$4; for(i=6;i<=10;i++) {n1=split(\$i,a,\":\"); printf \" \"a[2]};  for(i=11;i<=NF;i++) {n1=split(\$i,a,\":\"); printf \" \"a[1]\" \"a[2]}; printf \"\n\"}' > ${c}.b0.q0.rc.txt
        python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py -o ${c}.b0.q0.rc.filtered.txt ${c}.b0.q0.rc.txt
    "
done
rm mutant.b0.q0.txt mutant.b0.q0.filtered.txt
for c in ${chr[@]}; do
    cat ${c}.b0.q0.rc.txt >> mutant.b0.q0.txt
    cat ${c}.b0.q0.rc.filtered.txt >> mutant.b0.q0.filtered.txt
done

## WT bam
for s in WT_1 WT_19; do
    for c in ${chr[@]}; do
        bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${s}_${c}.log -eo ${s}_${c}.err "
            samtools view -b ../${s}/${s}.sorted.bt2.bam $c > ${s}_${c}.bam;
            samtools index ${s}_${c}.bam
            bam-readcount -b 0 -q 0 -w 0 -f $refcur -l tmp.bed ${s}_${c}.bam | awk '{printf \$1\" \"\$2\" \"\$3\" \"\$4; for(i=6;i<=10;i++) {n1=split(\$i,a,\":\"); printf \" \"a[2]};  for(i=11;i<=NF;i++) {n1=split(\$i,a,\":\"); printf \" \"a[1]\" \"a[2]}; printf \"\n\"}' > ${s}_${c}.b0.q0.rc.txt
            python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py -o ${s}_${c}.b0.q0.rc.filtered.txt ${s}_${c}.b0.q0.rc.txt
        "
    done
done
rm WT_1.b0.q0.txt WT_1.b0.q0.filtered.txt WT_19.b0.q0.txt WT_19_.b0.q0.filtered.txt
for c in ${chr[@]}; do
    cat WT_1_${c}.b0.q0.rc.txt >> WT_1.b0.q0.txt
    cat WT_1_${c}.b0.q0.rc.filtered.txt >> WT_1.b0.q0.filtered.txt
    cat WT_19_${c}.b0.q0.rc.txt >> WT_19.b0.q0.txt
    cat WT_19_${c}.b0.q0.rc.filtered.txt >> WT_19.b0.q0.filtered.txt
done


## MUT  bam
for s in MUT_11_1 MUT_15; do
    while read r; do
        p=$(echo $r | tr ' ' '_')
        bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${s}_${p}.log -eo ${s}_${p}.err "
            echo $r \
            | sed 's/ /\t/g' \
            | bam-readcount -b 0 -q 0 -w 0 -f $refcur -l /dev/stdin ../${s}/${s}.sorted.bt2.bam | awk '{printf \$1\" \"\$2\" \"\$3\" \"\$4; for(i=6;i<=10;i++) {n1=split(\$i,a,\":\"); printf \" \"a[2]};  for(i=11;i<=NF;i++) {n1=split(\$i,a,\":\"); printf \" \"a[1]\" \"a[2]}; printf \"\n\"}' > ${s}_${p}.b0.q0.rc.txt
                python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py -o ${s}_${p}.b0.q0.rc.filtered.txt ${s}_${p}.b0.q0.rc.txt
        "
    done < /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.ranges
done
cat MUT_11_1_[Cc]*b0.q0.rc.txt | sort -k 1,1 -k2,2n > MUT_11_1.b0.q0.txt
cat MUT_11_1_[Cc]*b0.q0.rc.filtered.txt | sort -k 1,1 -k2,2n > MUT_11_1.b0.q0.filtered.txt &
cat MUT_15_[Cc]*b0.q0.rc.txt | sort -k 1,1 -k2,2n > MUT_15.b0.q0.txt &
cat MUT_15_[Cc]*b0.q0.rc.txt | sort -k 1,1 -k2,2n > MUT_15.b0.q0.filtered.txt &




bam-readcount -b 30 -q 10 -w 0 -f $refcur -l tmp.bed ${c}.bam | awk '{printf $1 " " $2 " " $3 " " $4" " ; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > ${c}.rc.txt


################################################################################
################# STEP 8: CNV identification for transposons ###################
################################################################################

CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/cnv/'
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
for s in WT_1 WT_19 MUT_11_1 MUT_15 ; do
    cd $CWD; mkdir $s; cd $s
    bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
        samtools sort -@40 -n -O SAM ${INDIR}${s}/${s}.sorted.bt2.bam | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools view -S -b -@ 40 -o ${s}.bam
        # Extract the discordant paired-end alignments.
        samtools view -b -@ 40 -F 1294 ${s}.bam > ${s}.discordants.unsorted.bam

        # Extract the split-read alignments
        samtools view -h -@ 40 ${s}.bam | python2 /opt/share/software/packages/lumpy-sv-0.2.13/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb -@ 40 -o ${s}.splitters.unsorted.bam

        # Sort both alignments
        samtools sort -@ 40 ${s}.discordants.unsorted.bam -o ${s}.discordants.bam
        samtools sort -@ 40 ${s}.splitters.unsorted.bam -o ${s}.splitters.bam

        samtools view -@ 40 ${s}.bam | tail -n+100000 | python2 /opt/share/software/packages/lumpy-sv-0.2.13/scripts/pairend_distro.py  -r 150 -X 4 -N 10000 -o ${s}.lib1.histo
    "
    echo $s
done

cd $CWD/WT_1
mean=338
sd=140
s='WT_1'
bsub -q multicore40 -n 1 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
    lumpy -mw 4 -tt 0 -pe id:${s},bam_file:${s}.discordants.bam,histo_file:${s}.lib1.histo,mean:${mean},stdev:${sd},read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:${s},bam_file:${s}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > ${s}.vcf
"
cd $CWD/WT_19
mean=331
sd=137
s='WT_19'
bsub -q multicore40 -n 1 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
    lumpy -mw 4 -tt 0 -pe id:${s},bam_file:${s}.discordants.bam,histo_file:${s}.lib1.histo,mean:${mean},stdev:${sd},read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:${s},bam_file:${s}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > ${s}.vcf
"
cd $CWD/MUT_11_1
mean=340
sd=141
s='MUT_11_1'
bsub -q multicore40 -n 1 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
    lumpy -mw 4 -tt 0 -pe id:${s},bam_file:${s}.discordants.bam,histo_file:${s}.lib1.histo,mean:${mean},stdev:${sd},read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:${s},bam_file:${s}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > ${s}.vcf
"
cd $CWD/MUT_15
mean=327
sd=135
s='MUT_15'
bsub -q multicore40 -n 1 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
    lumpy -mw 4 -tt 0 -pe id:${s},bam_file:${s}.discordants.bam,histo_file:${s}.lib1.histo,mean:${mean},stdev:${sd},read_length:150,min_non_overlap:150,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:${s},bam_file:${s}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > ${s}.vcf
"

## Lumpy didn't identify any CNV. Further the identified indels were also quite noise
## ==================================================================================


# Haplotype-specific insertion and deletion regions
cat /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out \
| awk '{if($11=="INS"){print $1"\t"$2-10"\t"$3+10}}' > germline_insertions_buff10.bed
cat /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out \
| awk '{if($11=="DEL"){print $1"\t"$2-10"\t"$3+10}}' > germline_deletions_buff10.bed
cat germline_insertions_buff10.bed germline_deletions_buff10.bed | bedtools sort -i - | bedtools merge -i /dev/stdin > germline_indels_buff10.bed


# Using GATK for CNV detection
## Step 1: Recalibrate base-quality
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/cnv/
VCF=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snp.indel.vcf.gz
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
ploidy=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.gatk_ploidy.tsv
mapbed=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.genmap.E0_K51.map.unique.bed

#for s in WT_1 WT_19 MUT_11_1 MUT_15; do
#    cd ${CWD}/${s}
#    bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${s}.log -eo ${s}.err "
#        samtools view -h -@50 ${INDIR}/${s}/${s}.sorted.bt2.bam | awk '{if(\$1==\"@RG\"){print \$0\"\tPL:ILLUMINA\tSM:1\"} else {print \$0}}' | samtools view -b -@50 > ${s}.sorted.bt2.RG.bam
#    "
#done
#

# Spark variants allow parallelisation
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd ${CWD}/${s}
    bsub -q multicore40 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${s}.log -eo ${s}.err "
        export PATH=/srv/netscratch/dep_mercier/grp_schneeberger/bin/java/jdk-10.0.2/bin:$PATH
        gatk --java-options \"-Xms2G -Xmx5G -XX:ParallelGCThreads=2 -XX:ActiveProcessorCount=4 -XX:HeapBaseMinAddress=64g\" \
            BaseRecalibratorSpark \
            --spark-master local[10] \
            -I ${INDIR}/${s}/${s}.sorted.bt2.bam \
            -R $refcur \
            --known-sites $VCF \
            -O ${s}.BSQR.recal.table
         gatk --java-options \"-Xms2G -Xmx5G -XX:ParallelGCThreads=2 -XX:ActiveProcessorCount=4 -XX:HeapBaseMinAddress=64g\" \
            ApplyBQSRSpark \
            --spark-master local[10] \
            -R $refcur \
            -I ${INDIR}/${s}/${s}.sorted.bt2.bam \
            --bqsr-recal-file ${s}.BSQR.recal.table \
            -O ${s}.bsqr_recal.bam
    "
done

cd $CWD
gatk PreprocessIntervals \
    -R $refcur \
    --padding 0 \
    -imr OVERLAPPING_ONLY \
    -O cur.preprocessed.interval_list
gatk AnnotateIntervals \
    --mappability-track $mapbed \
    -L cur.preprocessed.interval_list \
    -R $refcur \
    -imr OVERLAPPING_ONLY \
    -O cur.annotated.tsv

for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd ${CWD}/${s}
    bsub -q multicore40 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${s}.log -eo ${s}.err "
        gatk CollectReadCounts \
            -L ../cur.preprocessed.interval_list \
            -R $refcur \
            -imr OVERLAPPING_ONLY \
            -I ${s}.bsqr_recal.bam \
            --format TSV \
            -O ${s}.readcounts.tsv

    "
done
gatk FilterIntervals \
    -L cur.preprocessed.interval_list \
    --annotated-intervals cur.annotated.tsv \
    -I WT_1/WT_1.readcounts.tsv \
    -I WT_19/WT_19.readcounts.tsv \
    -I MUT_11_1/MUT_11_1.readcounts.tsv \
    -I MUT_15/MUT_15.readcounts.tsv \
    -imr OVERLAPPING_ONLY \
    -O cur.cohort.gc.genmap.filtered.interval_list

gatk DetermineGermlineContigPloidy \
    -L cur.cohort.gc.genmap.filtered.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -I WT_1/WT_1.readcounts.tsv \
    -I WT_19/WT_19.readcounts.tsv \
    -I MUT_11_1/MUT_11_1.readcounts.tsv \
    -I MUT_15/MUT_15.readcounts.tsv \
    --contig-ploidy-priors $ploidy \
    --output . \
    --output-prefix ploidy \
    --verbosity DEBUG


    gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L scatter-sm/twelve_1of2.interval_list \
        -I cvg/HG00096.tsv -I cvg/HG00268.tsv -I cvg/HG00419.tsv -I cvg/HG00759.tsv \
        -I cvg/HG01051.tsv -I cvg/HG01112.tsv -I cvg/HG01500.tsv -I cvg/HG01565.tsv \
        -I cvg/HG01583.tsv -I cvg/HG01595.tsv -I cvg/HG01879.tsv -I cvg/HG02568.tsv \
        -I cvg/HG02922.tsv -I cvg/HG03006.tsv -I cvg/HG03052.tsv -I cvg/HG03642.tsv \
        -I cvg/HG03742.tsv -I cvg/NA18525.tsv -I cvg/NA18939.tsv -I cvg/NA19017.tsv \
        -I cvg/NA19625.tsv -I cvg/NA19648.tsv -I cvg/NA20502.tsv -I cvg/NA20845.tsv \
        --contig-ploidy-calls ploidy-calls \
        --annotated-intervals twelveregions.annotated.tsv \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output cohort24-twelve \
        --output-prefix cohort24-twelve_1of2 \
        --verbosity DEBUG
##__WGS parameters that increase the sensitivity of calling from Mark Walker
#    --class-coherence-length 1000.0 \
#    --cnv-coherence-length 1000.0 \
#    --enable-bias-factors false \
#    --interval-psi-scale 1.0E-6 \
#    --log-mean-bias-standard-deviation 0.01 \
#    --sample-psi-scale 1.0E-6 \
    gatk PostprocessGermlineCNVCalls \
        --model-shard-path cohort24-twelve/cohort24-twelve_1of2-model \
        --model-shard-path cohort24-twelve/cohort24-twelve_2of2-model \
        --calls-shard-path cohort24-twelve/cohort24-twelve_1of2-calls \
        --calls-shard-path cohort24-twelve/cohort24-twelve_2of2-calls \
        --allosomal-contig chrX --allosomal-contig chrY \
        --contig-ploidy-calls ploidy-calls \
        --sample-index 19 \
        --output-genotyped-intervals genotyped-intervals-cohort24-twelve-NA19017.vcf.gz \
        --output-genotyped-segments genotyped-segments-cohort24-twelve-NA19017.vcf.gz \
        --sequence-dictionary ref/Homo_sapiens_assembly38.dict

## Couldn't finish the GATK pipeline as there were constant issues. Maybe try again later.
## =======================================================================================

# Find TE insertion/deletion sites using TEPID
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/
BCLIST=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
yhindex=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.X15_01_65525S

## Merge and map scDNA reads
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD
    mkdir $s ; cd $s
    bsub -q multicore40 -n 1 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
        rm ${s}.R1.fastq.gz ${s}.R2.fastq.gz
        echo Merging fastq.gz
        cat ${BCLIST}/${s}/barcodes/*/*_dedup_R1.fastq.gz > ${s}.R1.fastq.gz
        cat ${BCLIST}/${s}/barcodes/*/*_dedup_R2.fastq.gz > ${s}.R2.fastq.gz

        echo Running tepid map
        source activate tepid
        tepid-map \
            -x $refcur \
            -y $yhindex \
            -p 40 \
            -s 400 \
            -n $s \
            -1 ${s}.R1.fastq.gz \
            -2 ${s}.R2.fastq.gz \
            -z
        source deactivate
    "
done
grep -v "#" /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.ann.gff3 \
 | egrep -v 'Low|Simple|RNA|other|Satellite'  \
 | cut -f 1,4,5,7,9 \
 | sed 's/;/\t/g ; s/ID=//g ; s/Subfamily=//g ; s/Family=//g ; s/Class=//g' > TMP.bed
paste <(cut -f 1-6 TMP.bed) <(cut -f7- TMP.bed | tr '\t' '/') \
 > cur.TE.bed


for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD
    mkdir $s ; cd $s
    bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=100000]" -M 120000 -oo ${s}.log -eo ${s}.err "
        echo Getting discordant reads
        samtools view -b -@ 40 -F 1294 ${s}.bam \
        | samtools sort -n -@ 40 -O BAM -o ${s}.discordants.bam -

        /bin/bash
        echo Running tepid discover
        source activate tepid
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/tepid/bin/python /srv/netscratch/dep_mercier/grp_schneeberger/software/TEPID/Scripts/tepid-discover \
            -k \
            -D ${s}.discordants.bam \
            -p 40 \
            -n $s \
            -c ${s}.bam \
            -s ${s}.split.bam \
            -t ../cur.TE.bed
        source deactivate
    "
done
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD; cd $s
    {
        samtools sort -@ 10 -O BAM  ${s}.split.bam > ${s}.split.sorted.bam
        samtools index -@ 10 ${s}.split.sorted.bam
        samtools sort -@ 10 -O BAM  ${s}.discordants.bam > ${s}.discordants.sorted.bam
        samtools index -@ 10 ${s}.discordants.sorted.bam
    } &
done


## Count number of Mutant specific Deletions and Insertions
### Using code in /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_candidate_for_mutant_phenotype.py

## Find TE insertion/deletion sites using TEfinder
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/tefinder
RPout=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.out.gff
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'
cd $CWD
grep -v -iE '(Motif\:[ATGC]+\-rich)|(Motif\:\([ATGC]+\)n)' $RPout > TEs.gtf
awk -F '\t' '{print $9}' TEs.gtf | awk -F '"' '{print $2}' | sort | uniq > List_of_TEs.txt

for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD;
    mkdir $s; cd $s
#    rm -r tefinder
    bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=100000]" -M 120000 -oo ${s}.log -eo ${s}.err "
        # TEfinder_mg has minor modifications to ensure proper parallelisation
        TEfinder_mg \
        -picard $picard \
        -fis 270 \
        -k 30 \
        -intermed yes \
        -maxHeapMem 2500 \
        -threads 40 \
        -workingdir tefinder \
        -alignment $INDIR/$s/${s}.sorted.bt2.bam \
        -fa $refcur\
        -gtf ../TEs.gtf \
        -te ../List_of_TEs.txt
    "
done

### Get Discordant reads BAM. This is being done manually because TEfinder jobs were crashing for this step
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD; cd $s
    {
        sams=$(ls tefinder/*/*_DiscordantPairs.bam)
        rm read_names.txt
        for sam in ${sams[@]}; do
            samtools view -@12 $sam | cut -f 1  | sort | uniq >> read_names.txt
        done
        sort read_names.txt | uniq  > read_names.sorted.uniq.txt
#        rm read_names.txt
        samtools view -H -@12 Alignments.bam > all_discordants.sam
        samtools view -@12 Alignments.bam | grep -Ff read_names.txt >> all_discordants.sam
        samtools sort -O BAM -@ 12 all_discordants.sam \
        > all_discordants.bam
        samtools index -@ 12 all_discordants.bam
    } &
done

for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD; cd $s
    bedtools sort -i tefinder/TEinsertions.bed > TEinsertions.sorted.bed
done
### Filtering done here:
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_candidate_for_mutant_phenotype.py

### Tried to find branch specific TE insertions, but could not find any

## Find TE insertion/deletion sites using SPLITREADER (V2) https://github.com/baduelp/public
## Using splitreader is too difficult as the documentation is incomplete
#refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
#CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/
#
#cd $CWD
#grep '>' TE_library.fa \
#| sed 's/>//g; s/#.*$//g' \
#| sort \
#| uniq \
#> TE_list.txt
#egrep -v 'Low|Simple|RNA|other|Satellite' TE_annotation.gff | sed 's/Subfamily/Name/g' > TE_annotation_edited.gff
#python -c "
#import sys
#sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/')
#from get_candidate_for_mutant_phenotype_funcs import generate_TEfiles_splitreader
#generate_TEfiles_splitreader('TE_annotation_edited.gff', 'TE-information-all.txt')
#"
#for s in WT_1 WT_19 MUT_11_1 MUT_15; do
#    cd ${CWD}/splitreader
#    mkdir $s ; cd $s
#    bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${s}.log -eo ${s}.err "
##        /bin/bash
##        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/SPLITREADER-beta2.5_part1_mg.sh \
##        ${s} \
##        part1 \
##        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/${s}/ \
##        .sorted.bt2 \
##        ${s} \
##        ${CWD}/splitreader/${s} \
##        ${CWD}/splitreader
##
#        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/SH/SPLITREADER-beta2.5_part2_mg.sh \
#        ${s} \
#        ${s} \
#        ${CWD}/splitreader/${s} \
#        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta \
#        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/splitreader/TE_annotation_edited.gff
#    "
#    break
#done



################################################################################
####################### STEP 10: Check telomere lenght #########################
################################################################################
# Calculate and compare the length of the telomere in the pollen/leaf/fruit WGS data. If the telomere length in the pollen and other tissues is different than that would mean that apricot loses telomeres while growing which are replenished when the gamete (pollen) are formed. If the lenght of the pollen and leaf/fruit tissue is same, then it would mean that either the tree do not lose telomere or that the telomeres are not replenished in the pollen.
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/telomer_length/
telfa=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/telomeric_sequence.fa
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta

## Align pollen reads to telomere sequence and cur genome
INDIR=/biodata/dep_mercier/grp_schneeberger/reads/Apricot/combined_scDNA/
for s in A B; do
    cd $CWD
    mkdir pollen; cd pollen
    bsub -q ioheavy -n 10 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo pollen_${s}.log -eo pollen_${s}.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 10 \
        -x $telfa \
        -1 ${INDIR}4350_and_4279_all_reads_${s}_R1_allruns_615_617_619.fastq.gz \
        -2 ${INDIR}4350_and_4279_all_reads_${s}_R2_allruns_615_617_619.fastq.gz \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools sort -@10 -O BAM - \
        > pollen_${s}_vs_telo.bam
    samtools index -@10 pollen_${s}_vs_telo.bam
    "
    bsub -q ioheavy -n 10 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo pollen_${s}.log -eo pollen_${s}.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 10 \
        -x $refcur \
        -1 ${INDIR}4350_and_4279_all_reads_${s}_R1_allruns_615_617_619.fastq.gz \
        -2 ${INDIR}4350_and_4279_all_reads_${s}_R2_allruns_615_617_619.fastq.gz \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools sort -@10 -O BAM - \
        > pollen_${s}_vs_cur.bam
    samtools index -@10 pollen_${s}_vs_cur.bam
    "
done

## Align leaf re-sequencing reads to
INDIR=/biodata/dep_mercier/grp_schneeberger/reads/Apricot/Reseq-RP.P-B.P/data/seqs/
for s in A B C D; do
    cd $CWD
    mkdir leaf; cd leaf
    bsub -q ioheavy -n 10 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo leaf_${s}.log -eo leaf_${s}.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 10 \
        -x $telfa \
        -1 ${INDIR}3649_${s}_run531_*_L008_R1_001.fastq.gz \
        -2 ${INDIR}3649_${s}_run531_*_L008_R2_001.fastq.gz \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools sort -@10 -O BAM - \
        > leaf_${s}_vs_telo.bam
    samtools index -@10 leaf_${s}_vs_telo.bam
    "
    bsub -q ioheavy -n 10 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo leaf_${s}.log -eo leaf_${s}.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 10 \
        -x $refcur \
        -1 ${INDIR}3649_${s}_run531_*_L008_R1_001.fastq.gz \
        -2 ${INDIR}3649_${s}_run531_*_L008_R2_001.fastq.gz \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools sort -@10 -O BAM - \
        > leaf_${s}_vs_cur.bam
    samtools index -@10 leaf_${s}_vs_cur.bam
    "
done

## Align scDNA-WGS reads to telomere
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD
    mkdir scdna; cd scdna
    bsub -q ioheavy -n 10 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo scdna_${s}.log -eo scdna_${s}.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 10 \
        -x $telfa \
        -1 ${INDIR}/${s}/${s}_fqfrom10xBam_bxCorrected_R1.fq.gz \
        -2 ${INDIR}/${s}/${s}_fqfrom10xBam_bxCorrected_R2.fq.gz \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools view -h -@10 -F4 - \
        | samtools sort -@10 -O BAM - \
        > scdna_${s}_vs_telo.bam
    samtools index -@10 scdna_${s}_vs_telo.bam
    "
done

## Align scDNA-WGS ONLY R2 reads to telomere and cur-genome
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd ${CWD}/scdna
    mkdir onlyR2; cd onlyR2
    reads=$(ls ${INDIR}/${s}/renamed_L50_${s}*_R2_001.fastq.gz | tr '\n' ',')
    bsub -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo scdna_${s}_telo.log -eo scdna_${s}_telo.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 20 \
        -x $telfa \
        -U $reads \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools view -h -@20 -F4 - \
        | samtools sort -@20 -O BAM - \
        > scdna_${s}_vs_telo.onlyr2.bam
    samtools index -@20 scdna_${s}_vs_telo.onlyr2.bam
    "

    bsub -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo scdna_${s}_cur.log -eo scdna_${s}_cur.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 20 \
        -x $refcur \
        -U $reads \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools view -h -@20 -F4 - \
        | samtools sort -@20 -O BAM - \
        > scdna_${s}_vs_cur.onlyr2.bam
    samtools index -@10 scdna_${s}_vs_cur.onlyr2.bam
    "
done

## Get average of read depth in the genome and in the telomere
### Pollen
INDIR=/biodata/dep_mercier/grp_schneeberger/reads/Apricot/combined_scDNA/
for s in A B; do
    cd $CWD/pollen
    bsub -q normal -n 2 -R "span[hosts=1] rusage[mem=15000]" -M 20000 -oo pollen_${s}_cov.log -eo pollen_${s}_cov.err "
    syri3.8
    hometools bamcov pollen_${s}_vs_telo.bam pollen_${s}_vs_telo.cov &
    hometools bamcov pollen_${s}_vs_cur.bam pollen_${s}_vs_cur.cov
    "
done

### Leaf
for s in A B C D; do
    cd $CWD/leaf
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 10000 -oo leaf_${s}_cov.log -eo leaf_${s}_cov.err "
    syri3.8
    hometools bamcov leaf_${s}_vs_telo.bam leaf_${s}_vs_telo.cov
    hometools bamcov leaf_${s}_vs_cur.bam leaf_${s}_vs_cur.cov
    "
done

### scdna-wgs leaf
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd ${CWD}/scdna
    bsub -q normal -n 2 -R "span[hosts=1] rusage[mem=15000]" -M 20000 -oo scdna_${s}_cov.log -eo scdna_${s}_cov.err "
    syri3.8
    hometools bamcov scdna_${s}_vs_telo.bam scdna_${s}_vs_telo.cov
    "
done

### scdna-wgs leaf onlyR2
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd ${CWD}/scdna/onlyR2
    bsub -q normal -n 2 -R "span[hosts=1] rusage[mem=15000]" -M 20000 -oo scdna_${s}_cov.log -eo scdna_${s}_cov.err "
    syri3.8
    hometools bamcov scdna_${s}_vs_telo.onlyr2.bam scdna_${s}_vs_telo.onlyr2.cov &
    hometools bamcov scdna_${s}_vs_cur.onlyr2.bam scdna_${s}_vs_cur.onlyr2.cov
    "
done


################################################################################
##### STEP 9: Check unmapped reads for mutant specific pathogen infection ######
################################################################################
# Filter unampped reads
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/patho_inf/
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD
    {
    samtools view -@10 -hf 13 -O SAM ${INDIR}/${s}/${s}.sorted.bt2.bam \
    | samtools sort -@10 -n -O SAM - \
    | samtools fasta -@10 -1 ${s}_R1.fa.gz -2 ${s}_R2.fa.gz -s /dev/null -0 /dev/null -
    } &
done

for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD
    bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 22000 -oo ${s}.log -eo ${s}.err "
#        megahit -m 20000000000 -t 20  -1 ${s}_R1.fa.gz -2 ${s}_R2.fa.gz -o ${s}_megahit_assembly
        cd-hit -i ${s}_megahit_assembly/final.contigs.fa -o ${s}_megahit_assembly/final.contigs.cdhit.fa
    "
done
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD
    cd-hit -i ${s}_megahit_assembly/final.contigs.fa -o ${s}_megahit_assembly/final.contigs.cdhit.fa > ${s}_cdhit.log 2>&1 &
done

for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD/${s}_megahit_assembly/
    for i in {00..48}; do
        bsub -q normal -n 4 -R "span[hosts=1] rusage[mem=5000]" -M 10000 -oo nt${i}.log -eo nt${i}.err "
            blastn -query final.contigs.cdhit.fa \
                -db /opt/share/blastdb/ncbi/nt.${i} \
                -out nt_${i}.oblast \
                -outfmt '6 delim=  qaccver saccver staxid pident qcovs qcovhsp qcovus length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
                -num_threads 4 \
                -perc_identity 90 \
                -subject_besthit \
                -evalue 0.000001 &
        "
    done
done

for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD/${s}_megahit_assembly/
    cat nt_{00..48}.oblast \
    | awk '{if($7>=80){print $0}}' \
    > ntall.oblast
done

################################################################################
############### STEP 9: Group cells to layers using mutations ##################
################################################################################
rf */cells_readcount/*/*readcount.bt2.txt > cells_readcount_paths.txt

<<Comment
# This part is not being used as of now, but might be important later.
rnabamdir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells/'
for sample in ${samples[@]}; do
    {
        bam-readcount -w 0 \
            -f $refcur \
            -l good_candidates.regions \
            ${rnabamdir}/${sample}/${sample}/outs/possorted_genome_bam.bam |
            awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}' \
                >${sample}/${sample}_good_candidates_rna_read_counts.txt
    } &
done

rf */cells_readcount/*/*readcount.bt2.txt >cells_readcount_paths.txt
awk '{n1=split($1, a, "/"); print $1"\t"a[11]";"a[13]}' cells_readcount_paths.txt >cells_readcount_paths_id.txt
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_variants_in_cells.py -frc cells_readcount_paths_id.txt good_candidates.regions

Comment


# <editor-fold desc="OBSELETE: LOSS OF HETEROZYGOSITY IDENTIFICATION">
## The following did not work and did not result in any good candidates
## Gene Conversion identification
#cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
#samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
#for sample in ${samples[@]}; do
#    cd ${cwd}/${sample}
#    {
#        vcftools --vcf indels/samtools/bt2/${sample}_mq1.bt2.vcf --remove-indels --recode --recode-INFO-all --out ${sample}_mq1_onlysnps.bt2
#        hometools vcfdp ${sample}_mq1_onlysnps.bt2.recode.vcf -o ${sample}_mq1_onlysnps.bt2.recode.vcf.dp
#        awk '{print $6+$7+$8+$9}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp | sort -n | uniq -c | hometools plthist -o ${sample}_mq1_onlysnps.read_depth.pdf -x read_depth -y frequency -xlim 0 300 &
#        awk '{print $8+$9}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp | sort -n | uniq -c | hometools plthist -o ${sample}_mq1_onlysnps.allele_depth.pdf -x allele_depth -y frequency -xlim 0 300 &
#        awk '{print ($8+$9)/($6+$7+$8+$9)}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp | sort -n | uniq -c | hometools plthist -o ${sample}_mq1_onlysnps.allele_freq.pdf -x allele_freq -y frequency -xlim 0 1 &
#    } &
#done
#
#for sample in ${samples[@]}; do
#    cd ${cwd}/${sample}
#    {
#        awk '{if (($6+$7+$8+$9)>=60 && ($6+$7+$8+$9)<=180) print }' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp >${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.dp
#        awk '{print $8+$9}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.dp | sort -n | uniq -c | hometools plthist -o ${sample}_mq1_onlysnps.allele_depth.dp_60_180.pdf -x allele_depth -y frequency -xlim 0 300 &
#        awk '{print ($8+$9)/($6+$7+$8+$9)}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.dp | sort -n | uniq -c | hometools plthist -o ${sample}_mq1_onlysnps.allele_freq.dp_60_180.pdf -x allele_freq -y frequency -xlim 0 1 &
#    } &
#done
#
#for sample in ${samples[@]}; do
#    cd ${cwd}/${sample}
#    {
#        awk '{if (($8+$9)/($6+$7+$8+$9)>=0.3 && ($8+$9)/($6+$7+$8+$9)<=0.6) print }' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.dp >${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.af_03_06.dp
#        awk '{print $8+$9}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.af_03_06.dp | sort -n | uniq -c | hometools plthist -o ${sample}_mq1_onlysnps.allele_depth.dp_60_180.af_03_06.pdf -x allele_depth -y frequency -xlim 0 300 &
#    } &
#done
#for sample in ${samples[@]}; do
#    cd ${cwd}/${sample}
#    {
#        awk '{if (($8+$9)>=25 && ($8+$9)<=90) print }' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.af_03_06.dp >${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.af_03_06.ad_25_90.dp
#        awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.af_03_06.ad_25_90.dp >${sample}_mq1_onlysnps.bt2.recode.vcf.dp_60_180.af_03_06.ad_25_90.bed
#    } &
#done
#
#cd $cwd
#ls */*ad_25_90.bed | xargs multiIntersectBed -header -names MUT_11_1 MUT_15 WT_19 WT_1 -i >intersect_mq1_onlysnps.bt2.recode.vcf.dp_60_180.af_03_06.ad_25_90.bed
#
#for sample in ${samples[@]}; do
#    cd ${cwd}/${sample}
#    awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' ${sample}_mq1_onlysnps.bt2.recode.vcf.dp >${sample}_mq1_onlysnps.bt2.recode.vcf.bed &
#done
#
#ls */*recode.vcf.bed | xargs multiIntersectBed -header -names MUT_11_1 MUT_15 WT_19 WT_1 -i >intersect_mq1_onlysnps.bt2.recode.vcf.bed
#
#picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'
#refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
#for sample in ${samples[@]}; do
#    {
#        cd ${cwd}/${sample}/indels/gatk_hc/bt2
#        java \
#            -jar ${picard} MergeVcfs \
#            I variants.list \
#            O ${sample}.vcf
#        gatk SelectVariants \
#            -R $refcur \
#            -V ${sample}.vcf \
#            --select-type-to-include SNP \
#            -O ${sample}.snp.vcf
#        cd ${cwd}/${sample}
#        cp indels/gatk_hc/bt2/${sample}.snp.vcf ${sample}.gatk_hc.snps.bt2.vcf
#        vcf2bed --snvs <${sample}.gatk_hc.snps.bt2.vcf >${sample}.gatk_hc.snps.bt2.bed
#    } &
#done
#
#ls */*gatk_hc.snps.bt2.bed | xargs multiIntersectBed -header -names MUT_11_1 MUT_15 WT_19 WT_1 -i >intersect_gatk_hc.snps.bt2.bed
#
#awk '{if($4==3) print}' intersect_mq1_onlysnps.bt2.recode.vcf.bed >intersect_mq1_onlysnps.bt2.recode.vcf.candidates.bed
#awk '{if($4==3) print}' intersect_gatk_hc.snps.bt2.bed >intersect_gatk_hc.snps.bt2.candidates.bed
#bedtools intersect -a intersect_mq1_onlysnps.bt2.recode.vcf.candidates.bed -b intersect_gatk_hc.snps.bt2.candidates.bed >common_candidates.bed
#awk '{print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' common_candidates.bed >common_candidates.regions
#for sample in ${samples[@]}; do
#    {
#        cd ${cwd}/${sample}
#        bedtools intersect -a ${sample}.gatk_hc.snps.bt2.bed -b ../common_candidates.bed >${sample}.gatk_hc.snps.bt2.candidates.bed
#        bam-readcount -w 0 -q 1 -f $refcur -l ../common_candidates.regions ${sample}.sorted.RG.bt2.bam |
#            awk '{n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}' \
#                >common_candidates.read_count.txt
#    } &
#done
# </editor-fold>


####################################################################
############ Step 5: Indel positions to be filtered out
####################################################################

## Calling indels using Manta, samtools, and GATK HaplotypeCaller

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'

# Manta
manta_config='/srv/netscratch/dep_mercier/grp_schneeberger/software/manta-1.6.0.centos6_x86_64/bin/configManta.py'
python='/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python'

cd $cwd
for sample in ${samples[@]}; do
    cd $cwd/$sample
    mkdir indels
    cd indels
    mkdir manta
    cd manta
    echo "IN $sample"
    for sample2 in ${samples[@]}; do
        if [[ $sample == $sample2 ]]; then
            continue
        else

            ## Submit BT2 based jobs
            bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${sample}_${sample2}_bt2.log -eo ${sample}_${sample2}_bt2.err "
        $python \
          $manta_config \
            --tumourBam ${cwd}/${sample}/${sample}.sorted.RG.bt2.bam \
            --normalBam ${cwd}/${sample2}/${sample2}.sorted.RG.bt2.bam \
            --referenceFasta $refcur \
            --runDir ${sample}_vs_${sample2} \
            --scanSizeMb 2
        cd ${sample}_vs_${sample2}
        $python \
          ./runWorkflow.py \
            -m local \
            -j 40 \
            -g 50
        gunzip -c results/variants/candidateSmallIndels.vcf.gz > candidateSmallIndels.vcf
        vcf2bed --do-not-sort --insertions <candidateSmallIndels.vcf> candidateSmallIndels.insertions.bed
        vcf2bed --do-not-sort --deletions <candidateSmallIndels.vcf> candidateSmallIndels.deletions.bed
      "
        fi
    done
done

# Samtools
cur_contig_list='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.contig_list'
ploidy='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/marker_generated/ploidy'
cd $cwd
for sample in ${samples[@]}; do
    cd $cwd/$sample
    cd indels
    mkdir samtools
    cd samtools

    # Run job for bt2 alignments
    bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=30000]" -M 36000 -o samtools_bt2.log -e samtools_bt2.err "
    cat $cur_contig_list \
    | sed 's/\\n/ /g' \
    | xargs -d '\n' -n 1 -P 40 -I {} /bin/bash -c ' \
        bcftools mpileup -q \$2 -d 800 --threads 1 -f \$3 -r \$1 \${4} > \${1}_mq\${2}.bt2.pileup \
      ' -- {} 1 $refcur ${cwd}/${sample}/${sample}.sorted.RG.bt2.bam

    bcftools concat --threads 38 *_mq1.bt2.pileup \
    | bcftools call --threads 38 -m -v -Ov --ploidy-file ${ploidy} \
    > ${sample}_mq1.bt2.vcf
    vcftools --vcf ${sample}_mq1.bt2.vcf --keep-only-indels --recode --recode-INFO-all --out ${sample}_mq1_onlyindels.bt2.vcf

    vcf2bed --do-not-sort --insertions <${sample}_mq1_onlyindels.bt2.vcf.recode.vcf> ${sample}_mq1_onlyindels.bt2.vcf.recode.insertions.bed
    vcf2bed --do-not-sort --deletions <${sample}_mq1_onlyindels.bt2.vcf.recode.vcf> ${sample}_mq1_onlyindels.bt2.vcf.recode.deletions.bed
    "
done

# GATK HaplotypeCaller
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.formatted.ranges'
jv="--java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=2'"
cd $cwd
for sample in ${samples[@]}; do
    cd $cwd/$sample
    cd indels
    mkdir gatk_hc
    cd gatk_hc

    # Step 1
#    while read r; do
#        r2=$( echo $r | sed 's/:/_/g' | sed 's/-/_/g'  )
#        bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}.log -eo ${r2}.err "
#            gatk $jv \
#            HaplotypeCaller \
#            -R $refcur \
#            -I ${cwd}/${sample}/${sample}.sorted.bt2.bam \
#            -O ${r2}.vcf.gz \
#            --native-pair-hmm-threads 1 \
#            -L $r
#        "
#    done < $regions


#     #Step 2
    mq_value='35.00'
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo indel.log -eo indel.err "
        ls CUR*.vcf.gz > variants.list
        java -XX:ParallelGCThreads=1 \
            -jar ${picard} MergeVcfs \
            -I variants.list \
            -O TMP.vcf
        gatk $jv SelectVariants \
           -R $refcur \
           -V TMP.vcf \
           --select-type-to-include INDEL \
           -O TMP.indel.vcf
        gatk $jv VariantFiltration \
            -R $refcur \
            -V TMP.indel.vcf \
            -O TMP.indel.filter.vcf \
            --filter-expression \"QUAL < 0 || QD < 2.00 || FS > 60.000 || SOR > 3.000  || MQ < $mq_value || MQRankSum < -10.00 || ReadPosRankSum < -6.000 || ReadPosRankSum > 4.000\" \
            --filter-name \"indel_filter\" \
            --filter-expression \"DP < 50 || DP > 300\" \
            --filter-name \"DP_filter\"
        grep -E '^#|PASS' TMP.indel.filter.vcf > indels.vcf
        vcf2bed --do-not-sort --insertions <indels.vcf> insertions.bed
        vcf2bed --do-not-sort --deletions <indels.vcf> deletions.bed
        vcf2bed --do-not-sort <TMP.indel.vcf > indels.unfilt.bed
#        rm TMP.*
   "
     cd ..
done
# Comparisons of indels predicted by GATK also shows that there are no good MUTATION SPECIFIC varations

while read r; do
    r2=$( echo $r | sed 's/:/_/g' | sed 's/-/_/g'  )
    {
        gatk --java-options '-Xms5G -Xmx20G -XX:ParallelGCThreads=5' \
        HaplotypeCaller \
        -R $refcur \
        -I ${cwd}/${sample}/${sample}.sorted.bt2.bam \
        -O ${r2}.vcf.gz \
        --native-pair-hmm-threads 1 \
        -L $r 2> ${r2}.err
    } &
done < TMP.regions

jv="--java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=2'"
for sample in WT_1 MUT_11_1 MUT_15; do
for sample in MUT_11_1 ; do
    cd $cwd/$sample/indels/gatk_hc
    # Step 2
    mq_value='35.00'
    ls CUR*.vcf.gz > variants.list
    java -XX:ParallelGCThreads=1 \
        -jar ${picard} MergeVcfs \
        -I variants.list \
        -O TMP.vcf
    gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=2' SelectVariants \
       -R $refcur \
       -V TMP.vcf \
       --select-type-to-include INDEL \
       -O TMP.indel.vcf
    gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=2' VariantFiltration \
        -R $refcur \
        -V TMP.indel.vcf \
        -O TMP.indel.filter.vcf \
        --filter-expression 'QUAL < 0 || QD < 2.00 || FS > 60.000 || SOR > 3.000  || MQ < $mq_value || MQRankSum < -10.00 || ReadPosRankSum < -6.000 || ReadPosRankSum > 4.000' \
        --filter-name "indel_filter" \
        --filter-expression "DP < 50 || DP > 300" \
        --filter-name "DP_filter"
    grep -E '^#|PASS' TMP.indel.filter.vcf > indels.vcf
    vcf2bed --do-not-sort --insertions <indels.vcf> insertions.bed
    vcf2bed --do-not-sort --deletions <indels.vcf> deletions.bed
    vcf2bed --do-not-sort <TMP.indel.vcf > indels.unfilt.bed
     cd ..
done

while read r; do
    r2=$( echo $r | sed 's/:/_/g' | sed 's/-/_/g'  )
    gunzip -d -c ${r2}.vcf.gz | awk 'OFS="\t" {if($1=="#CHROM"){$10="WT_1"; print $0} else {print $0}}' | bgzip > TMP_${r2}.vcf.gz
    mv TMP_${r2}.vcf.gz ${r2}.vcf.gz
done < TMP.regions



# Merge all the called indel calls. This will give us all regions where identified
# sSNVs would be of low confidence because of alignment errors.
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
for sample in ${samples[@]}; do
    echo $sample
    cd $cwd/$sample/indels
    rm insertions.list deletions.list

    ls manta/*/*insertions.bed >>insertions.list
    ls samtools/*/*insertions.bed >>insertions.list
    ls gatk_hc/*/*insertions.bed >>insertions.list
    while read r; do
        cat $r >>insertions.all.bed
    done <insertions.list
    cut -f1,2,3 insertions.all.bed >insertions.bed
    sortBed -faidx /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.fai -i insertions.bed >insertions.sorted.bed
    mergeBed -i insertions.sorted.bed >insertions.sorted.merged.bed

    ls manta/*/*deletions.bed >>deletions.list
    ls samtools/*/*deletions.bed >>deletions.list
    ls gatk_hc/*/*deletions.bed >>deletions.list
    while read r; do
        cat $r >>deletions.all.bed
    done <deletions.list
    cut -f1,2,3 deletions.all.bed >deletions.bed
    sortBed -faidx /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.fai -i deletions.bed >deletions.sorted.bed
    mergeBed -i deletions.sorted.bed >deletions.sorted.merged.bed

    ## Add buffer of 20 BP around indel positions
    cat insertions.sorted.merged.bed |
        awk '{printf $1"\t"} {if($2<20) {printf 0"\t"} else {printf $2-20"\t"}} {print $3+20}' \
            >insertions.sorted.merged.padded_20.bed
    cat deletions.sorted.merged.bed |
        awk '{printf $1"\t"} {if($2<20) {printf 0"\t"} else {printf $2-20"\t"}} {print $3+20}' \
            >deletions.sorted.merged.padded_20.bed

    rm *.all.bed *.sorted.bed insertions.bed deletions.bed
done

## Step 4c: Select candidates that are in highly unqiue (low mappability) regions

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
unimapbed='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.genmap.E0_K51.map.unique.bedgraph'
repeatbed='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/repeat_region/cur.genome.v1.fasta.out.padded_20.bed'
# inbed='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.genmap.E0_K51.map.unique.bedgraph'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
cd $cwd
## Use python3.7 environment
for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
    echo $sample

    FOR MM2 SNPs
    Select candidates that are in uniquely mapping regions
    cat multi_cell_${sample}_filtered_SNPs_candidate.sorted.bed |
        awk '{if($2>=25 && $3>=26){print $1"\t"$2-25"\t"$3-25"\t"$4"\t"$5"\t"$6"\t"$7}}' \
            >multi_cell_${sample}_filtered_SNPs_candidate.sorted.pos_adjusted.bed

    bedtools intersect \
        -a multi_cell_${sample}_filtered_SNPs_candidate.sorted.pos_adjusted.bed \
        -b $unimapbed |
        awk '{print $1"\t"$2+25"\t"$3+25"\t"$4"\t"$5"\t"$6"\t"$7}' \
            >multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.bed

    # Remove candidates near repeatitive regions
    bedtools intersect \
        -v \
        -a multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.bed \
        -b $repeatbed \
        >multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.bed

    # Remove candidates near indels
    bedtools intersect \
        -v \
        -a multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.bed \
        -b indels/deletions.sorted.merged.padded_20.bed \
        >no_del.bed
    bedtools intersect \
        -v \
        -a no_del.bed \
        -b indels/insertions.sorted.merged.padded_20.bed \
        >multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed

    # FOR BT2 SNPs
    # Select candidates that are in uniquely mapping regions
    cat multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.bed |
        awk '{if($2>=25 && $3>=26){print $1"\t"$2-25"\t"$3-25"\t"$4"\t"$5"\t"$6"\t"$7}}' \
            >multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.pos_adjusted.bed

    bedtools intersect \
        -a multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.pos_adjusted.bed \
        -b $unimapbed |
        awk '{print $1"\t"$2+25"\t"$3+25"\t"$4"\t"$5"\t"$6"\t"$7}' \
            >multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.bed

    # Remove candidates near repeatitive regions
    bedtools intersect \
        -v \
        -a multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.bed \
        -b $repeatbed \
        >multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.bed

    # Remove candidates near indels
    bedtools intersect \
        -v \
        -a multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.bed \
        -b indels/deletions.sorted.merged.padded_20.bed \
        >no_del_bt2.bed
    bedtools intersect \
        -v \
        -a no_del_bt2.bed \
        -b indels/insertions.sorted.merged.padded_20.bed \
        >multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed

    cut -f6 multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed |
        sort -n |
        uniq -c |
        hometools plthist -o ${sample}_candidates_allele_counts.pdf -W 6 -H 3 -x alternate allele count -y number of candidates -ylog -t ${sample}_mm2 &

    cut -f6 multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed |
        sort -n |
        uniq -c |
        hometools plthist -o ${sample}_bt2_candidates_allele_counts.pdf -W 6 -H 3 -x alternate allele count -y number of candidates -ylog -t ${sample}_bt2 &

    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_mutation_changes_histogram.py \
        -f multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed \
        multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed \
        -rc 4 \
        -qc 5 \
        -samples ${sample}_MM2 ${sample}_BT2 \
        -t Mutation changes in $sample \
        -W 8 \
        -ymax 40 \
        -o mutation_changes_${sample}.png &
done

# Step 4d: Perform read-level filtering

# run using python 3.7
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
for sample in ${samples[@]}; do
    cd ${cwd}/${sample}
    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_candidates_based_on_read_alignment_features.py \
        multi_cell_${sample}_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed \
        ${sample}.sorted.RG.bam \
        $refcur \
        -o ${sample} &

    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_candidates_based_on_read_alignment_features.py \
        multi_cell_bt2_${sample}_bt2_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.non_indel.bed \
        ${sample}.sorted.RG.bt2.bam \
        $refcur \
        -o ${sample}_bt2 &
done


################################################################################
################################################################################
################################################################################
################################################################################

## STEP 5a : STRELKA BASED sSNVs
# Calling variants using Strelka: For each sample, call
# variant in that sample (considering it as tumor) compared
# to all other samples (considered as normals). Variants
# identified in the sample in all comparisons are selected
# for further processing.

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/strelka/'
stlka_config='/srv/netscratch/dep_mercier/grp_schneeberger/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")

## Run strelka on all combinations of samples
for sample in ${samples[@]}; do
    cd $cwd
    mkdir $sample
    cd $sample
    echo "IN $sample"
    for sample2 in ${samples[@]}; do
        if [[ $sample == $sample2 ]]; then
            continue
        else
            ## Submit BT2 based jobs
            bsub -q multicore40 -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${sample}_${sample2}_bt2.log -eo ${sample}_${sample2}_bt2.err "
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python \
          $stlka_config \
            --tumourBam ${indir}/${sample}/${sample}.sorted.bt2.bam \
            --normalBam ${indir}/${sample2}/${sample2}.sorted.bt2.bam \
            --referenceFasta $refcur \
            --runDir ${sample}_vs_${sample2} \
            --scanSizeMb 2
        cd ${sample}_vs_${sample2}
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python \
          ./runWorkflow.py \
            -m local \
            -j 40 \
            -g 50
      "
        fi
    done
done

## STEP 5b : Mutect2 BASED sSNVs
# First call mutations in all individual samples, RojoPassion,
# and the parents data in tumor-only mode. Then for each sample
# create a Panel_of_Normals using all vcf files (other than the
# focal sample) identified above. Finally, call variants in the
# comparative mode using sample data (from cell reads) and the
# corresponding Panel_of_Normals.

# Add RG tag in the bam files for the Currot and OrangeRed samples
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/'
samples=('currot' 'orangered')
for sample in ${samples[@]}; do
    cd ${cwd}/$sample
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 8000 -oo RG.log -eo RG.err "
      java -XX:ParallelGCThreads=1 -jar $picard AddOrReplaceReadGroups \
        I=${sample}.dedup.sorted.bam \
        O=${sample}.dedup.sorted.RG.bam \
        LB=${sample} PL=illumina PU=unit1 SM=${sample}
      samtools index ${sample}.dedup.sorted.RG.bam
    "
done

# Add RG tag in the bam files for Rojo Passion merged samples also submit jobs for calling variants
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/rp_leaf_illumina/merged_samples/'
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'
cd $cwd
# Step 1
bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 8000 -oo add_RG.log -eo add_RG.err "
  # Merge samples to get representative tree sequence
  samtools merge -O BAM -@ 80 merged.bam ../*/*.DUPmarked.deduped.bam
  samtools index -@60 merged.bam

  java -XX:ParallelGCThreads=1 -jar $picard AddOrReplaceReadGroups \
    I=merged.bam \
    O=merged.RG.bam \
    LB=RP PL=illumina PU=unit1 SM=RP
  samtools index merged.RG.bam
  while read r; do
    r2=\$( echo \$r | sed 's/:/_/g' | sed 's/-/_/g'  )
    bsub -q multicore20 -n 1 -R \"span[hosts=1] rusage[mem=5000]\" -M 5000 -oo \${r2}_mm2.log -eo \${r2}_mm2.err \"
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
        Mutect2 \
          -R $refcur \
          -I merged.RG.bam \
          -tumor RP \
          -O \${r2}.vcf.gz \
          --native-pair-hmm-threads 1 \
          -L \$r
    \"
  done < $regions
"
# Step 2
bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo sort.log -eo sort.err "
  bcftools concat -O z CUR*.vcf.gz > all_variants.vcf.gz
  java -jar $picard SortVcf I=all_variants.vcf.gz O=all_variants.sorted.vcf.gz
"

#---------------- DETECT VARIANTS ---------------------------------
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/mutect2/'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'
# Get variations in tumor-only mode for scdna-seq samples
cd $cwd
for sample in ${samples[@]}; do
    cd $cwd
    mkdir $sample
    cd $sample
    mkdir tumor_only_mode
    cd tumor_only_mode

    # Step 1
    while read r; do
        r2=$(echo $r | sed 's/:/_/g' | sed 's/-/_/g')
        bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}_mm2.log -eo ${r2}_mm2.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
        Mutect2 \
        -R $refcur \
        -I ${indir}/${sample}/${sample}.sorted.RG.bam \
        -tumor $sample \
        -O ${r2}_mm2.vcf.gz \
        --native-pair-hmm-threads 1 \
        -L $r
    "

        bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}_bt2.log -eo ${r2}_bt2.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
        Mutect2 \
        -R $refcur \
        -I ${indir}/${sample}/${sample}.sorted.RG.bt2.bam \
        -tumor $sample \
        -O ${r2}_bt2.vcf.gz \
        --native-pair-hmm-threads 1 \
        -L $r
    "
    done <$regions

    # Step 2
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo sort.log -eo sort.err "
    bcftools concat -O z CUR*_mm2.vcf.gz > all_variants.mm2.vcf.gz
    java -jar $picard SortVcf I=all_variants.mm2.vcf.gz O=all_variants.mm2.sorted.vcf.gz

    bcftools concat -O z CUR*_bt2.vcf.gz > all_variants.bt2.vcf.gz
    java -jar $picard SortVcf I=all_variants.bt2.vcf.gz O=all_variants.bt2.sorted.vcf.gz
  "
done

# Get variations in tumor-only mode for the Currot and OrangeRed samples
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/'
samples=('currot' 'orangered')
for sample in ${samples[@]}; do
    cd ${cwd}/$sample
    mkdir mutect2_var_call
    cd mutect2_var_call

    # Step 1
    while read r; do
        r2=$(echo $r | sed 's/:/_/g' | sed 's/-/_/g')
        bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}_mm2.log -eo ${r2}_mm2.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
        Mutect2 \
          -R $refcur \
          -I ${cwd}/${sample}/${sample}_cur_v1.1.dedup.sorted.RG.bam \
          -tumor $sample \
          -O ${sample}.${r2}.vcf.gz \
          --native-pair-hmm-threads 1 \
          -L $r
    "
    done <$regions

    # Step 2
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo sort.log -eo sort.err "
    bcftools concat -O z ${sample}.CUR*.vcf.gz > all_variants.vcf.gz
    java -jar $picard SortVcf I=all_variants.vcf.gz O=all_variants.sorted.vcf.gz
  "
done

# Call variations using Mutect2
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/mutect2/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'
samples=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
cd $cwd

rf */tumor_only_mode/all_variants.mm2.sorted.vcf.gz >variants_mm2.args
rf */tumor_only_mode/all_variants.bt2.sorted.vcf.gz >variants_bt2.args
echo /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/rp_leaf_illumina/merged_samples/all_variants.sorted.vcf.gz | tee -a variants_mm2.args variants_bt2.args
echo /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/mutect2_var_call/all_variants.sorted.vcf.gz | tee -a variants_mm2.args variants_bt2.args
echo /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/mutect2_var_call/all_variants.sorted.vcf.gz | tee -a variants_mm2.args variants_bt2.args

# Create Panel Of Normals
for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
    for dir in mm2 bt2; do
        mkdir $dir
        cd $dir
        bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo pon.log -eo pon.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=1' \
        CreateSomaticPanelOfNormals \
        -vcfs ${cwd}/variants_${dir}.args \
        -O pon_${dir}.vcf.gz
    "
        cd ..
    done
done

# Get variants
for sample in ${samples[@]}; do
    cd $cwd
    cd $sample
    for dir in mm2 bt2; do
        echo $sample $dir
        tool=''
        if [[ $dir == 'bt2' ]]; then
            tool='.bt2'
        fi
        cd $dir
        for sample2 in ${samples[@]}; do
            if [[ $sample == $sample2 ]]; then
                continue
            else

                mkdir ${sample}_vs_${sample2}
                cd ${sample}_vs_${sample2}
                while read r; do
                    if [[ $r != "CUR"* ]]; then
                        continue
                    else
                        r2=$(echo $r | sed 's/:/_/g' | sed 's/-/_/g')
                        bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}.log -eo ${r2}.err "
              gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=1' \
                Mutect2 \
                  -R $refcur \
                  -I ${indir}/${sample}/${sample}.sorted.RG${tool}.bam \
                  -tumor $sample \
                  -I ${indir}/${sample2}/${sample2}.sorted.RG${tool}.bam \
                  -normal ${sample2} \
                  -O ${r2}.vcf.gz \
                  --panel-of-normals ${cwd}/${sample}/${dir}/pon_${dir}.vcf.gz \
                  --native-pair-hmm-threads 1 \
                  -L $r
            "
                    fi
                done <$regions
                cd ..
            fi
        done
        cd ..
    done
done

#######################################################################################
######################## END OF NECESSARY COMMANDS ####################################
#######################################################################################
## Testing commands
cd $cwd
samtools depth -a -q 30 -Q 10 -d 0 WT_1/WT_1.sorted.bam WT_19/WT_19.sorted.bam MUT_11_1/MUT_11_1.sorted.bam MUT_15/MUT_15.sorted.bam >read_depth_q30_Q10.txt

#-----------------------------------------------------------------------
## Test different sSNP identification methods using simulated BAM files
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/CUR1G.fasta'
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/CUR1G'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/'
cd $cwd

bowtie2 --end-to-end \
    --very-sensitive \
    --threads 40 \
    -x $curidxbt2 \
    -1 WT_1_R1.fastq.gz \
    -2 WT_1_R2.fastq.gz \
    -S CUR1G.bt2.sam

samtools sort -O BAM CUR1G.bt2.sam >CUR1G.bt2.sorted.bam
samtools index CUR1G.bt2.sorted.bam
samtools view -b -f 2 CUR1G.bt2.sorted.bam >CUR1G.bt2.sorted.onlymapped.bam

# Add RG information
java -jar $picard AddOrReplaceReadGroups I=CUR1G.bt2.sorted.onlymapped.bam O=CUR1G.bt2.sorted.onlymapped.RG.bam LB=WT_1 PL=illumina PU=unit1 SM=1

# Check validity of BAM file
java -jar $picard ValidateSamFile I=CUR1G.bt2.sorted.onlymapped.RG.bam O=validation.report M=VERBOSE

samtools index CUR1G.bt2.sorted.onlymapped.RG.bam

# Get all positions with no alternate base and read-coverage between 100-200 bps, SNPs would be added within these positions only
bam-readcount -w 0 -f $refcur CUR1G.bt2.sorted.onlymapped.RG.bam |
    awk '{if($4>=100 && $4<=200) {n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' \
        >bam_read_counts_b30_q20.bt2.txt

python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/get_conserved_positions.py bam_read_counts_b30_q20.bt2.txt

addsnv.py -v mutate_loci_bam_read_counts_b30_q20.bt2.txt \
    -f CUR1G.bt2.sorted.onlymapped.RG.bam \
    -r ../CUR1G.fasta \
    -o CUR1G.bt2.sorted.onlymapped.RG.spiked.bam \
    -p 60 \
    --picardjar /srv/netscratch/dep_mercier/grp_schneeberger/software/picard_2.25.0/picard.jar \
    --mindepth 100 \
    --maxdepth 200 \
    --tagreads \
    --aligner bowtie2 \
    --alignopts end-to-end:,very-sensitive:,threads:40,bowtie2ref:CUR1G \
    --seed 1

cd WT_1
samtools sort -O BAM -@ 50 CUR1G.bt2.sorted.onlymapped.RG.spiked.bam >CUR1G.spiked.sorted.bam
samtools index CUR1G.spiked.sorted.bam

## ABOVE STEPS ARE REPEATED FOR WT_19 SAMPLE AS WELL

# Using py2.6mg environment run STRELKA
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/'
cd $cwd
python /srv/netscratch/dep_mercier/grp_schneeberger/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam WT_19/CUR1G.bt2.sorted.onlymapped.RG.bam \
    --tumourBam WT_1/CUR1G.spiked.sorted.bam \
    --referenceFasta CUR1G.fasta \
    --runDir ./mut_in_wt_1 \
    --scanSizeMb 2
cd mut_in_wt_1
python ./runWorkflow.py \
    -m local \
    -j 60 \
    -g 250

# run strelka for simulated WT_19
python /srv/netscratch/dep_mercier/grp_schneeberger/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --tumourBam CUR1G.spiked.sorted.bam --normalBam ../WT_1/CUR1G.bt2.sorted.onlymapped.RG.bam --referenceFasta ../CUR1G.fasta --runDir ./mut_in_wt_19 --scanSizeMb 2
python ./runWorkflow.py \
    -m local \
    -j 60 \
    -g 250

# Running mutect2 on WT_1
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_1
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
gatk Mutect2 \
    -R $refcur \
    -I CUR1G.spiked.sorted.RG.bam \
    -I ../WT_19/CUR1G.spiked.sorted.RG.bam \
    -tumor-sample WT_1 \
    -normal WT_19 \
    --native-pair-hmm-threads 40 \
    -O WT_1_vs_WT_19.vcf.gz

# Running mutect2 on WT_19
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_19
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
gatk Mutect2 \
    -R $refcur \
    -I CUR1G.spiked.sorted.RG.bam \
    -I ../WT_1/CUR1G.spiked.sorted.RG.bam \
    -tumor-sample WT_19 \
    -normal WT_1 \
    --native-pair-hmm-threads 40 \
    -O WT_19_vs_WT_1.vcf.gz
#-----------------------------------------------------------------------