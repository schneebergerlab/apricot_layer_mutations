#####################################################################
### Testing directly using assemblers
#####################################################################
reads='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_hifi/4784_A_run521_HIFI.fastq'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/kmer_analysis'
cd $cwd

bsub -q multicore40 -n 40 -R "rusage[mem=50000] span[hosts=1]" -M 50000 -oo jf.log -eo jf.err "
  jellyfish count -C -o hifi -m 21 -t 40 -s 5G $reads;
  jellyfish histo -h 250000 -o hifi.histo hifi;
"

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/test/'
reads='/srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4784/4784_A_run521_HIFI.fastq.gz'
cd $cwd

########## RUN HIFIASM ############################
hifiasm='/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/hifiasm'
bsub -q multicore40 -n 40 -R "rusage[mem=200000] span[hosts=1]" -M 200000 -oo s1.log -eo s1.err "
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    $hifiasm -t 40 -o hifiasm $reads
"
cat hifiasm.p_utg.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' >hifiasm.p_utg.fasta

######### CANU IN PACBIO-HIFI MODE
# Ran on Dell-node
cd $cwd
canu -p canu \
  -d canu \
  genomeSize=240m \
  useGrid=false \
  -pacbio-hifi \
  $reads \
  2>canu.log

######### CANU WITH TRIO-BINNING
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/test/canu_trio'
reads='/srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4784/4784_A_run521_HIFI.fastq.gz'
curR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/3710_A_run527_GAGGATGG_S175_L005_R1_001.fastq.gz'
curR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/3710_A_run527_GAGGATGG_S175_L005_R2_001.fastq.gz'
oraR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/3710_B_run527_GATTCATC_S176_L005_R1_001.fastq.gz'
oraR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/3710_B_run527_GATTCATC_S176_L005_R2_001.fastq.gz'
cd $cwd
canu -haplotype \
  -p canu_trio \
  -d canu_trio \
  genomeSize=240m \
  -haplotypecur $curR1 $curR2 \
  -haplotypeora $oraR1 $oraR2 \
  useGrid=false \
  -pacbio-raw \
  $reads \
  2>canu.log

#####################################################################
### Step 1: Phase reads using gamete-binning and trio-binning
#####################################################################

## Gamete Binning
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/'
refgenome='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/campoy_sun_genome_biology_2020/manually_curated.fasta'
reads='/srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4784/4784_A_run521_HIFI.fastq.gz'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/campoy_sun_genome_biology_2020/'
long_read_genotyper='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/campoy_sun_genome_biology_2020/long_read_genotyper'

cd $cwd
bsub -q multicore40 -n 40 -R "rusage[mem=50000] span[hosts=1]" -M 60000 -oo s1.log -eo s1.err "
#  minimap2 -ax asm20 -t 40 $refgenome $reads \
#  | samtools sort -@ 40 -O BAM - \
#  > apricot_hifi.sorted.bam
#  samtools index apricot_hifi.sorted.bam
#  samtools view -O BAM -@ 50  -F 256 -F 2048 apricot_hifi.sorted.bam > apricot_hifi_F256_F2048.sorted.bam
#  samtools index apricot_hifi_F256_F2048.sorted.bam
  $long_read_genotyper --sam apricot_hifi_F256_F2048.sorted.sam \
    --marker ${indir}s4_phased_markers.txt \
    --marker2 ${indir}phased_s2_genotype_contig_seq_del_like.txt \
    --phase ${indir}zphase_contigs_linksage.txt \
    --ims 0.9 \
    -o apricot_hifi \
    > apricot_hifi_phasing.log
"

## Trio-binning
# Use reads phased using canu -haplotype mode in the tests above
# Symbolic links are created in the /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/trio_binning

## Assemble gamete-binning based linkage group -------------------------------------------
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/lg_assemblies/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads/'
cd $cwd
for i in {1..8}; do
  for s in PPP MMM; do
    mkdir ${i}_${s}_assembly
    cd ${i}_${s}_assembly
    bsub -q multicore40 -n 40 -R "rusage[mem=100000] span[hosts=1]" -M 120000 -oo s1.log -eo s1.err "
      hifiasm -t 40 -r 4 -o ${i}_${s} ${indir}/${i}.txt_${s}_pbreads.fa
    "
    cd ..
  done
done

## K-mer haplotyping accuracy measurement using meryl ------------------------------------------------

# Create meryl database for the parental reads
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/'
adapter='/srv/netscratch/dep_mercier/grp_schneeberger/projects/hyper_co/data/reads/shqMergedAdapters_Primers_representative_rc.fa'
meryl='/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_jose/meryl/build/bin/meryl'
samples=('currot' 'orangered')
for sample in ${samples[@]}; do
  cd ${indir}/${sample}
  bsub -q multicore40 -n 40 -R "rusage[mem=50000] span[hosts=1]" -M 60000 -oo s1.log -eo s1.err "
    # Trim reads
    skewer -r 0.1 -d 0.05 -k 8 -q 20 -l 75 -m pe \
      -t 40 \
      -x $adapter \
      -z -o ${sample}_ql \
      *_L005_R1_001.fastq.gz \
      *_L005_R2_001.fastq.gz
    $meryl k=21 count threads=40 memory=150 output ${sample}_R1.meryl ${sample}_ql-trimmed-pair1.fastq.gz
    $meryl k=21 count threads=40 memory=150 output ${sample}_R2.meryl ${sample}_ql-trimmed-pair2.fastq.gz
    $meryl union-sum output ${sample}.meryl ${sample}_R1.meryl ${sample}_R2.meryl
    "
done

# Create meryl data for the phased reads
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads'
cd $cwd
for i in {1..8}; do
  for s in PPP MMM; do
    bsub -q multicore40 -n 10 -R "rusage[mem=50000] span[hosts=1]" -M 60000 -oo ${i}_${s}.log -eo ${i}_${s}.err "
      $meryl k=21 count threads=10 memory=50 output ${i}_${s}.meryl ${i}.txt_${s}_pbreads.fa
    "
  done
done

# Create meryl data for the phased reads
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/trio_binning'
cd $cwd
for s in cur ora; do
  bsub -q multicore40 -n 10 -R "rusage[mem=50000] span[hosts=1]" -M 60000 -oo ${s}.log -eo ${s}.err "
    $meryl k=21 count threads=10 memory=50 output ${s}.meryl haplotype-${s}.fasta.gz
  "
done

## Generate hapmers
# For gamete-binning phased reads
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads'
oramer='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered.meryl'
curmer='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot.meryl'
export MERQURY='/srv/netscratch/dep_mercier/grp_schneeberger/software/merqury-1.1'
cd $cwd
for i in {1..8}; do
  for s in PPP MMM; do
    mkdir ${i}_${s}_hapmers
    cd ${i}_${s}_hapmers
    bsub -q normal -n 1 -R "rusage[mem=5000] span[hosts=1]" -M 6000 -oo ${i}_${s}.log -eo ${i}_${s}.err "
      meryl histogram ../${i}_${s}.meryl > ${i}_${s}.hist
      export MERQURY='/srv/netscratch/dep_mercier/grp_schneeberger/software/merqury-1.1'
      $MERQURY/trio/hapmers.sh $curmer $oramer ../${i}_${s}.meryl
    "
    cd ..
  done
done

# For trio-binning phased reads
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/trio_binning'
cd $cwd
for s in cur ora; do
  mkdir ${s}_hapmers
  cd ${s}_hapmers
  bsub -q multicore40 -n 10 -R "rusage[mem=50000] span[hosts=1]" -M 60000 -oo ${s}.log -eo ${s}.err "
    export MERQURY='/srv/netscratch/dep_mercier/grp_schneeberger/software/merqury-1.1'
    meryl histogram ../${s}.meryl > ${s}.hist
    $MERQURY/trio/hapmers.sh $curmer $oramer ../${s}.meryl
  "
  cd ..
done

## K-mer haplotyping accuracy measurement using KMC ------------------------------------------------
# kmc version needed: K-Mer Counter (KMC) ver. 3.1.1 (2019-05-19)
# /srv/netscratch/dep_mercier/grp_schneeberger/software/kmc3//kmc
#
## kmer count in parental samples
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/
cd $cwd
samples=('currot' 'orangered')
for sample in ${samples[@]}; do
  cd $sample
  ls *trimmed-pair*.fastq.gz >fastq_list
  mkdir kmc_21mers_tmp
  bsub -q multicore20 -n 20 -R "rusage[mem=51200]" -M 61200 -eo kmc_21mers.err -oo kmc_21mers.log "
    kmc -k21 \
        -m50 \
        -fq \
        -t20 \
        -ci1 -cx500 \
        -n1000 \
        -r \
        @fastq_list ${sample}_21mers ./kmc_21mers_tmp/
    kmc_tools  -t1 transform ${sample}_21mers histogram ${sample}_kmc_21mers.histo dump ${sample}_kmc_21mers.txt -ci3 -cx200
  "
  cd ..
done

## find parent-specific kmers
#
# Currot-specific: coverage min = 6, coverage max = 52
fg=currot
bg=orangered
cd $cwd/${fg}
bsub -q normal -R "rusage[mem=8000]" -M 8000 -e ${fg}_only_kmers.err -o ${fg}_only_kmers.log "
  kmc_tools -t1 simple ${fg}_21mers ../${bg}/${bg}_21mers kmers_subtract ${fg}_only_kmers_specific_to_${bg} -ci6 -cx52
"
# Currot-specific: coverage min = 6, coverage max = 52
fg=orangered
bg=currot
cd $cwd/${fg}
bsub -q normal -R "rusage[mem=8000]" -M 8000 -e ${fg}_only_kmers.err -o ${fg}_only_kmers.log "
  kmc_tools -t1 simple ${fg}_21mers ../${bg}/${bg}_21mers kmers_subtract ${fg}_only_kmers_specific_to_${bg} -ci6 -cx52
"

## kmer count for clustered haplotigs/reads from gamete-binning
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads/'
cd $cwd
for i in {1..8}; do
  for s in PPP MMM; do
    mkdir kmc_21mers_tmp_${i}_${s}
    bsub -q multicore20 -n 10 -R "rusage[mem=51200] span[hosts=1]" -M 61200 -eo kmc_21mers_${i}_${s}.err -oo kmc_21mers_${i}_${s}.log "
      kmc -k21 \
          -r \
          -m50 \
          -fm \
          -t10 \
          -ci1 -cx500 \
          -n1000 \
          ${i}.txt_${s}_pbreads.fa ${i}_${s}_21mers ./kmc_21mers_tmp_${i}_${s}
    kmc_tools  -t1 transform ${i}_${s}_21mers histogram ${i}_${s}_kmc_21mers.histo dump ${i}_${s}_kmc_21mers.txt -ci3 -cx200
    "
  done
done

## kmer count for clustered haplotigs/reads from trio-binning
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/trio_binning/'
cd $cwd
for s in cur ora; do
  mkdir kmc_21mers_${s}
  bsub -q multicore40 -n 10 -R "rusage[mem=51200] span[hosts=1]" -M 61200 -oo ${s}_kmc.log -eo ${s}_kmc.err "
    kmc -k21 \
      -r \
      -m50 \
      -fm \
      -t10 \
      -ci1 -cx500 \
      -n1000 \
      haplotype-${s}.fasta.gz ${s}_21mers ./kmc_21mers_${s}
    kmc_tools -t1 transform ${s}_21mers histogram ${s}_kmc_21mers.histo dump ${s}_kmc_21mers.txt -ci3 -cx500
  "
done

## get intersection of kmers in reads and parents
# For Gamete-binning
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads/'
samples=('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered_only_kmers_specific_to_currot' '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_only_kmers_specific_to_orangered')
cd $cwd
for i in {1..8}; do
  for s in PPP MMM; do
    for sample in ${samples[@]}; do
      mkdir intersect_${i}_${s}
      cd intersect_${i}_${s}
      if [[ $sample == ${samples[1]} ]]; then
        parent='currot'
      else
        parent='orangered'
      fi
      bsub -q normal -R "rusage[mem=10000]" -M 10000 -eo ${i}_${s}_${parent}.err -oo ${i}_${s}_${parent}.log "
          kmc_tools -t1 simple $sample ../${i}_${s}_21mers -ci5 -cx50 \
            intersect ${i}_${s}_${parent}_21mers;
          kmc_tools -t1 transform ${i}_${s}_${parent}_21mers histogram ${i}_${s}_${parent}_21mers.histo dump ${i}_${s}_${parent}_21mers.txt -ci1;
           wc -l ${i}_${s}_${parent}_21mers.txt > ${i}_${s}_${parent}_21mers.count
        "
      cd ..
    done
  done
done

# For trio-binning
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/trio_binning/'
samples=('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered_only_kmers_specific_to_currot' '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_only_kmers_specific_to_orangered')
cd $cwd
for s in cur ora; do
  for sample in ${samples[@]}; do
    mkdir intersect_${s}
    cd intersect_${s}
    if [[ $sample == ${samples[1]} ]]; then
      parent='currot'
    else
      parent='orangered'
    fi
    bsub -q normal -R "rusage[mem=10000]" -M 10000 -eo ${s}_${parent}.err -oo ${s}_${parent}.log "
          kmc_tools -t1 simple $sample ../${s}_21mers -ci5 -cx50 \
            intersect ${s}_${parent}_21mers;
          kmc_tools -t1 transform ${s}_${parent}_21mers histogram ${s}_${parent}_21mers.histo dump ${s}_${parent}_21mers.txt -ci1;
           wc -l ${s}_${parent}_21mers.txt > ${s}_${parent}_21mers.count
        "
    cd ..
  done
done

## Create phasing quality summary plots
using commands in /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/compare_triobin_and_gametebin.py

# CONCLUSION READ PHASING: The purged assembly from the genome biology paper is not good enough to phase the HiFi reads properly.
# So read phasing of trio-binning is used. Also, it is not sure whether the reads are mapped correctly, so gamete-binning linkage
# group cannot be used at the read level.

#####################################################################
### Step 2: Assemble trio-binning phased reads
#####################################################################

## Assemble phased reads

# Using Hifiasm
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'

cd $cwd
for s in cur ora; do
  mkdir $s
  cd $s
  bsub -q multicore40 -n 40 -R "rusage[mem=20000] span[hosts=1]" -M 20000 -oo ${s}.log -eo ${s}.err "
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
#    hifiasm -t 40 -o $s ${indir}/haplotype-${s}.fasta.gz
#    cat $s.p_ctg.gfa | grep '^S' | cut -f2,3 | awk '{print \">\"\$1\"\\n\"\$2}' > $s.p_ctg.fasta
    cat hifiasm.p_utg.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' >hifiasm.p_utg.fasta
    hometools seqsize $s.p_ctg.fasta | tail +2 > $s.p_ctg.fasta.chrsize

    hifiasm -t 40 -o ${s}_unk ${indir}/haplotype-${s}.fasta.gz ${indir}/haplotype-unknown.fasta.gz
    cat ${s}_unk.p_ctg.gfa | grep '^S' | cut -f2,3 | awk '{print \">\"\$1\"\\n\"\$2}' > ${s}_unk.p_ctg.fasta
    cat cur_unk.p_utg.gfa | grep '^S' | cut -f2,3 | awk '{print ">"$1"\n"$2}' > TEST_cur_unk.p_utg.fasta
    hometools seqsize ${s}_unk.p_ctg.fasta | tail +2 > ${s}_unk.p_ctg.fasta.chrsize
  "
  cd ..
done

# Using Canu: Use the assembly generated during testing
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/canu_assembly/
ln -s ../../test/canu_trio/canu_trio/canu_trio-haplotypecur/canu_trio-haplotypecur.contigs.fasta .
ln -s ../../test/canu_trio/canu_trio/canu_trio-haplotypecur/canu_trio-haplotypecur.unitigs.fasta .

ln -s ../../test/canu_trio/canu_trio/canu_trio-haplotypeora/canu_trio-haplotypeora.contigs.fasta .
ln -s ../../test/canu_trio/canu_trio/canu_trio-haplotypeora/canu_trio-haplotypeora.unitigs.fasta .

## Contig-filtering: Remove contigs for which the read depth (from remapping of PB reads) is too low

# Hifiasm test

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
# I used wtB reads to test assembly, but it is better to use the parents reads rather than RP. In later analysis wtB is not used.
wtBR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/rp_leaf_illumina/wt_B/wt_B_ql-trimmed-pair1.fastq.gz'
wtBR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/rp_leaf_illumina/wt_B/wt_B_ql-trimmed-pair2.fastq.gz'

for s in cur ora; do
  cd $s
  t=20
  bsub -q multicore40 -n 20 -R "rusage[mem=40000] span[hosts=1]" -M 40000 -oo ${s}_unk.log -eo ${s}_unk.err "
    # Align PACBIO reads and get read-depth
    minimap2 -ax map-pb -t 10 ${s}_unk.p_ctg.fasta ${indir}/haplotype-${s}.fasta.gz ${indir}/haplotype-unknown.fasta.gz \
    | samtools sort -O BAM -@ 10 \
    > ${s}_unk.sorted.bam
    samtools index ${s}_unk.sorted.bam
    samtools view -F 256 -F 2048 -h -@ 10 ${s}_unk.sorted.bam \
    | grep -v 'SA:' \
    | samtools view -O BAM -@ 10 - \
    > ${s}_unk_F256_F2048_noSA.sorted.bam
    samtools view -q 10 -O BAM -@ 10 ${s}_unk_F256_F2048_noSA.sorted.bam > ${s}_unk_q10_F256_F2048_noSA.sorted.bam
    samtools index ${s}_unk_q10_F256_F2048_noSA.sorted.bam
    samtools depth -a ${s}_unk_q10_F256_F2048_noSA.sorted.bam > ${s}_unk.depth

#    # Align Illumina reads from WT_B and get read-depth
#    minimap2 -ax sr -t $t ${s}_unk.p_ctg.fasta ${wtBR1} ${wtBR2} \
#    | samtools sort -O BAM -@ $t - \
#    > ${s}_unk_wtB.sorted.bam

    samtools sort -O BAM -@ $t -n ${s}_unk_wtB.sorted.bam  \
    | samtools fixmate -@ $t -c -m -O SAM - - \
    | samtools sort -O BAM -@ $t - \
    | samtools markdup -@ $t -S -s -O SAM - - \
    | samtools view -F 1024 -F 256 -F 2048 -h -O SAM -@ $t - \
    | samtools sort -O SAM -@ $t - \
    | grep -v 'SA:' \
    | samtools view -O BAM -h -@ $t - \
    > ${s}_unk_wtB_F256_F1024_F2048.deduped.sorted.bam

    samtools view -q 10 -O BAM -@ 60 ${s}_unk_wtB_F256_F1024_F2048.deduped.sorted.bam > ${s}_unk_wtB_q10_F256_F1024_F2048.deduped.sorted.bam
    samtools index ${s}_unk_wtB_q10_F256_F1024_F2048.deduped.sorted.bam
    samtools depth -a ${s}_unk_wtB_q10_F256_F1024_F2048.deduped.sorted.bam > ${s}_unk_wtB.depth

  "
  cd ..
done


# Assemblies generated using hifiasm had better assembly quality and are selected for further processing.
# Testing Cur_run3 which uses cur_unk.p_utg as the primary contig. This is because in the contig assembly some duplicated regions were collapsed.
# Finally, I decided to use the unitig assembly output from Hifiasm.. so,I deleted all the commands corresponding to contig assembly processing to keep this file clean.

# Hifiasm test
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
curR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair1.fastq.gz'
curR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair2.fastq.gz'

bsub -q multicore40 -n 40 -R "rusage[mem=40000] span[hosts=1]" -M 40000 -oo pb_unk.log -eo pb_unk.err "
  # Align PACBIO reads and get read-depth
  minimap2 -ax map-pb -t 40 -R '@RG\tID:trio_cur\tSM:trio_cur' cur_unk.p_utg.fasta ${indir}/haplotype-cur.fasta.gz \
    | samtools sort -O BAM -@ 40 \
    > cur_unk_cur.sorted.bam
  samtools index cur_unk_cur.sorted.bam

  minimap2 -ax map-pb -t 40 -R '@RG\tID:trio_ora\tSM:trio_ora' cur_unk.p_utg.fasta ${indir}/haplotype-ora.fasta.gz \
      | samtools sort -O BAM -@ 40 \
      > cur_unk_ora.sorted.bam
  samtools index cur_unk_ora.sorted.bam

  minimap2 -ax map-pb -t 40 -R '@RG\tID:trio_unk\tSM:trio_unk' cur_unk.p_utg.fasta ${indir}/haplotype-unknown.fasta.gz \
      | samtools sort -O BAM -@ 40 \
      > cur_unk_unk.sorted.bam
  samtools index cur_unk_unk.sorted.bam

  samtools merge -@40 -f -O BAM cur_unk.sorted.bam cur_unk_cur.sorted.bam cur_unk_unk.sorted.bam
  samtools index cur_unk.sorted.bam

  samtools view -F 256 -F 2048 -h -@ 40 cur_unk.sorted.bam \
  | grep -v 'SA:' \
  | samtools view -O BAM -@ 40 - \
  > cur_unk_F256_F2048_noSA.sorted.bam
  samtools view -q 10 -O BAM -@ 40 cur_unk_F256_F2048_noSA.sorted.bam > cur_unk_q10_F256_F2048_noSA.sorted.bam
  samtools index cur_unk_q10_F256_F2048_noSA.sorted.bam
  samtools depth -a cur_unk_q10_F256_F2048_noSA.sorted.bam > cur_unk.depth
"

bsub -q multicore40 -n 40 -R "rusage[mem=40000] span[hosts=1]" -M 40000 -oo ill_unk.log -eo ill_unk.err "
  # Align Illumina reads from WT_B and get read-depth
  minimap2 -ax sr -t 40 -R '@RG\tID:cursr\tSM:cursr' cur_unk.p_utg.fasta ${curR1} ${curR2} \
  | samtools sort -O BAM -@ 40 - \
  > cur_unk_cursr.sorted.bam

  samtools sort -O BAM -@ 40 -n cur_unk_cursr.sorted.bam  \
  | samtools fixmate -@ 40 -c -m -O SAM - - \
  | samtools sort -O BAM -@ 40 - \
  | samtools markdup -@ 40 -S -s -O SAM - - \
  | samtools view -F 1024 -F 256 -F 2048 -h -O SAM -@ 40 - \
  | samtools sort -O SAM -@ 40 - \
  | grep -v 'SA:' \
  | samtools view -O BAM -h -@ 40 - \
  > cur_unk_cursr_F256_F1024_F2048.deduped.sorted.bam

  samtools view -q 10 -O BAM -@ 40 cur_unk_cursr_F256_F1024_F2048.deduped.sorted.bam > cur_unk_cursr_q10_F256_F1024_F2048.deduped.sorted.bam
  samtools index cur_unk_cursr_q10_F256_F1024_F2048.deduped.sorted.bam
#  samtools depth -a cur_unk_cursr_q10_F256_F1024_F2048.deduped.sorted.bam > cur_unk_cursr.depth
  samtools depth -a cur_unk_cursr_F256_F1024_F2048.deduped.sorted.bam > cur_unk_cursr_noMQfilter.depth
"
# Get assembly stats using read_depth_summary function in the analysis_plots.py file
# CUR + HiFi_mapping
df = read_depth_summary('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk.depth',
                        '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk.p_utg.fasta.chrsize',
                        '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk_rd_stats.pdf')
df.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk.stats', sep='\t', index=False, header = ['contig', 'size', 'mean_read_depth', 'cummulative_assembly_size'])

# Cur + Cur_illumina
df = read_depth_summary('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk_cursr.depth',
                   '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk.p_utg.fasta.chrsize',
                   '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk_cursr_rd_stats.pdf', min=100)
df.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk_cursr.stats', sep='\t', index=False, header = ['contig', 'size', 'mean_read_depth', 'cummulative_assembly_size'])

# Cur + Cur_illumina_noMQ_filter
df = read_depth_summary('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk_cursr_noMQfilter.depth',
                   '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk.p_utg.fasta.chrsize',
                   '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk_cursr_noMQfilter_rd_stats.pdf', min=100)
df.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/cur_unk_cursr_noMQfilter.stats', sep='\t', index=False, header = ['contig', 'size', 'mean_read_depth', 'cummulative_assembly_size'])


## Get corrected Cur assembly
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
curR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair1.fastq.gz'
curR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair2.fastq.gz'


tail +2 cur_unk_cursr_noMQfilter.stats | awk '{if($3<30 || $2<100000){print $1}}' > candidate_save.contigs.v1.txt
hometools getchr -v -F candidate_save.contigs.v1.txt -o filtered.contigs.v1.fasta cur_unk.p_utg.fasta
hometools getchr -F candidate_save.contigs.v1.txt -o candidate_save.contigs.v1.fasta cur_unk.p_utg.fasta

hometools seqsize candidate_save.contigs.v1.fasta > candidate_save.contigs.v1.fasta.chrsize
hometools seqsize filtered.contigs.v1.fasta > filtered.contigs.v1.fasta.chrsize

minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' filtered.contigs.v1.fasta ${curR1} ${curR2} \
  | samtools sort -O BAM -@ 70 - \
  > filtered.contigs.v1.cursr.sorted.bam
samtools index filtered.contigs.v1.cursr.sorted.bam

# Check overlap of candidate contigs
cd minimap2_check
while read r; do
  bsub -q short -n 1 -R "rusage[mem=8000] span[hosts=1]" -M 10000 -oo ${r}.log -eo ${r}.err "
    hometools getchr -o ${r}.fasta ../cur_unk.p_utg.fasta --chrs ${r}
    minimap2 -x asm10 -c --eqx -t 1 ../filtered.contigs.v1.fasta ${r}.fasta > ${r}.fasta.paf
  "
  echo $r
done < ../candidate_save.contigs.v1.txt
indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/minimap2_check/'
minimap2_filt(indir, indir + 'contigs.overlap.pdf', indir + 'contigs.low_overlap.txt')

# get unmapped reads
samtools view -u -@ 60 -f 4  -F 264 filtered.contigs.v1.cursr.sorted.bam  > v1.tmp1.bam
samtools view -u -@ 60 -f 8  -F 260 filtered.contigs.v1.cursr.sorted.bam  > v1.tmp2.bam
samtools view -u -@ 60 -f 12 -F 256 filtered.contigs.v1.cursr.sorted.bam  > v1.tmp3.bam

samtools merge -u - v1.tmp[123].bam | samtools sort -n -@ 60 -O bam -o unmapped.bam -
samtools fastq -@ 60 -1 unmapped_R1.fastq -2 unmapped_R2.fastq unmapped.bam

# map unmapped reads to the candidate save contigs
minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' candidate_save.contigs.v1.fasta unmapped_R1.fastq unmapped_R2.fastq \
  | samtools sort -O BAM -@ 70 - \
  > candidate.contigs.v1.unmapped_cursr.sorted.bam
samtools index candidate.contigs.v1.unmapped_cursr.sorted.bam
samtools depth -a candidate.contigs.v1.unmapped_cursr.sorted.bam > candidate.contigs.v1.unmapped_cursr.sorted.depth

indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/'
read_coverage(indir + 'candidate.contigs.v1.unmapped_cursr.sorted.depth', indir + 'candidate.contigs.v1.unmapped_cursr.sorted.depth.pdf', indir + 'minimap2_check/contigs.low_overlap.txt')


# Select bad contigs (low coverage) manually from the candidate.contigs.v1.unmapped_cursr.sorted.depth.pdf
# Saved in remove_contigs.v1.txt
cat contigs.txt remove_contigs.v1.txt | sort | uniq -c | grep "1 " | cut -d' ' -f8 > saved.contigs.v1.txt
hometools getchr -F saved.contigs.v1.txt -o saved.contigs.v1.fasta cur_unk.p_utg.fasta

# Merge the contigs
cat filtered.contigs.v1.fasta saved.contigs.v1.fasta > cur.v1.fasta


# Second round of CUR assembly correction
minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' cur.v1.fasta ${curR1} ${curR2} \
  | samtools sort -O BAM -@ 70 - \
  > cur.v1.sorted.bam
samtools index cur.v1.sorted.bam
t=60
samtools sort -O BAM -@ 60 -n cur.v1.sorted.bam  \
| samtools fixmate -@ 60 -c -m -O SAM - - \
| samtools sort -O BAM -@ $t - \
| samtools markdup -@ $t -S -s -O SAM - - \
| samtools view -F 1024 -F 256 -F 2048 -h -O SAM -@ $t - \
| samtools sort -O SAM -@ $t - \
| grep -v 'SA:' \
| samtools view -O BAM -h -@ $t - \
> cur.v1.F256_F1024_F2048.deduped.sorted.bam
samtools depth -a cur.v1.F256_F1024_F2048.deduped.sorted.bam > cur.v1.noMQfilter.depth

indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/'
df = read_depth_summary(indir + 'cur.v1.noMQfilter.depth',
                        indir + 'cur.v1.fasta.chrsize',
                        indir + 'cur.v1.cursr.pdf', min=100)
df.to_csv(indir + 'cur.v1.cursr.stats', sep='\t', index=False, header = ['contig', 'size', 'mean_read_depth', 'cummulative_assembly_size'])

tail +2 cur.v1.cursr.stats | awk '{if($3<34 || $2<100000){print $1}}' > candidate_save.contigs.v2.txt
hometools getchr -v -F candidate_save.contigs.v2.txt -o filtered.contigs.v2.fasta cur.v1.fasta
hometools getchr -F candidate_save.contigs.v2.txt -o candidate_save.contigs.v2.fasta cur_unk.p_utg.fasta

hometools seqsize candidate_save.contigs.v2.fasta > candidate_save.contigs.v2.fasta.chrsize
hometools seqsize filtered.contigs.v2.fasta > filtered.contigs.v2.fasta.chrsize

minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' filtered.contigs.v2.fasta ${curR1} ${curR2} \
  | samtools sort -O BAM -@ 70 - \
  > filtered.contigs.v2.cursr.sorted.bam
samtools index filtered.contigs.v2.cursr.sorted.bam


# Check overlap of candidate contigs
cd minimap2_check2
while read r; do
  bsub -q short -n 1 -R "rusage[mem=8000] span[hosts=1]" -M 10000 -oo ${r}.log -eo ${r}.err "
    hometools getchr -o ${r}.fasta ../cur.v1.fasta --chrs ${r}
    minimap2 -x asm10 -c --eqx -t 1 ../filtered.contigs.v2.fasta ${r}.fasta > ${r}.fasta.paf
  "
  echo $r
done < ../candidate_save.contigs.v2.txt
indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/minimap2_check2/'
minimap2_filt(indir, indir + 'contigs.overlap.pdf', indir + 'contigs.low_overlap.txt')


# get unmapped reads
samtools view -u -@ 60 -f 4  -F 264 filtered.contigs.v2.cursr.sorted.bam  > v2.tmp1.bam
samtools view -u -@ 60 -f 8  -F 260 filtered.contigs.v2.cursr.sorted.bam  > v2.tmp2.bam
samtools view -u -@ 60 -f 12 -F 256 filtered.contigs.v2.cursr.sorted.bam  > v2.tmp3.bam

samtools merge -u - v2.tmp[123].bam | samtools sort -n -@ 60 -O bam -o unmapped2.bam -
samtools fastq -@ 60 -1 unmapped2_R1.fastq -2 unmapped2_R2.fastq unmapped2.bam

# map unmapped reads to the candidate save contigs
minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' candidate_save.contigs.v2.fasta unmapped2_R1.fastq unmapped2_R2.fastq \
  | samtools sort -O BAM -@ 70 - \
  > candidate.contigs.v2.unmapped_cursr.sorted.bam
samtools index candidate.contigs.v2.unmapped_cursr.sorted.bam
samtools depth -a candidate.contigs.v2.unmapped_cursr.sorted.bam > candidate.contigs.v2.unmapped_cursr.sorted.depth

indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/'
read_coverage(indir + 'candidate.contigs.v2.unmapped_cursr.sorted.depth', indir + 'candidate.contigs.v2.unmapped_cursr.sorted.depth.pdf', indir + 'minimap2_check2/contigs.low_overlap.txt', y_cut=20)

# Select bad contigs (low coverage) manually from the candidate.contigs.v1.unmapped_cursr.sorted.depth.pdf
# Saved in remove_contigs.v1.txt
cut -f1 minimap2_check2/contigs.low_overlap.txt > contigs2.txt
cat contigs2.txt remove_contigs.v2.txt | sort | uniq -c | grep "1 " | cut -d' ' -f8 > saved.contigs.v2.txt
hometools getchr -F saved.contigs.v2.txt -o saved.contigs.v2.fasta cur_unk.p_utg.fasta

# Merge the contigs
cat filtered.contigs.v2.fasta saved.contigs.v2.fasta > cur.v2.fasta
hometools seqsize cur.v2.fasta > cur.v2.fasta.chrsize

# Third round of CUR assembly correction
minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' cur.v2.fasta ${curR1} ${curR2} \
  | samtools sort -O BAM -@ 70 - \
  > cur.v2.sorted.bam
samtools index cur.v2.sorted.bam
t=70
samtools sort -O BAM -@ $t -n cur.v2.sorted.bam  \
| samtools fixmate -@ $t -c -m -O SAM - - \
| samtools sort -O BAM -@ $t - \
| samtools markdup -@ $t -S -s -O SAM - - \
| samtools view -F 1024 -F 256 -F 2048 -h -O SAM -@ $t - \
| samtools sort -O SAM -@ $t - \
| grep -v 'SA:' \
| samtools view -O BAM -h -@ $t - \
> cur.v2.F256_F1024_F2048.deduped.sorted.bam
samtools depth -a cur.v2.F256_F1024_F2048.deduped.sorted.bam > cur.v2.noMQfilter.depth

indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/'
df = read_depth_summary(indir + 'cur.v2.noMQfilter.depth',
                        indir + 'cur.v2.fasta.chrsize',
                        indir + 'cur.v2.cursr.pdf', min=100)
df.to_csv(indir + 'cur.v2.cursr.stats', sep='\t', index=False, header = ['contig', 'size', 'mean_read_depth', 'cummulative_assembly_size'])


tail +2 cur.v2.cursr.stats | awk '{if($3<30){print $1}}' > candidate_save.contigs.v3.txt
hometools getchr -v -F candidate_save.contigs.v3.txt -o filtered.contigs.v3.fasta cur.v2.fasta
hometools getchr -F candidate_save.contigs.v3.txt -o candidate_save.contigs.v3.fasta cur.v2.fasta
hometools seqsize filtered.contigs.v3.fasta > filtered.contigs.v3.fasta.chrsize

minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' filtered.contigs.v3.fasta ${curR1} ${curR2} \
  | samtools sort -O BAM -@ 70 - \
  > filtered.contigs.v3.cursr.sorted.bam

## Save contigs that were removed in the previous step but have high coverage for reads not mapping with cur.v2.fasta
samtools view -u -@ 60 -f 4  -F 264 cur.v2.sorted.bam  > v3.tmp1.bam
samtools view -u -@ 60 -f 8  -F 260 cur.v2.sorted.bam  > v3.tmp2.bam
samtools view -u -@ 60 -f 12 -F 256 cur.v2.sorted.bam  > v3.tmp3.bam
samtools merge -u - v3.tmp[123].bam | samtools sort -n -@ 60 -O bam -o unmapped3.bam -
samtools fastq -@ 60 -1 unmapped3_R1.fastq -2 unmapped3_R2.fastq unmapped3.bam
hometools getchr -F remove_contigs.v2.txt -o remove_contigs.v2.fasta cur_unk.p_utg.fasta
hometools seqsize remove_contigs.v2.fasta > remove_contigs.v2.fasta.chrsize
minimap2 -ax sr -t 70 -R '@RG\tID:cursr\tSM:cursr' remove_contigs.v2.fasta unmapped3_R1.fastq unmapped3_R2.fastq \
  | samtools sort -O BAM -@ 70 - \
  > remove_contigs.v2.sorted.bam
samtools index remove_contigs.v2.sorted.bam
samtools depth -a remove_contigs.v2.sorted.bam > remove_contigs.v2.sorted.depth
indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/'
df = read_depth_summary(indir + 'remove_contigs.v2.sorted.depth',
                        indir + 'remove_contigs.v2.fasta.chrsize',
                        indir + 'remove_contigs.v2.unmapped3.pdf', min=100)
df.to_csv(indir + 'remove_contigs.v2.unmapped3.stats', sep='\t', index=False, header = ['contig', 'size', 'mean_read_depth', 'cummulative_assembly_size'])
awk  'NR>1 {if($3>10){print $1}}' remove_contigs.v2.unmapped3.stats > saved.contigs.v3.txt
hometools getchr -F saved.contigs.v3.txt -o saved.contigs.v3.fasta cur_unk.p_utg.fasta
hometools seqsize saved.contigs.v3.fasta > saved.contigs.v3.fasta.chrsize


cat saved.contigs.v3.fasta >> candidate_save.contigs.v3.fasta
hometools seqsize candidate_save.contigs.v3.fasta > candidate_save.contigs.v3.fasta.chrsize



################################################################################
# Filter orangered assembly
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
oraR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered_ql-trimmed-pair1.fastq.gz'
oraR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered_ql-trimmed-pair2.fastq.gz'

bsub -q multicore40 -n 40 -R "rusage[mem=40000] span[hosts=1]" -M 40000 -oo ill_unk.log -eo ill_unk.err "
  # Align Illumina reads and get read-depth
  minimap2 -ax sr -t 40 -R '@RG\tID:orasr\tSM:orasr' ora_unk.p_utg.fasta ${oraR1} ${oraR2} \
  | samtools sort -O BAM -@ 40 - \
  > ora_unk_orasr.sorted.bam

  samtools sort -O BAM -@ 40 -n ora_unk_orasr.sorted.bam  \
  | samtools fixmate -@ 40 -c -m -O SAM - - \
  | samtools sort -O BAM -@ 40 - \
  | samtools markdup -@ 40 -S -s -O SAM - - \
  | samtools view -F 1024 -F 256 -F 2048 -h -O SAM -@ 40 - \
  | samtools sort -O SAM -@ 40 - \
  | grep -v 'SA:' \
  | samtools view -O BAM -h -@ 40 - \
  > ora_unk_orasr_F256_F1024_F2048.deduped.sorted.bam

  samtools index ora_unk_orasr_F256_F1024_F2048.deduped.sorted.bam
  samtools depth -a ora_unk_orasr_F256_F1024_F2048.deduped.sorted.bam > ora_unk_orasr_noMQfilter.depth
"

indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/'
# Cur + Cur_illumina_noMQ_filter
df = read_depth_summary(indir + 'ora_unk_orasr_noMQfilter.depth',
                        indir + 'ora_unk.p_utg.fasta.chrsize',
                        indir + 'ora_unk_orasr_noMQfilter_rd_stats.pdf', min=100)
df.to_csv(indir + 'ora_unk_orasr_noMQfilter.stats', sep='\t', index=False, header = ['contig', 'size', 'mean_read_depth', 'cummulative_assembly_size'])


tail +2 ora_unk_orasr_noMQfilter.stats | awk '{if($3<30 || $2<100000){print $1}}' > candidate_save.contigs.v1.txt
hometools getchr -v -F candidate_save.contigs.v1.txt -o filtered.contigs.v1.fasta ora_unk.p_utg.fasta
hometools getchr -F candidate_save.contigs.v1.txt -o candidate_save.contigs.v1.fasta ora_unk.p_utg.fasta

hometools seqsize candidate_save.contigs.v1.fasta > candidate_save.contigs.v1.fasta.chrsize
hometools seqsize filtered.contigs.v1.fasta > filtered.contigs.v1.fasta.chrsize

minimap2 -ax sr -t 70 -R '@RG\tID:orasr\tSM:orasr' filtered.contigs.v1.fasta ${oraR1} ${oraR2} \
  | samtools sort -O BAM -@ 70 - \
  > filtered.contigs.v1.orasr.sorted.bam
samtools index filtered.contigs.v1.orasr.sorted.bam

# Check overlap of candidate contigs
cd minimap2_check
while read r; do
  bsub -q short -n 1 -R "rusage[mem=8000] span[hosts=1]" -M 10000 -oo ${r}.log -eo ${r}.err "
    hometools getchr -o ${r}.fasta ../ora_unk.p_utg.fasta --chrs ${r}
    minimap2 -x asm10 -c --eqx -t 1 ../filtered.contigs.v1.fasta ${r}.fasta > ${r}.fasta.paf
  "
  echo $r
done < ../candidate_save.contigs.v1.txt
indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/minimap2_check/'
minimap2_filt(indir, indir + 'contigs.overlap.pdf', indir + 'contigs.low_overlap.txt')

samtools view -u -@ 60 -f 4  -F 264 filtered.contigs.v1.orasr.sorted.bam  > v1.tmp1.bam
samtools view -u -@ 60 -f 8  -F 260 filtered.contigs.v1.orasr.sorted.bam  > v1.tmp2.bam
samtools view -u -@ 60 -f 12 -F 256 filtered.contigs.v1.orasr.sorted.bam  > v1.tmp3.bam

samtools merge -u - v1.tmp[123].bam | samtools sort -n -@ 60 -O bam -o unmapped.bam -
samtools fastq -@ 60 -1 unmapped_R1.fastq -2 unmapped_R2.fastq unmapped.bam

# map unmapped reads to the candidate save contigs
minimap2 -ax sr -t 70 -R '@RG\tID:orasr\tSM:orasr' candidate_save.contigs.v1.fasta unmapped_R1.fastq unmapped_R2.fastq \
  | samtools sort -O BAM -@ 70 - \
  > candidate.contigs.v1.unmapped_orasr.sorted.bam
samtools depth -a candidate.contigs.v1.unmapped_orasr.sorted.bam > candidate.contigs.v1.unmapped_orasr.sorted.depth

indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/'
read_coverage(indir + 'candidate.contigs.v1.unmapped_orasr.sorted.depth', indir + 'candidate.contigs.v1.unmapped_orasr.sorted.depth.pdf', indir + 'minimap2_check/contigs.low_overlap.txt')



################################################################################
# OLD CODE: Was trying to use BLAST to find contigs that are overlapping with good contigs, but later switched to minimap for this.

cd blast
xargs -a ../bad_contigs.txt -n 1 -P 6 -I {} bash -c '
  hometools getchr -o ${1}.fasta ../cur_unk.p_utg.fasta --chrs ${1}
  blastn -query ${1}.fasta \
   -task megablast \
   -db ../good_contigs.fasta \
   -out ${1}.oblast \
   -outfmt 7 \
   -max_target_seqs 10 \
   -max_hsps 500 \
   -perc_identity 95 \
   -evalue 0.001 \
   -num_threads 20 \
   >> blastall_blastn.log
' -- {}

from get_bad_contig_coverage_in_good_contigs.py import get_stats, get_stats_plot
D = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/blast'
B = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/bad_contigs.fasta.chrsize'
P = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/bad_contigs.fasta' + ".blast.sum"

get_stats('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/blast', '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/bad_contigs.fasta.chrsize',
'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/bad_contigs.fasta.blast.sum')
get_stats_plot('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/bad_contigs.fasta.blast.sum',
'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/bad_contigs.fasta.blast.sum.pdf')



makeblastdb -dbtype nucl -in filtered.contigs.v3.fasta -input_type fasta
cd blast2
nohup xargs -a ../contigs_cursr_DP_lower25_long.txt -n 1 -P 6 -I {} bash -c '
#  hometools getchr -o ${1}.fasta ../cur_unk.p_utg.fasta --chrs ${1}
  blastn -query ${1}.fasta \
   -task megablast \
   -db ../filtered.contigs.v3.fasta \
   -out ${1}.oblast \
   -outfmt 7 \
   -max_target_seqs 10 \
   -max_hsps 500 \
   -perc_identity 90 \
   -evalue 0.0001 \
   -num_threads 40 \
   >> blastall_blastn.log
' -- {} &

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/'
get_stats(indir + 'blast2',
  indir + 'candidate_save.contigs.v3.fasta.chrsize',
  indir + 'blast2/bad_contigs.fasta.blast.sum')
get_stats_plot(indir + 'blast2/bad_contigs.fasta.blast.sum',
  indir + 'blast2/bad_contigs.fasta.blast.sum.length.pdf',
  indir + 'blast2/bad_contigs.fasta.blast.sum.overlap.pdf')


cd blast3
while read r; do
  hometools getchr -o ${r}.fasta ../cur_unk.p_utg.fasta --chrs ${r}
  bsub -q multicore40 -n 40 -R "rusage[mem=25000] span[hosts=1]" -M 25000 -oo ${r}.log -eo ${r}.err "
    blastn -query ${r}.fasta \
     -task dc-megablast \
     -db ../filtered.contigs.v3.fasta \
     -out ${r}.oblast \
     -outfmt 7 \
     -max_target_seqs 100 \
     -max_hsps 100 \
     -perc_identity 95 \
     -evalue 0.0001 \
     -num_threads 60 \
     >> blastall_blastn.log
  "
  echo $r
done < ../contigs_cursr_DP_lower25_long.txt

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/'
get_stats(indir + 'blast3',
  indir + 'candidate_save.contigs.v3.fasta.chrsize',
  indir + 'blast3/bad_contigs.fasta.blast.sum')
get_stats_plot(indir + 'blast3/bad_contigs.fasta.blast.sum',
  indir + 'blast3/bad_contigs.fasta.blast.sum.length.pdf',
  indir + 'blast3/bad_contigs.fasta.blast.sum.overlap.pdf')





cd trash
minimap2 -ax sr -t 60 -R '@RG\tID:cursr\tSM:cursr' variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_Currot_on_Original_v1.0.fasta ${curR1} ${curR2} \
  | samtools sort -O BAM -@ 60 - \
  > cur.v1.1.sorted.bam
  samtools index cur.v1.1.sorted.bam

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 50 -R '@RG\tID:trio_cur\tSM:trio_cur' variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_Currot_on_Original_v1.0.fasta ${indir}/haplotype-cur.fasta.gz \
    | samtools sort -O BAM -@ 50 \
    > cur.v1.1.cur.sorted.bam
  samtools index cur.v1.1.cur.sorted.bam


# Canu test
curR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair1.fastq.gz'
curR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair2.fastq.gz'
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/canu_assembly/
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
for s in ora; do
  t=20
  bsub -q multicore40 -n 20 -R "rusage[mem=40000] span[hosts=1]" -M 40000 -oo ${s}.log -eo ${s}.err "
    # Align PACBIO reads and get read-depth
#    minimap2 -ax map-pb -t 10 canu_trio-haplotype${s}.contigs.fasta ${indir}/haplotype-${s}.fasta.gz ${indir}/haplotype-unknown.fasta.gz \
#    | samtools sort -O BAM -@ 10 \
#    > ${s}_contigs.sorted.bam
#    samtools index ${s}_contigs.sorted.bam
#    samtools view -F 256 -F 2048 -h -@ 10 ${s}_contigs.sorted.bam \
#    | grep -v 'SA:' \
#    | samtools view -O BAM -@ 10 - \
#    > ${s}_contigs_F256_F2048_noSA.sorted.bam
#    samtools index ${s}_unk_F256_F2048_noSA.sorted.bam
#    samtools depth -a ${s}_contigs_F256_F2048_noSA.sorted.bam > ${s}_contigs.depth


    # Align Illumina reads from WT_B and get read-depth
    minimap2 -ax sr -t $t canu_trio-haplotype${s}.contigs.fasta ${curR1} ${curR2} \
    | samtools sort -O BAM -@ $t - \
    > ${s}_contigs_cursr.sorted.bam

    samtools sort -O BAM -@ $t -n ${s}_contigs_cursr.sorted.bam  \
    | samtools fixmate -@ $t -c -m -O SAM - - \
    | samtools sort -O BAM -@ $t - \
    | samtools markdup -@ $t -S -s -O SAM - - \
    | samtools view -F 1024 -F 256 -F 2048 -h -O SAM -@ $t - \
    | samtools sort -O SAM -@ $t - \
    | grep -v 'SA:' \
    | samtools view -O BAM -h -@ $t - \
    > ${s}_contigs_cursr_F256_F1024_F2048.deduped.sorted.bam

    samtools index ${s}_contigs_cursr_F256_F1024_F2048.deduped.sorted.bam
    samtools view -q 10 -O BAM -@ $t ${s}_contigs_cursr_F256_F1024_F2048.deduped.sorted.bam > ${s}_contigs_cursr_F256_F1024_F2048.deduped.sorted.bam
    samtools index ${s}_contigs_cursr_F256_F1024_F2048.deduped.sorted.bam
    samtools depth -a ${s}_contigs_cursr_F256_F1024_F2048.deduped.sorted.bam > ${s}_contigs_cursr.depth
  "
done


