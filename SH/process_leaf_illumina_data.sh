refgenome='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/rp.diploid.fasta.mm2_Xsr.idx'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/leaf_illumina/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/'
samples=( 'rojo_wt_A' 'rojo_wt_B' 'rojo_mut_C' 'rojo_mut_D' )

cd $cwd

for sample in ${samples[@]}; do
  bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo $sample.log -eo $sample.err "
    skewer -r 0.1 -d 0.05 -k 8 -q 20 -l 75 -m pe -t 20 -x /srv/netscratch/dep_mercier/grp_schneeberger/projects/hyper_co/data/reads/shqMergedAdapters_Primers_representative_rc.fa ${indir}/${sample}_R1.fastq.gz ${indir}/${sample}_R2.fastq.gz -z -o ${sample}_ql
    minimap2 -ax sr --eqx -t 20 $refgenome ${sample}_ql-trimmed-pair1.fastq.gz ${sample}_ql-trimmed-pair2.fastq.gz | samtools view -b - | samtools sort -O BAM -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam
    "
    done
