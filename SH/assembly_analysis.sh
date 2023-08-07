# Commands for analysing the currot and oragenred assemblies

# <editor-fold desc="Align contigs to the final assemblies">
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemblyplots/
# Align ora contigs to the ora genome
curcon=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/cur/cur.genomic_contigs.fasta
curasm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/cur/cur.genome.v1.fasta
cd $cwd
bsub -q multicore20 -n 10 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo cur.log -eo cur.err "
    minimap2 -ax asm5 -t 10 --secondary=no \
          $curasm \
          $curcon \
    | samtools sort -@ 10 -O BAM - \
    > cur.contig.sorted.bam
    samtools index cur.contig.sorted.bam
"
# Align ora contigs to the ora genome
oracon=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/ora/ora.genomic_contigs.fasta
oraasm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/ora/ora.genome.v1.fasta
cd $cwd
bsub -q multicore20 -n 10 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ora.log -eo ora.err "
    minimap2 -ax asm5 -t 10 --secondary=no \
          $oraasm \
          $oracon \
    | samtools sort -@ 10 -O BAM - \
    > ora.contig.sorted.bam
    samtools index ora.contig.sorted.bam
"
# </editor-fold>


# <editor-fold desc="Re-Run repeatmasker to get centromeric regions">
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemblyplots/centromerecheck/
# Copy the centromeric and rRNA sequence
cd $cwd
cp /netscratch/dep_mercier/grp_schneeberger/projects/Athal/MA_lines/4.annotation/45srDNA_RepeatMasker/Mock-A/rep1/NAR_bioRxiv/rDNA_NaishCEN_telomeres.fa .
# run Repeatmasker
curasm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/cur/cur.genome.v1.fasta
oraasm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/ora/ora.genome.v1.fasta

cd $cwd; mkdir -p cur; cd $cwd/cur
bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo cur.log -eo cur.err "
    RepeatMasker -lib ../rDNA_NaishCEN_telomeres.fa -gff -dir . -pa 40 $curasm
"
cd $cwd; mkdir -p ora; cd $cwd/ora
bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ora.log -eo ora.err "
    RepeatMasker -lib ../rDNA_NaishCEN_telomeres.fa -gff -dir . -pa 40 $oraasm
"
# </editor-fold>

