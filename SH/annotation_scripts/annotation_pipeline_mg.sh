#Gene annotation
#
#Step 0: prepare the data and tools (see the settings in the file annotation.cofig)
#Working_dir/
#abinitio/
#annotation.config*
#evaluation/
#EVM_PASA/
#noncoding/
#protein/
#reference/
#repeat/
#RNAseq/
#scipio/
#version/
#

################################################################################
######################## step 1: Get repeat  regions ###########################
################################################################################

cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/
for sample in cur ora; do
    cd ${cwd}/${sample}/repeat
    mkdir RepeatMasker
    cd RepeatMasker
    ref=../../reference/${sample}.genome.v1.fasta
    bsub -q bigmem -n 40 -R "span[hosts=1] rusage[mem=30000]" -M 35000 -oo rm.log -eo rm.err "
        rm -r RM_*
        BuildDatabase -name $sample -engine ncbi $ref 
#        /srv/netscratch/dep_mercier/grp_schneeberger/bin/RepeatModeler/RepeatModeler-2.0.1/BuildDatabase -name $sample -engine rmblast $ref
        RepeatModeler -engine ncbi \
          -pa 40 \
          -database $sample \
          -LTRStruct
        cp RM_*/consensi.fa.classified consensi.fa.classified
        RepeatMasker -lib ./consensi.fa.classified -gff -dir . -pa 40 $ref
    "
done

# When I try to run RepeatModeler it resulted in too many Unknown annotations. To
# overcome this, I think I can use the libraries created by Wen-Biao as input to
# run RepeatMasker.
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/
for sample in cur ora; do
    cd ${cwd}/${sample}/repeat
    mkdir RepeatMasker_jiao_lib
    cd RepeatMasker_jiao_lib
    ref=../../reference/${sample}.genome.v1.fasta
    #    bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=30000]" -M 35000 -oo rm.log -eo rm.err "
    rm -r RM_*
    cp /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/jiao_folder/annotation/v1/Phased_Currot/repeat/RepeatMasker/consensi.fa.classified currot_consensi.fa.classfied
    cp /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/jiao_folder/annotation/v1/Phased_OrangeRed/repeat/RepeatMasker/consensi.fa.classified orangered_consensi.fa.classfied
    #        cp RM_*/consensi.fa.classified consensi.fa.classified
    cat currot_consensi.fa.classfied orangered_consensi.fa.classfied >consensi.fa.classified
    RepeatMasker -lib ./consensi.fa.classified -gff -dir . -pa 40 $ref &
    #    "
done

################################################################################
############################ step 2: run pipeline ##############################
################################################################################
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur
nohup python ../../scripts/evm.pasa.integrate.pipeline.py -f ./annotation.config &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora
nohup python ../../scripts/evm.pasa.integrate.pipeline.py -f ./annotation.config &

################################################################################
################### step 3 : Improve annotation using PASA #####################
################################################################################

# Get transcripts using genome-guided trinity
for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/trinity
    for r in C D abps 4031; do
        bsub -q bigmem -n 20 -R "span[hosts=1] rusage[mem=110000]" -M 120000 -oo ${r}.log -eo ${r}.err "
      rm -r ${r}_trinity
      Trinity --genome_guided_bam ../RNAseq/${r}.srt.bam \
         --genome_guided_max_intron 100000 \
         --max_memory 100G --CPU 20 \
         --output ${r}_trinity \
         --full_cleanup
      # seqclean abps_trinity.fasta -c 10
    "
    done
done

# Get transcripts using denovo trinity
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/data/RNAseq/trinity_DN
for r in C D abps; do
    bsub -q bigmem -n 40 -R "span[hosts=1] rusage[mem=110000]" -M 120000 -oo ${r}.log -eo ${r}.err "
    rm -r ${r}_trinity_DN
    Trinity --seqType fq \
      --left ../${r}_1.fastq.gz --right ../${r}_2.fastq.gz \
      --CPU 40 --max_memory 100G \
      --output ${r}_trinity_DN \
      --full_cleanup
  "
done

# Only use our Rojo Passion sequenced reads (library: 4031) for running PASA
# Get transcripts using denovo trinity all_reads
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/data/RNAseq/trinity_DN
bsub -q bigmem -n 40 -R "rusage[mem=40000]" -M 50000 -oo 4031_DN.log -eo 4031_DN.err "
  Trinity --seqType fq \
    --single ../4031.fastq.gz \
    --CPU 40 --max_memory 35G \
    --output 4031_trinity_DN \
    --full_cleanup
"

# The manual curations scripts seems to optimised for the EVM output GFF.
# So, I use that for the gene-model correction by manual curation, and then can
# RUN PASA at the end to add UTR and isoform information.

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/data/RNAseq/trinity_DN/
for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/EVM_PASA/PASA

    # Get accessions
    /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/accession_extractor.pl <${indir}4031_trinity_DN.Trinity.fasta >tdn.accs

    # Concatenate the DN and GG annotations
    cat ${indir}4031_trinity_DN.Trinity.fasta ../../trinity/4031_trinity/Trinity-GG.fasta >transcripts.fasta
    seqclean transcripts.fasta -c 1 -n 10000 -o transcripts.fasta.clean

    # Delete the reference gmap index if GMAP version was changed

    rm -r __pasa_${sample}_pasa.sqlite_SQLite_chkpts pasa_run.log.dir transcripts.fasta.clean.fai 11.ooc blat_out_dir __pasa_cur_pasa.sqlite_SQLite_chkpts.cmds_log blat.spliced_alignments.gff3 pasa.log gmap.spliced_alignments.gff3 ${sample}_pasa.sqlite
    # Run PASA assembly pipeline
    nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AlignAssembly.cfg -C -R -g ../../reference/${sample}.genome.v1.fasta --ALIGNERS blat,gmap --TRANSDECODER --CPU 10 -T -t transcripts.fasta.clean -u transcripts.fasta --stringent_alignment_overlap 30.0 |& tee pasa.log &

    nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/scripts/build_comprehensive_transcriptome.dbi -c pasa.AlignAssembly.cfg -t transcripts.fasta.clean &

    nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/${sample}.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots ../evm.all.gff3 &
done

# PASA run 2
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/PASA
nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/cur.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots cur_pasa.sqlite.gene_structures_post_PASA_updates.19892.gff3 &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/PASA
nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/ora.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots ora_pasa.sqlite.gene_structures_post_PASA_updates.19893.gff3 &

# PASA run 3
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/PASA
nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/cur.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots cur_pasa.sqlite.gene_structures_post_PASA_updates.3995.gff3 &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/PASA
nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/ora.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots ora_pasa.sqlite.gene_structures_post_PASA_updates.4008.gff3 &

# Clean up annotation file, make gene IDs and generate common accompanying files

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/PASA
/srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_simple_acc_resetter.pl cur_pasa.sqlite.gene_structures_post_PASA_updates.17477.gff3 ../../reference/cur.genome.v1.fasta >cur.PASA_out.gff3 &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/PASA
/srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_simple_acc_resetter.pl ora_pasa.sqlite.gene_structures_post_PASA_updates.17490.gff3 ../../reference/ora.genome.v1.fasta >ora.PASA_out.gff3 &

for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/EVM_PASA/PASA

    /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.PASA_out.gff3 ../../reference/${sample}.genome.v1.fasta prot >${sample}.PASA_out.proteins.fasta &

    /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.PASA_out.gff3 ../../reference/${sample}.genome.v1.fasta CDS >${sample}.PASA_out.CDS.fasta &

    /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.PASA_out.gff3 ../../reference/${sample}.genome.v1.fasta cDNA >${sample}.PASA_out.cDNA.fasta &

    /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.PASA_out.gff3 ../../reference/${sample}.genome.v1.fasta gene >${sample}.PASA_out.gene.fasta &

    /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_to_gtf_format.pl ${sample}.PASA_out.gff3 ../../reference/${sample}.genome.v1.fasta >${sample}.PASA_out.gtf &
done

################################################################################
###################### step 4 : Get TE based annotations #######################
################################################################################

for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/repeat/RepeatMasker
    # perl ../../../../scripts/repeat.classfied.gff3.pl \
    # ./${sample}.genome.v1.fasta.out.gff \
    # ./${sample}.genome.v1.fasta.out \
    # ${sample}.genome.v1.fasta.ann.gff3 \
    # repeat.ann.stats
    # egrep -v 'Low|Simple|RNA|other|Satellite' ${sample}.genome.v1.fasta.ann.gff3 | cut -f 1,4,5,9 > ${sample}.genome.TE.bed

    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/repeat/TErelated
    # nohup perl ../../../../scripts/remove.TErelated.genes.pl \
    # ../../EVM_PASA/evm.annotation.protein.fasta \
    # ../../EVM_PASA/evm.annotation.gene.fasta \
    # ../RepeatMasker/${sample}.genome.TE.bed \
    # ../../EVM_PASA/evm.all.gff3 ./ &

    # /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/scripts/select_first_isoform_from_gff.py annotation.genes.gff annotation.genes.iso1.gff
    {
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl annotation.genes.gff ../../reference/${sample}.genome.v1.fasta gene >annotation.genes.gff.gene.fasta
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl annotation.genes.gff ../../reference/${sample}.genome.v1.fasta prot >annotation.genes.gff.prot.fasta
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl annotation.genes.gff ../../reference/${sample}.genome.v1.fasta CDS >annotation.genes.gff.cds.fasta
        sed -i "s/>evm.model./>/g" annotation.genes.gff.gene.fasta
        sed -i "s/>evm.model./>/g" annotation.genes.gff.prot.fasta
        sed -i "s/>evm.model./>/g" annotation.genes.gff.cds.fasta
    } &
    # /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/scripts/get_cds_bed_from_gff.py  annotation.genes.iso1.gff ${sample}.cds.bed ${sample}
done

################################################################################
### step 4b : Get BUSCO stats for the annotaions before doing manual curation ###
################################################################################
for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation
    mkdir -p BUSCO/run1
    cd BUSCO/run1
    # hometools gfftrans ../../../repeat/TErelated/annotation.genes.gff ../../../reference/${sample}.genome.v1.fasta
    # ln -s ../../../repeat/TErelated/annotation.genes.iso1.gff.prot.fasta .
    ln -s ../../../repeat/TErelated/annotation.genes.gff.prot.fasta .
    source activate BUSCO
    nohup /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/BUSCO/bin/busco \
        -m prot \
        -i annotation.genes.gff.prot.fasta \
        -o ${sample}_run1_nonTE_gff \
        -l eudicots_odb10 \
        -f \
        -c 30 &
    conda deactivate
done

################################################################################
####################### step 5 : Evaluate and Update ###########################
# find wrong annotation (mis-merge, mis-split, mis-exon, missing) and update them
################################################################################

# STEP:5a Evaluate the Cur/Ora annotation using the Ora/Cur annotation as the reference
# Blast Currot and OrangeRed gene annotations against the OrangeRed and Currot assemblies

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation
mkdir blast
cd blast
# makeblastdb -dbtype nucl -input_type fasta -in ../../reference/cur.genome.v1.fasta
nohup perl ../../../../scripts/assembly.evaluation.by.blastn.pl ../../reference/cur.genome.v1.fasta ../../../ora/repeat/TErelated/annotation.genes.gff.gene.fasta ./ ora.blastn.cur &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/evaluation
mkdir blast
cd blast
# makeblastdb -dbtype nucl -input_type fasta -in ../../reference/ora.genome.v1.fasta
nohup perl ../../../../scripts/assembly.evaluation.by.blastn.pl ../../reference/ora.genome.v1.fasta ../../../cur/repeat/TErelated/annotation.genes.gff.gene.fasta ./ cur.blastn.ora &

# Run Orthofinder between the Currot and OrangeRed Protein Sequences
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation
mkdir orthofinder
cd orthofinder
mkdir prot_seq
ln -sf ../../../repeat/TErelated/annotation.genes.gff.prot.fasta prot_seq/cur.prot.fasta
ln -sf ../../../../ora/repeat/TErelated/annotation.genes.gff.prot.fasta prot_seq/ora.prot.fasta
nohup orthofinder -t 40 -a 40 -f prot_seq &

# Get files for performing annotation comparison
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation
mkdir mis
cd mis
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/scripts/get_cds_bed_from_gff.py ../../repeat/TErelated/annotation.genes.gff cur.gene.bed CUR
sed -i 's/evm.model.//g' cur.gene.bed
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/scripts/get_cds_bed_from_gff.py ../../../ora/repeat/TErelated/annotation.genes.gff ora.gene.bed ORA
sed -i 's/evm.model.//g' ora.gene.bed
ln -sf ../orthofinder/prot_seq/OrthoFinder/Results_Apr07/Orthogroups/Orthogroups.txt .
ln -sf ../blast/ora.blastn.cur.besthit.out .
ln -sf ../../../ora/evaluation/blast/cur.blastn.ora.besthit.out .
ln -sf ../../repeat/TErelated/annotation.genes.gff.prot.fasta cur.prot.fa
ln -sf ../../../ora/repeat/TErelated/annotation.genes.gff.prot.fasta ora.prot.fa
zcat ../orthofinder/prot_seq/OrthoFinder/Results_Apr07/WorkingDirectory/Blast0_1.txt.gz | awk '{if ($3>60) print}' | sort -k1,1 -k2,2 -k7,7n >cur.ora.blastp.i60.srt.txt
zcat ../orthofinder/prot_seq/OrthoFinder/Results_Apr07/WorkingDirectory/Blast1_0.txt.gz | awk '{if ($3>60) print}' | sort -k1,1 -k2,2 -k7,7n >ora.cur.blastp.i60.srt.txt
../../../../scripts/match_blastp_ids.py cur.ora.blastp.i60.srt.txt ../orthofinder/prot_seq/OrthoFinder/Results_Apr07/WorkingDirectory/SequenceIDs.txt cur.ora.blastp.i60.srt.renamed.txt
../../../../scripts/match_blastp_ids.py ora.cur.blastp.i60.srt.txt ../orthofinder/prot_seq/OrthoFinder/Results_Apr07/WorkingDirectory/SequenceIDs.txt ora.cur.blastp.i60.srt.renamed.txt
awk '{if ($3=="transcript") print}' ../../RNAseq/stringtie.merged.gtf >cur.rna.transcript.gtf
awk '{if ($3=="transcript") print}' ../../../ora/RNAseq/stringtie.merged.gtf >ora.rna.transcript.gtf
awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' ora.blastn.cur.besthit.out | sort -k1,1 -k2,2n -k3,3n >ora.blastn.cur.besthit.out2
awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' cur.blastn.ora.besthit.out | sort -k1,1 -k2,2n -k3,3n >cur.blastn.ora.besthit.out2

# Evaluate mis-assemblies
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation/mis
nohup python -u ../../../../scripts/annotation.evaluate.find-mis.py -g Orthogroups.txt -n ora.blastn.cur.besthit.out2 -p cur.ora.blastp.i60.srt.renamed.txt -x ora.prot.fa -y cur.prot.fa -s ora.gene.bed -q cur.gene.bed -m cur.blastn.ora.besthit.out2 -r cur.rna.transcript.gtf -o ./cur >cur.log &
nohup python -u ../../../../scripts/annotation.evaluate.find-mis.py -g Orthogroups.txt -n cur.blastn.ora.besthit.out2 -p ora.cur.blastp.i60.srt.renamed.txt -x cur.prot.fa -y ora.prot.fa -s cur.gene.bed -q ora.gene.bed -m ora.blastn.cur.besthit.out2 -r ora.rna.transcript.gtf -o ./ora >ora.log &

# I tried using PASA before doing these manual curation, but Wen-Biao's script are
# not compatible with the PASA generated GFF. So, now I perform the manual curation
# first and then later I can add UTR information using PASA to the final GFF.
#nohup python -u ../../../../scripts/annotation.evaluate.find-mis_mg.py -g Orthogroups.txt -n ora.blastn.cur.besthit.out2 -p cur.ora.blastp.i60.srt.renamed.txt -x ora.prot.fa -y cur.prot.fa -s ora.gene.bed -q cur.gene.bed -m cur.blastn.ora.besthit.out2 -r cur.rna.transcript.gtf -o ./cur > cur.log &
#nohup python -u ../../../../scripts/annotation.evaluate.find-mis_mg.py -g Orthogroups.txt -n cur.blastn.ora.besthit.out2 -p ora.cur.blastp.i60.srt.renamed.txt -x cur.prot.fa -y ora.prot.fa -s cur.gene.bed -q ora.gene.bed -m ora.blastn.cur.besthit.out2 -r ora.rna.transcript.gtf -o ./ora > ora.log

# STEP:5b Update the Cur/Ora annotation using the Ora/Cur annotation as the reference
declare -A sample_ref=(['cur']='ora' ['ora']='cur')
for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/
    mkdir -p scipio/run1/split
    cd scipio/run1
    perl ../../../../../scripts/util/split.fa.N.pl ../../../../${sample_ref[${sample}]}/repeat/TErelated/annotation.genes.gff.prot.fasta split ${sample_ref[${sample}]}.split 500
    nohup python ../../../../../scripts/run.scipio.py -i ./split/ -o ./split/ -r ../../../reference/${sample}.genome.v1.fasta -t 40 >run.log &
done

for sample in cur ora; do
    {
        cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/
        mkdir -p update/run1
        cd update/run1

        ln -sf ../../../../cur/evaluation/mis/${sample}/alt-genes.to.be.updated.added.srt.txt .
        ln -sf ../../../repeat/TErelated/annotation.genes.TE.gff annotation.preV1.gff
        perl ../../../../../scripts/util/gff.sort.pl ../../../abinitio/abinitio.4evm.gff SNAP SNAP.ann.gff
        perl ../../../../../scripts/util/gff.sort.pl ../../../abinitio/abinitio.4evm.gff GlimmerHMM GlimmerHMM.ann.gff
        perl ../../../../../scripts/util/gff.sort.pl ../../../abinitio/abinitio.4evm.gff Augustus Augustus.ann.gff
        cut -f1,2,3,4 ../../../../cur/evaluation/mis/${sample_ref[${sample}]}.blastn.${sample}.besthit.out2 >ref.blastn.bed
        ln -sf ../../../reference/${sample}.genome.v1.fasta .
        ln -sf ../../../../cur/evaluation/mis/${sample_ref[${sample}]}.prot.fa .
        ln -sf ../../../../cur/evaluation/mis/${sample_ref[${sample}]}.gene.bed .
        cat ../../scipio/run1/split/${sample_ref[${sample}]}.split.*.fa.scipio.gff >scipio.gff

        nohup python -u ../../../../../scripts/annotation.correct.update.py \
            -u alt-genes.to.be.updated.added.srt.txt \
            -g annotation.preV1.gff \
            -o . \
            -s scipio.gff \
            -a Augustus.ann.gff \
            -n SNAP.ann.gff \
            -l GlimmerHMM.ann.gff \
            -b ref.blastn.bed \
            -f ${sample}.genome.v1.fasta \
            -p ${sample_ref[${sample}]}.prot.fa \
            -i ${sample_ref[${sample}]}.gene.bed
        python ../../../../../scripts/annotation.gene.ID.update.py -i updated.highConf.gff -o ./ -v v1.0 -a ${sample} -g ${sample}.genome.v1.fasta
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v1.0.gff ${sample}.genome.v1.fasta gene >${sample}.protein-coding.genes.v1.0.gff.gene.fasta
        sed -i 's/\.1//g' ${sample}.protein-coding.genes.v1.0.gff.gene.fasta # Set sequence IDs to gene IDs
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v1.0.gff ${sample}.genome.v1.fasta prot >${sample}.protein-coding.genes.v1.0.gff.prot.fasta
        sed -i 's/\.1//g' ${sample}.protein-coding.genes.v1.0.gff.prot.fasta # Set sequence IDs to gene IDs

    } &
done

## Repeat the evaluation and update step using the outputs of Run1

# Blast Currot and OrangeRed gene annotations against the OrangeRed and Currot assemblies
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation
mkdir blast/run2
cd blast/run2
nohup perl ../../../../../scripts/assembly.evaluation.by.blastn.pl ../../../reference/cur.genome.v1.fasta ../../../../ora/evaluation/update/run1/ora.protein-coding.genes.v1.0.gff.gene.fasta ./ ora.blastn.cur &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/evaluation
mkdir blast/run2
cd blast/run2
nohup perl ../../../../../scripts/assembly.evaluation.by.blastn.pl ../../../reference/ora.genome.v1.fasta ../../../../cur/evaluation/update/run1/cur.protein-coding.genes.v1.0.gff.gene.fasta ./ cur.blastn.ora &

# Run Orthofinder between the Currot and OrangeRed Protein Sequences
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation/orthofinder
mkdir run2
ln -sf ../../update/run1/cur.protein-coding.genes.v1.0.gff.prot.fasta run2/.
ln -sf ../../../../ora/evaluation/update/run1/ora.protein-coding.genes.v1.0.gff.prot.fasta run2/.
nohup orthofinder -t 40 -a 40 -f run2 &

# Get files for performing annotation comparison
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation/mis/run2
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/scripts/get_cds_bed_from_gff.py ../../update/run1/cur.protein-coding.genes.v1.0.gff cur.gene.bed CUR
sed -i 's/\.1//g' cur.gene.bed
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/scripts/get_cds_bed_from_gff.py ../../../../ora/evaluation/update/run1/ora.protein-coding.genes.v1.0.gff ora.gene.bed ORA
sed -i 's/\.1//g' ora.gene.bed
ln -sf ../../orthofinder/run2/OrthoFinder/Results_Apr12/Orthogroups/Orthogroups.txt .
ln -sf ../../blast/run2/ora.blastn.cur.besthit.out .
ln -sf ../../../../ora/evaluation/blast/run2/cur.blastn.ora.besthit.out .
ln -sf ../../update/run1/cur.protein-coding.genes.v1.0.gff.prot.fasta cur.prot.fa
ln -sf ../../../../ora/evaluation/update/run1/ora.protein-coding.genes.v1.0.gff.prot.fasta ora.prot.fa
zcat ../../orthofinder/run2/OrthoFinder/Results_Apr12/WorkingDirectory/Blast0_1.txt.gz | awk '{if ($3>60) print}' | sort -k1,1 -k2,2 -k7,7n >cur.ora.blastp.i60.srt.txt
zcat ../../orthofinder/run2/OrthoFinder/Results_Apr12/WorkingDirectory/Blast1_0.txt.gz | awk '{if ($3>60) print}' | sort -k1,1 -k2,2 -k7,7n >ora.cur.blastp.i60.srt.txt
../../../../../scripts/match_blastp_ids.py cur.ora.blastp.i60.srt.txt ../../orthofinder/run2/OrthoFinder/Results_Apr12/WorkingDirectory/SequenceIDs.txt cur.ora.blastp.i60.srt.renamed.txt
../../../../../scripts/match_blastp_ids.py ora.cur.blastp.i60.srt.txt ../../orthofinder/run2/OrthoFinder/Results_Apr12/WorkingDirectory/SequenceIDs.txt ora.cur.blastp.i60.srt.renamed.txt
awk '{if ($3=="transcript") print}' ../../../RNAseq/stringtie.merged.gtf >cur.rna.transcript.gtf
awk '{if ($3=="transcript") print}' ../../../../ora/RNAseq/stringtie.merged.gtf >ora.rna.transcript.gtf
awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' ora.blastn.cur.besthit.out | sort -k1,1 -k2,2n -k3,3n >ora.blastn.cur.besthit.out2
awk '{print $2"\t"$7"\t"$8"\t"$1"\t"$9"\t"$10}' cur.blastn.ora.besthit.out | sort -k1,1 -k2,2n -k3,3n >cur.blastn.ora.besthit.out2

# Evaluate mis-assemblies
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/evaluation/mis/run2
nohup python -u ../../../../../scripts/annotation.evaluate.find-mis.py -g Orthogroups.txt -n ora.blastn.cur.besthit.out2 -p cur.ora.blastp.i60.srt.renamed.txt -x ora.prot.fa -y cur.prot.fa -s ora.gene.bed -q cur.gene.bed -m cur.blastn.ora.besthit.out2 -r cur.rna.transcript.gtf -o ./cur >cur.log &
nohup python -u ../../../../../scripts/annotation.evaluate.find-mis.py -g Orthogroups.txt -n cur.blastn.ora.besthit.out2 -p ora.cur.blastp.i60.srt.renamed.txt -x cur.prot.fa -y ora.prot.fa -s cur.gene.bed -q ora.gene.bed -m ora.blastn.cur.besthit.out2 -r ora.rna.transcript.gtf -o ./ora >ora.log &

# STEP:5b Run2: Update the Cur/Ora annotation using the Ora/Cur annotation as the reference
declare -A sample_ref=(['cur']='ora' ['ora']='cur')
for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/
    mkdir -p scipio/run2/split
    cd scipio/run2
    perl ../../../../../scripts/util/split.fa.N.pl ../../../../${sample_ref[${sample}]}/evaluation/update/run1/${sample_ref[${sample}]}.protein-coding.genes.v1.0.gff.prot.fasta split ${sample_ref[${sample}]}.split 500
    nohup python ../../../../../scripts/run.scipio.py -i ./split/ -o ./split/ -r ../../../reference/${sample}.genome.v1.fasta -t 40 >run.log &
done

for sample in cur ora; do
    {
        cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/
        mkdir -p update/run2
        cd update/run2

        ln -sf ../../../../cur/evaluation/mis/run2/${sample}/alt-genes.to.be.updated.added.srt.txt .
        ln -sf ../run1/${sample}.protein-coding.genes.v1.0.gff annotation.preV1.gff
        perl ../../../../../scripts/util/gff.sort.pl ../../../abinitio/abinitio.4evm.gff SNAP SNAP.ann.gff
        perl ../../../../../scripts/util/gff.sort.pl ../../../abinitio/abinitio.4evm.gff GlimmerHMM GlimmerHMM.ann.gff
        perl ../../../../../scripts/util/gff.sort.pl ../../../abinitio/abinitio.4evm.gff Augustus Augustus.ann.gff
        cut -f1,2,3,4 ../../../../cur/evaluation/mis/run2/${sample_ref[${sample}]}.blastn.${sample}.besthit.out2 >ref.blastn.bed
        ln -sf ../../../reference/${sample}.genome.v1.fasta .
        ln -sf ../../../../cur/evaluation/mis/run2/${sample_ref[${sample}]}.prot.fa .
        ln -sf ../../../../cur/evaluation/mis/run2/${sample_ref[${sample}]}.gene.bed .
        cat ../../scipio/run2/split/${sample_ref[${sample}]}.split.*.fa.scipio.gff >scipio.gff
        sed -i 's/%2D/-/g' scipio.gff

        nohup python -u ../../../../../scripts/annotation.correct.update.py \
            -u alt-genes.to.be.updated.added.srt.txt \
            -g annotation.preV1.gff \
            -o . \
            -s scipio.gff \
            -a Augustus.ann.gff \
            -n SNAP.ann.gff \
            -l GlimmerHMM.ann.gff \
            -b ref.blastn.bed \
            -f ${sample}.genome.v1.fasta \
            -p ${sample_ref[${sample}]}.prot.fa \
            -i ${sample_ref[${sample}]}.gene.bed

        cat ../../../repeat/TErelated/annotation.genes.TE.gff |
            sed 's/evm\.model\.//g' |
            sed 's/evm\.TU\.//g' |
            sed 's/EVM%20prediction%20//g' \
                >annotation.preV1.TE.gff

        nohup perl ../../../../../scripts/util/get.TE.gff.pl annotation.preV1.TE.gff TE.gff >nohup.out
        cat TE.gff updated.highConf.gff >annotation.preV2.TE.gff

        nohup perl ../../../../../scripts/util/gff.sort.pl ./annotation.preV2.TE.gff tool annotation.preV2.TE.srt.gff >nohup.out

        python ../../../../../scripts/annotation.gene.ID.update.py -i annotation.preV2.TE.srt.gff -o ./ -v v2.0 -a ${sample} -g ${sample}.genome.v1.fasta

        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v2.0.gff ${sample}.genome.v1.fasta gene >${sample}.protein-coding.genes.v2.0.gff.gene.fasta
        sed -i 's/\.1//g' ${sample}.protein-coding.genes.v2.0.gff.gene.fasta # Set sequence IDs to gene IDs
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v2.0.gff ${sample}.genome.v1.fasta prot >${sample}.protein-coding.genes.v2.0.gff.prot.fasta
        sed -i 's/\.1//g' ${sample}.protein-coding.genes.v2.0.gff.prot.fasta # Set sequence IDs to gene IDs
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v2.0.gff ${sample}.genome.v1.fasta CDS >${sample}.protein-coding.genes.v2.0.gff.cds.fasta
        sed -i 's/\.1//g' ${sample}.protein-coding.genes.v2.0.gff.cds.fasta # Set sequence IDs to gene IDs
    } &
done

################################################################################
################### step 6 : Remove invalid gene models ########################
# inframe stop, CDS length not 3 times
################################################################################
for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/validate
    {
        nohup perl ../../../../scripts/gene.model.validate.pl \
            ../update/run2/${sample}.protein-coding.genes.v2.0.gff \
            ../update/run2/${sample}.protein-coding.genes.v2.0.gff.gene.fasta \
            ../update/run2/${sample}.protein-coding.genes.v2.0.gff.prot.fasta \
            ../update/run2/${sample}.protein-coding.genes.v2.0.gff.cds.fasta \
            ./${sample}.protein-coding.genes.v2.1.gff >nohup.out

        grep -i ^${sample}.g ${sample}.protein-coding.genes.v2.1.gff >${sample}.protein-coding.genes.v2.1.only_chr.gff
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v2.1.only_chr.gff ../../reference/${sample}.genome.v1.fasta gene | awk '{if (/>/) print $1;else print}' >${sample}.v2.1.gene.fasta
        sed -i 's/\.1//g' ${sample}.v2.1.gene.fasta
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v2.1.only_chr.gff ../../reference/${sample}.genome.v1.fasta prot | awk '{if (/>/) print $1;else print}' >${sample}.v2.1.prot.fasta
        sed -i 's/\.1//g' ${sample}.v2.1.prot.fasta
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.protein-coding.genes.v2.1.only_chr.gff ../../reference/${sample}.genome.v1.fasta CDS | awk '{if (/>/) print $1;else print}' >${sample}.v2.1.cds.fasta
        sed -i 's/\.1//g' ${sample}.v2.1.cds.fasta
    } &
done

################################################################################
############## Step 7 : Fix differeneces between haplotypes ####################
# inframe stop, CDS length not 3 times
################################################################################
# Run SyRI on the assemblies
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run
hometools getchr --chrs CUR1G CUR2G CUR3G CUR4G CUR5G CUR6G CUR7G CUR8G -o cur.genome.v1.filtered.fasta ../../cur/reference/cur.genome.v1.fasta
hometools getchr --chrs ORA1G ORA2G ORA3G ORA4G ORA5G ORA6G ORA7G ORA8G -o ora.genome.v1.filtered.fasta ../../ora/reference/ora.genome.v1.fasta

minimap2 -ax asm5 -t 50 --eqx cur.genome.v1.filtered.fasta ora.genome.v1.filtered.fasta |
    samtools sort -O BAM -@ 50 - \
        >out.bam
samtools index -@ 50 out.bam
nohup python3 /srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri/syri/bin/syri \
    -c out.bam \
    -r cur.genome.v1.filtered.fasta \
    -q ora.genome.v1.filtered.fasta \
    -k -F B --nc 8 &

## Run syri with the --all parameter to also find SVs in the duplication
syri \
    -c out.bam \
    -r cur.genome.v1.filtered.fasta \
    -q ora.genome.v1.filtered.fasta \
    -k -F B --nc 8 --all --prefix all_sv. &


python3 $PATH_TO_PLOTSR syri.out refgenome qrygenome -H 8 -W 5

# Run Orthofinder
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/
ln -s ../../cur/evaluation/validate/cur.v2.1.prot.fasta orthofinder/.
ln -s ../../ora/evaluation/validate/ora.v2.1.prot.fasta orthofinder/.
nohup orthofinder -t 40 -a 40 -f ./orthofinder &

# Run Blastn
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/blastn/
mkdir ora_gene_cur_ref
cd ora_gene_cur_ref
makeblastdb -dbtype nucl -input_type fasta -in ../../syri_run/cur.genome.v1.filtered.fasta
nohup perl ../../../../scripts/assembly.evaluation.by.blastn.pl ../../syri_run/cur.genome.v1.filtered.fasta ../../../ora/evaluation/validate/ora.v2.1.gene.fasta ./ ora.gene.blastn.cur &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/blastn/
mkdir cur_gene_ora_ref
cd cur_gene_ora_ref
makeblastdb -dbtype nucl -input_type fasta -in ../../syri_run/ora.genome.v1.filtered.fasta
nohup perl ../../../../scripts/assembly.evaluation.by.blastn.pl ../../syri_run/ora.genome.v1.filtered.fasta ../../../cur/evaluation/validate/cur.v2.1.gene.fasta ./ cur.gene.blastn.ora &

# Get Haplotype differences in gene annotation
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/run1
ln -sf ../syri_run/cur.genome.v1.filtered.fasta cur.genome
ln -sf ../syri_run/ora.genome.v1.filtered.fasta ora.genome
nohup python ../../../scripts/haplotype.gene.diff.py \
    -g ../orthofinder/OrthoFinder/Results_Apr15/Orthogroups/Orthogroups.txt \
    -s ../syri_run/syri.out \
    -r ../../cur/evaluation/validate/cur.protein-coding.genes.v2.1.only_chr.gff \
    -a ../../ora/evaluation/validate/ora.protein-coding.genes.v2.1.only_chr.gff \
    -m ../syri_run/cur.chr.size \
    -n ../syri_run/ora.chr.size \
    -x CUR -y ORA \
    -j ../blastn/cur_gene_ora_ref/cur.gene.blastn.ora.besthit.out \
    -k ../blastn/ora_gene_cur_ref/ora.gene.blastn.cur.besthit.out \
    -o ./ \
    >nohup.out &

# Step 7B: Update2 : add the unannotated genes (based on the other haplotype specific genes, their blastn aginst the haplotype assembly and original abinitio annotation). V1.2
declare -A sample_ref=(['cur']='ora' ['ora']='cur')
declare -A other_ref=(['cur']='alt' ['ora']='ref')

for sample in cur ora; do
    {
        cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/update2
        ln -sf ../validate/${sample}.protein-coding.genes.v2.1.only_chr.gff .
        ln -sf ../../../haplodiff/run1/${other_ref[${sample}]}.spGene.noDel.blastn.srt.bed .
        ln -sf ../../../${sample_ref[${sample}]}/evaluation/validate/${sample_ref[${sample}]}.v2.1.prot.fasta .

        awk '{if ($3=="gene") print}' ../update/run2/SNAP.ann.gff | cut -f 1,4,5,7,9 | sort -k1,1 -k2,2n -k3,3n >ann/snap.bed
        awk '{if ($3=="gene") print}' ../update/run2/GlimmerHMM.ann.gff | cut -f 1,4,5,7,9 | sort -k1,1 -k2,2n -k3,3n >ann/glimmerhmm.bed
        big_id=$(echo $sample | tr '[:lower:]' '[:upper:]')
        cat ../update/run2/Augustus.ann.gff | grep "^$big_id" >augustus.gff
        gffread --keep-genes -F --keep-exon-attrs augustus.gff >augustus2.gff
        awk -F $'\t' '
        BEGIN {OFS="\t"}
        {if($3=="gene") {
            print $0;
            exon_count=0
        }
        else if ($3=="mRNA") print $0;
        else if ($3=="CDS") {
            split($9,ID,"=");
            print $1,$2,$3,$4,$5,$6,$7,$8,"ID=cds."ID[2]";"$9
        }
        else if ($3=="exon") {
            exon_count+=1;
            split($9,ID,"=");
            print $1,$2,$3,$4,$5,$6,$7,$8,"ID="ID[2]".exon"exon_count";"$9
        }
    }' augustus2.gff >augustus3.gff
        perl ../../../../scripts/util/augustus.gff.fix.pl augustus.gff augustus4.gff
        cp augustus4.gff ann/augustus.gff
        rm augustus.gff augustus2.gff augustus3.gff augustus4.gff
        awk '{if ($3=="gene") print}' ann/augustus.gff | cut -f 1,4,5,7,9 >ann/augustus.bed
        cp ../update/run2/SNAP.ann.gff ./ann/snap.gff
        cp ../update/run2/GlimmerHMM.ann.gff ./ann/glimmerhmm.gff
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ann/augustus.gff ../update/run2/${sample}.genome.v1.fasta | awk '{if (/>/) print $1;else print}' >ann/augustus.prot.fa
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ann/snap.gff ../update/run2/${sample}.genome.v1.fasta | awk '{if (/>/) print $1;else print}' >ann/snap.prot.fa
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ann/glimmerhmm.gff ../update/run2/${sample}.genome.v1.fasta | awk '{if (/>/) print $1;else print}' >ann/glimmerhmm.prot.fa
        sed -i 's/model/gene/g' ann/glimmerhmm.prot.fa
        sed -i 's/model/gene/g' ann/snap.prot.fa
        sed -i 's/model/gene/g' ann/augustus.prot.fa
        nohup perl ../../../../scripts/add.false.spGenes.pl ${other_ref[${sample}]}.spGene.noDel.blastn.srt.bed ./${sample}.protein-coding.genes.v2.1.only_chr.gff ${sample_ref[${sample}]}.v2.1.prot.fasta ./ann ./run1/ ${sample} >run1.log
    } &
done

for sample in cur ora; do
    {
        cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/update2/run1
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl ${sample}.new.protein-coding.genes.gff ../../update/run2/${sample}.genome.v1.fasta prot | awk '{if (/>/) print $1;else print}' >${sample}.new.protein-coding.genes.prot.fasta
    } &
done

################################################################################
################# Step 8: Remove low confidence annotations ####################
################################################################################

# Run orthofinder using the A.Thal and Cur/Ora proteins found using update2
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/orthofinder/athal_cur_ora/run1
# Use Araport 11 proteins
ln -s /srv/biodata/dep_mercier/grp_schneeberger/data/Athal/Araport11/Araport11_genes.201606.pep.longest_transcript.fasta ath.fa
ln -s ../../../../cur/evaluation/update2/run1/cur.new.protein-coding.genes.prot.fasta cur.prot.fa
ln -s ../../../../ora/evaluation/update2/run1/ora.new.protein-coding.genes.prot.fasta ora.prot.fa
cd ..
nohup orthofinder -t 40 -a 40 -f ./run1 &
gff3_to_proteins=/srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl

declare -A ref_type=(['cur']='ref' ['ora']='alt')
declare -A sample_ref=(['cur']='ora' ['ora']='cur')
for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/rmlowconf
    ln -sf ../../../haplodiff/orthofinder/athal_cur_ora/run1/OrthoFinder/Results_Apr22/Orthogroups/Orthogroups.txt ath.cur.ora.groups.txt
    cut -f 6-13 ../../../haplodiff/run1/${ref_type[${sample}]}.spGene.noDel.noLoF.blastn.srt.bed >${ref_type[${sample}]}.spGenes.noDel.noLoF.srt.bed
    big_id=$(echo $sample | tr '[:lower:]' '[:upper:]')
    cut -f 1,4,5 ../../RNAseq/rnaseq.4evm.gff |
        grep -v '#' | grep '^$big_id' | sort -k1,1 -k2,2n -k3,3n |
        bedtools merge -i - >${sample}.rna.transcript.merged.bed
    ln -sf ../update2/run1/${sample}.new.protein-coding.genes.gff .
    ln -sf ../../../${sample_ref[${sample}]}/evaluation/update2/run1/ungr.gene.can.be.added.txt .
    nohup perl ../../../../scripts/remove.lowConf.gene.model.pl ungr.gene.can.be.added.txt ${ref_type[${sample}]}.spGenes.noDel.noLoF.srt.bed ath.cur.ora.groups.txt ${sample}.rna.transcript.merged.bed ${sample}.new.protein-coding.genes.gff ${sample}.new.protein-coding.genes.rmlow.gff v3.rm.id.run1.txt >run1.log
    $gff3_to_proteins ${sample}.new.protein-coding.genes.rmlow.gff ../update/run2/${sample}.genome.v1.fasta prot | awk '{if (/>/) print $1;else print}' >${sample}.prot.fasta
    $gff3_to_proteins ${sample}.new.protein-coding.genes.rmlow.gff ../update/run2/${sample}.genome.v1.fasta gene | awk '{if (/>/) print $1;else print}' >${sample}.gene.fasta
    $gff3_to_proteins ${sample}.new.protein-coding.genes.rmlow.gff ../update/run2/${sample}.genome.v1.fasta CDS | awk '{if (/>/) print $1;else print}' >${sample}.cds.fasta
    nohup perl ../../../../scripts/gene.model.validate.pl ${sample}.new.protein-coding.genes.rmlow.gff ${sample}.gene.fasta ${sample}.prot.fasta ${sample}.cds.fasta ${sample}.protein-coding.v2.2.gff >nohup.out &
done

# WenBiao's pipeline had a step to filter out very long genes, but for my run
# there were not any very long (>150k) genes, so that step was not needed.

################################################################################
################### Step 9: Run PASA_pipeline to add UTRs ######################
################################################################################

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/data/RNAseq/trinity_DN/
for sample in cur ora; do
    {
        cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/EVM_PASA/pasa_on_mancur

        # rm -r __pasa_${sample}_pasa.sqlite_SQLite_chkpts pasa_run.log.dir transcripts.fasta.clean.fai 11.ooc blat_out_dir __pasa_cur_pasa.sqlite_SQLite_chkpts.cmds_log blat.spliced_alignments.gff3 pasa.log gmap.spliced_alignments.gff3 gmap.spliced_alignments.gff3.completed ${sample}_pasa.sqlite
        rm -r !(pasa.A*.cfg)

        # Get accessions
        /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/accession_extractor.pl <${indir}4031_trinity_DN.Trinity.fasta >tdn.accs

        # Concatenate the DN and GG annotations
        cat ${indir}4031_trinity_DN.Trinity.fasta ../../trinity/4031_trinity/Trinity-GG.fasta >transcripts.fasta
        seqclean transcripts.fasta -c 15 -n 10000 -o transcripts.fasta.clean

        # Delete the reference gmap index if GMAP version was changed

        # Run PASA assembly pipeline
        nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AlignAssembly.cfg -C -R -g ../../reference/${sample}.genome.v1.fasta --ALIGNERS blat,gmap --TRANSDECODER --CPU 30 -T -t transcripts.fasta.clean -u transcripts.fasta --stringent_alignment_overlap 30.0 |& tee pasa.log

        nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/scripts/build_comprehensive_transcriptome.dbi -c pasa.AlignAssembly.cfg -t transcripts.fasta.clean

        big_id=$(echo $sample | tr '[:lower:]' '[:upper:]')
        cat ../../evaluation/rmlowconf/${sample}.protein-coding.v2.2.gff ../../evaluation/update/run2/TE.gff |
            grep "^${big_id}" >${sample}.protein_coding.TE.gff

        source activate syri3.8
        hometools gffsort ${sample}.protein_coding.TE.gff ${sample}.protein_coding.TE.sorted.gff
        source deactivate

        source activate mgpy2.7
        python ../../../../scripts/annotation.gene.ID.update.py -i ${sample}.protein_coding.TE.sorted.gff -o ./ -v v2.3 -a ${sample} -g ../../reference/${sample}.genome.v1.fasta
        source deactivate

        mv ${sample}.genes.annotation.v2.3.gff ${sample}.genes.annotation.v2.3.gff3

    } &
done

for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/EVM_PASA/pasa_on_mancur
    python ../../../../scripts/add_TE_cds.py ${sample}.genes.annotation.v2.3.gff3 ${sample}.genes.annotation.v2.3.CDS.gff3
    nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/${sample}.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots ${sample}.genes.annotation.v2.3.CDS.gff3 >pasa_annot_update.log &
done

# PASA run 2
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur
nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/cur.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots cur_pasa.sqlite.gene_structures_post_PASA_updates.13304.gff3 &

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/pasa_on_mancur
nohup /srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl -c pasa.AnnotCompare.cfg -g ../../reference/ora.genome.v1.fasta -t transcripts.fasta.clean -A -L --annots ora_pasa.sqlite.gene_structures_post_PASA_updates.13311.gff3 &

# Clean up annotation file, make gene IDs and generate common accompanying files
gff3_to_proteins=/srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_to_proteins.pl

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur
/srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_simple_acc_resetter.pl cur_pasa.sqlite.gene_structures_post_PASA_updates.703.gff3 ../../reference/cur.genome.v1.fasta >cur.PASA_out.gff3

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/pasa_on_mancur
/srv/netscratch/dep_mercier/grp_schneeberger/bin/PASA/PASApipeline.v2.4.1/misc_utilities/gff3_file_simple_acc_resetter.pl ora_pasa.sqlite.gene_structures_post_PASA_updates.718.gff3 ../../reference/ora.genome.v1.fasta >ora.PASA_out.gff3

for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/EVM_PASA/pasa_on_mancur
    big_id=$(echo $sample | tr '[:lower:]' '[:upper:]')
    grep "^${big_id}" ${sample}.PASA_out.gff3 >${sample}.PASA_out.chr.gff3
    hometools gffsort ${sample}.PASA_out.chr.gff3 ${sample}.pasa_out.sort.gff3
    # Mark rRNA at specified rRNA cluster sites (run the script manually)
#    /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/mark_rrna_trna.py

    # Get protein coding genes
    awk '{
        if($3=="gene"){
            if($9~"trans" || $9~"rRNA") trans=1
            else {
                trans=0
                print $0
            }
        }
        else{
            if(trans==0) print $0
        }
    }' ${sample}.pasa_out.sort.rrna.gff3 \
    | sed 's/=%20/=/g' >${sample}.pasa_out.sort.protein_coding.gff3
    # Get rRNA genes
    awk '{
        if($3=="gene"){
            if($9~"rRNA"){
                rrna=1
                print $0
            }
            else rrna=0
        }
        else{
            if(rrna==1) print $0
        }
    }' ${sample}.pasa_out.sort.rrna.gff3 \
    | sed 's/=%20/=/g' > ${sample}.pasa_out.sort.only_rrna.gff3
    # Get TE genes
    awk '{
        if($3=="gene"){
            if($9~"trans"){
                trans=1
                print $0
            }
            else trans=0
        }
        else{
            if(trans==1) print $0
        }
    }' ${sample}.pasa_out.sort.rrna.gff3 \
    | sed 's/=%20/=/g' >${sample}.pasa_out.sort.TE.gff3

    # Add 3UTR to all transcripts that do not have it
    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/add_3UTR_to_gff_files.py ${sample}.pasa_out.sort.protein_coding.gff3 ${sample}.pasa_out.sort.protein_coding.3utr.gff3
    sed -i 's/Name=.*/Name=protein_coding_gene/g' ${sample}.pasa_out.sort.protein_coding.3utr.gff3

    $gff3_to_proteins ${sample}.pasa_out.sort.protein_coding.3utr.gff3 ../../reference/${sample}.genome.v1.fasta prot | awk '{if (/>/) print $1;else print}' >${sample}.pasa_out.prot.fasta &
    $gff3_to_proteins ${sample}.pasa_out.sort.protein_coding.3utr.gff3 ../../reference/${sample}.genome.v1.fasta gene | awk '{if (/>/) print $1;else print}' >${sample}.pasa_out.gene.fasta &
    $gff3_to_proteins ${sample}.pasa_out.sort.protein_coding.3utr.gff3 ../../reference/${sample}.genome.v1.fasta CDS | awk '{if (/>/) print $1;else print}' >${sample}.pasa_out.CDS.fasta &
    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/select_longest_transcripts.py ${sample}.pasa_out.prot.fasta ${sample}.pasa_out.longest_prot.fasta

    # Merge coding and rRNA genes
    cat ${sample}.pasa_out.sort.protein_coding.3utr.gff3 ${sample}.pasa_out.sort.only_rrna.gff3 \
    > ${sample}.pasa_out.3utr.gene_rrna.gff3
    hometools gffsort ${sample}.pasa_out.3utr.gene_rrna.gff3 ${sample}.pasa_out.3utr.gene_rrna.sort.gff3

    # Merge coding, rRNA, and TE genes
    cat ${sample}.pasa_out.sort.protein_coding.3utr.gff3 ${sample}.pasa_out.sort.TE.gff3 ${sample}.pasa_out.sort.only_rrna.gff3 >${sample}.pasa_out.3utr.gff3
    hometools gffsort ${sample}.pasa_out.3utr.gff3 ${sample}.pasa_out.3utr.sort.gff3

done

################################################################################
##################### Step 9: Save version and run BUSCO #######################
################################################################################

for sample in cur ora; do
    cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/${sample}/evaluation/BUSCO
    d=$(date +%d_%b_%Y)
    rm -r run${d}
    mkdir run${d}
    cd run${d}
    cp ../../../EVM_PASA/pasa_on_mancur/${sample}.pasa_out.longest_prot.fasta .
    conda activate busco
    nohup /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/BUSCO/bin/busco \
        -m prot \
        -i ${sample}.pasa_out.longest_prot.fasta \
        -o ${sample}_busco_out \
        -l eudicots_odb10 \
        -f \
        -c 30 &
    conda deactivate
done

################################################################################
############# Step 10: Functional annotation of predicted genes ################
################################################################################

# Currot
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/function
sed 's/\*//g' ../EVM_PASA/pasa_on_mancur/cur.pasa_out.prot.fasta >cur.pasa_out.prot.fasta
/srv/netscratch/dep_mercier/grp_schneeberger/private/manish/toolbox/Parser/FASTA/split_fasta.pl -n 25 -o split_prot_fasta -f cur.pasa_out.prot.fasta
export PATH=/srv/netscratch/dep_mercier/grp_schneeberger/bin/java/jdk-11.0.2/bin:$PATH
fins=$(ls ./split_prot_fasta/cur.pasa_out.prot.split.*fa)
finished=($(ls cur.pasa_out.prot.split.*.fa.tsv))
runfin=()
for fin in ${fins[@]}; do
    name=$(basename $fin)lea
#    echo $name
    if [[ ${finished[*]} =~ ${name}.tsv ]]; then
        continue
    else
        runfin+=($fin)
    fi
done

echo ${runfin[@]} |
    tr -d '\n' |
    xargs -n 1 -P 25 -d ' ' -I {} \
        bash -c '
            name=$(basename {})
            /srv/netscratch/dep_mercier/grp_schneeberger/bin/InterProScan/interproscan-5.48-83.0/interproscan_mg.sh -i {} -f tsv -o ${name}.tsv
    ' -- {}

# Run eggNOG-mapper
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/function
fins=$(ls ./split_prot_fasta/cur.pasa_out.prot.split.*fa)
for fin in ${fins[@]}; do
    name=$(basename $fin)
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -o eggnog.log -e eggnog.err -m "hpc001 hpc002 hpc004 hpc005 hpc006" "
        echo Analysis $name
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python \
        /srv/netscratch/dep_mercier/grp_schneeberger/software/eggnog-mapper/emapper.py \
        -i $fin \
        --itype proteins \
        -o $name \
        --output_dir eggnog_output \
        --temp_dir eggnog_tmp
        echo Finished $name
    "
done

# Merge the output of Eggnog-mapper into a single file
grep 'query' eggnog_output/cur.pasa_out.prot.split.1.fa.emapper.annotations >  cur.pasa_out.prot.fasta.eggnog_mapper.tsv
grep -h mRNA eggnog_output/*.annotations >> cur.pasa_out.prot.fasta.eggnog_mapper.tsv

# Number of functionally annotated genes
cut -f 1 cur.pasa_out.prot.fasta.eggnog_mapper.tsv | cut -d'.' -f1,2 | tail +2| sed 's/mRNA/Gene/g' | sort | uniq | wc -l

# Orangered
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/function
sed 's/\*//g' ../EVM_PASA/pasa_on_mancur/ora.pasa_out.prot.fasta > ora.pasa_out.prot.fasta
/srv/netscratch/dep_mercier/grp_schneeberger/private/manish/toolbox/Parser/FASTA/split_fasta.pl -n 25 -o split_prot_fasta -f ora.pasa_out.prot.fasta
export PATH=/srv/netscratch/dep_mercier/grp_schneeberger/bin/java/jdk-11.0.2/bin:$PATH
fins=$(ls ./split_prot_fasta/ora.pasa_out.prot.split.*fa)
finished=($(ls ora.pasa_out.prot.split.*.fa.tsv))
runfin=()
for fin in ${fins[@]}; do
    name=$(basename $fin)lea
    if [[ ${finished[*]} =~ ${name}.tsv ]]; then
        continue
    else
        runfin+=($fin)
    fi
done

echo ${runfin[@]} |
    tr -d '\n' |
    xargs -n 1 -P 25 -d ' ' -I {} \
        bash -c '
            name=$(basename {})
            /srv/netscratch/dep_mercier/grp_schneeberger/bin/InterProScan/interproscan-5.48-83.0/interproscan_mg.sh -i {} -f tsv -o ${name}.tsv
    ' -- {}

# Run eggNOG-mapper
fins=$(ls ./split_prot_fasta/ora.pasa_out.prot.split.*fa)
for fin in ${fins[@]}; do
    name=$(basename $fin)
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -o eggnog.log -e eggnog.err "
        echo Analysis $name
        /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/python \
        /srv/netscratch/dep_mercier/grp_schneeberger/software/eggnog-mapper/emapper.py \
        -i $fin \
        --itype proteins \
        -o $name \
        --output_dir eggnog_output \
        --temp_dir eggnog_tmp
        echo Finished $name
    "
done

# Merge the output of Eggnog-mapper into a single file
grep 'query' eggnog_output/ora.pasa_out.prot.split.1.fa.emapper.annotations >  ora.pasa_out.prot.fasta.eggnog_mapper.tsv
grep -h mRNA eggnog_output/*.annotations >> ora.pasa_out.prot.fasta.eggnog_mapper.tsv

# Number of functionally annotated genes
cut -f 1 ora.pasa_out.prot.fasta.eggnog_mapper.tsv | cut -d'.' -f1,2 | tail +2| sed 's/mRNA/Gene/g' | sort | uniq | wc -l



