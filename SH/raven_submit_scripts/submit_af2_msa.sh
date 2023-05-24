# Create files containing names fasta files where each fasta file is a mRNA sequence
# Number of mRNAs: 37285
# Number of mRNAs per file: 100
# Total number of files: 373
tar -tf protein_fasta.tar.gz > mrna.fa.list.txt
split -l 100 --numeric-suffixes=1 --additional-suffix=.txt -a 3 mrna.fa.list.txt fa.list.


# Run MSA jobs for individual fasta lists
## This can only submit three sbatch jobs at a time (job submit limit 300)
cd /ptmp/mgoel/cur_proteins
for b in 00{4..9} 0{10..99} {100..373}; do
    sbatch -J falist.${b} \
        -o output_%x.txt -e error_%x.txt \
        /raven/u/mgoel/apricot/scripts/SH/raven_submit_scripts/alphafold/jobscript-alphafold2-step_1-msa.sh \
        /raven/u/mgoel/apricot/cur_protein/fa.list.${b}.txt
done


## Code to check if any of the msa jobs are unfinished/crashed.
cd /u/mgoel/apricot/cur_protein
indir=/ptmp/mgoel/cur_proteins/af2_msa/
touch failed_mrna2.txt
while read m; do
    mrna=$(basename $m .fa)
    if [ -f ${indir}/${mrna}/features.pkl ]; then
        true
    else
        echo $m >> failed_mrna2.txt
    fi
done < failed_mrna.txt
mv failed_mrna2.txt failed_mrna.txt
split -l 12 --numeric-suffixes=1 --additional-suffix=.txt -a 3 failed_mrna.txt fa.list2.

cd /ptmp/mgoel/cur_proteins
for b in 00{1..2} ; do
    sbatch -J fa.list2.${b} \
        -o output_%x.txt -e error_%x.txt \
        /raven/u/mgoel/apricot/scripts/SH/raven_submit_scripts/alphafold/jobscript-alphafold2-step_1-msa.sh \
        /raven/u/mgoel/apricot/cur_protein/fa.list2.${b}.txt
done


# Submit prediction job with GPU nodes
## Create test job mRNA list
head -12 mrna.fa.list.txt > test.fa.list
## Run test job with 3 nodes and 4 tasks per node
cd /ptmp/mgoel/cur_proteins
sbatch -J test.predict \
    -o output_%x.txt -e error_%x.txt \
    /raven/u/mgoel/apricot/scripts/SH/raven_submit_scripts/alphafold/jobscript-alphafold2-step_2-prediction.sh \
    /raven/u/mgoel/apricot/cur_protein/test.fa.list

## mRNA for which MSA worked
cd /u/mgoel/apricot/cur_protein
indir=/ptmp/mgoel/cur_proteins/af2_msa/
touch mrna_msa.txt
while read m; do
    mrna=$(basename $m .fa)
    if [ -f ${indir}/${mrna}/features.pkl ]; then
        echo $m >> mrna_msa.txt
    else
        true
    fi
done < mrna.fa.list.txt

## Split input mRNA files
rm fa.list.*
split -l 960 --numeric-suffixes=1 --additional-suffix=.txt -a 3 mrna_msa.txt fa.list.

cd /ptmp/mgoel/cur_proteins
for b in 00{1..9} 0{10..39} ; do
    sbatch -J predict.${b} \
        -o out_%x.txt -e err_%x.txt \
        /raven/u/mgoel/apricot/scripts/SH/raven_submit_scripts/alphafold/jobscript-alphafold2-step_2-prediction.sh \
        /raven/u/mgoel/apricot/cur_protein/fa.list.${b}.txt
done

