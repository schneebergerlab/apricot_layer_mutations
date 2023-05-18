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
rm failed_mrna.txt; touch failed_mrna.txt
for i in  00{1..9} 0{10..99} {100..373}; do
    echo $i
    while read m; do
        mrna=$(basename $m .fa)
        if [ -f ${indir}/${mrna}/features.pkl ]; then
            true
        else
            echo $m >> failed_mrna.txt
        fi
    done < fa.list.${i}.txt
done
split -l 320 --numeric-suffixes=1 --additional-suffix=.txt -a 3 failed_mrna.txt fa.list2.

cd /ptmp/mgoel/cur_proteins
for b in 00{1..1}; do
    sbatch -J fa.list2.${b} \
        -o output_%x.txt -e error_%x.txt \
        /raven/u/mgoel/apricot/scripts/SH/raven_submit_scripts/alphafold/jobscript-alphafold2-step_1-msa.sh \
        /raven/u/mgoel/apricot/cur_protein/fa.list2.${b}.txt
done