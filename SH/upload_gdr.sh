cwd=/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/apricot_leaf/assembly_upload_gdr/
cd $cwd

# Using the genomes and final GFF files (uploaded to ENA). Repeatmasker output from /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker

# Get fasta files for annotation elements
# CURROT
## Get gene fasta (https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html)
agat_sp_extract_sequences.pl -g cur.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f cur.genome.v1.fasta -t gene -o cur.pasa_out.3utr.sort2.no_source.no_dup.gff3.gene.fasta

## Get mRNA fasta
agat_sp_extract_sequences.pl -g cur.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f cur.genome.v1.fasta  -t exon --merge -o cur.pasa_out.3utr.sort2.no_source.no_dup.gff3.mrna.fasta

# Get CDS fasta
agat_sp_extract_sequences.pl -g cur.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f cur.genome.v1.fasta -t cds -o cur.pasa_out.3utr.sort2.no_source.no_dup.gff3.CDS.fasta

# Get prot fasta
agat_sp_extract_sequences.pl -g cur.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f cur.genome.v1.fasta -t cds -p -o cur.pasa_out.3utr.sort2.no_source.no_dup.gff3.prot.fasta


# OrangeRed
## Get gene fasta (https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html)
agat_sp_extract_sequences.pl -g ora.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f ora.genome.v1.fasta -t gene -o ora.pasa_out.3utr.sort2.no_source.no_dup.gff3.gene.fasta

## Get mRNA fasta
agat_sp_extract_sequences.pl -g ora.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f ora.genome.v1.fasta  -t exon --merge -o ora.pasa_out.3utr.sort2.no_source.no_dup.gff3.mrna.fasta

# Get CDS fasta
agat_sp_extract_sequences.pl -g ora.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f ora.genome.v1.fasta -t cds -o ora.pasa_out.3utr.sort2.no_source.no_dup.gff3.CDS.fasta

# Get prot fasta
agat_sp_extract_sequences.pl -g ora.pasa_out.3utr.sort2.no_source.no_dup.gff3 -f ora.genome.v1.fasta -t cds -p -o ora.pasa_out.3utr.sort2.no_source.no_dup.gff3.prot.fasta


