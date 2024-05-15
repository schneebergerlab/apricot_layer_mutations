# <editor-fold desc="For Currot genome and GFF3">

## Remove the source information from the GFF3
awk '{$2="."; print $0}' FS="\t" OFS="\t" cur.pasa_out.3utr.sort.gff3 > cur.pasa_out.3utr.sort.no_source.gff3
## Remove duplicate entries
agat_sp_fix_features_locations_duplicated.pl --gff cur.pasa_out.3utr.sort.no_source.gff3 -o cur.pasa_out.3utr.sort.no_source.no_dup.gff3
# Convert to .embl format
EMBLmyGFF3 cur.pasa_out.3utr.sort.no_source.no_dup.gff3 cur.genome.v1.fasta \
        --data_class STD \
        --topology linear \
        --molecule_type "genomic DNA" \
        --transl_table 1  \
        --species 36596 \
        --locus_tag ARP \
        --project_id PRJEB71142 \
        --rg schneeberger \
        --author 'Manish Goel' \
        -o cur.genome.v1.embl
# Gzip the file
pigz -p 4 cur.genome.v1.embl
# Validate using webin
java -jar ../../../../pn29fi-dss-0003/software/bin_manish/webin/webin-cli-7.2.0.jar -username Webin-59528 -password Mpipzmpipz2021! -context genome -manifest cur_manifest.txt -validate
# Use the output of validation to remove any remaining duplicated entries using the GFF_rename_for_ENA.py

# Second round of cleaning
## Remove duplicate entries
agat_sp_fix_features_locations_duplicated.pl --gff cur.pasa_out.3utr.sort.no_source.no_dup.clean.gff3 -o cur.pasa_out.3utr.sort2.no_source.no_dup.gff3
# Convert to .embl format
EMBLmyGFF3 cur.pasa_out.3utr.sort2.no_source.no_dup.gff3 cur.genome.v1.fasta \
        --data_class STD \
        --topology linear \
        --molecule_type "genomic DNA" \
        --transl_table 1  \
        --species 36596 \
        --locus_tag ARP \
        --project_id PRJEB71142 \
        --rg schneeberger \
        --author 'Manish Goel' \
        -o cur.genome2.v1.embl
# Gzip the file
pigz -p 4 cur.genome2.v1.embl
# Validate using webin
java -jar ../../../../pn29fi-dss-0003/software/bin_manish/webin/webin-cli-7.2.0.jar -username Webin-59528 -password Mpipzmpipz2021! -context genome -manifest cur_manifest2.txt -validate

# Few positions were still problematic. Fixed them manually in cur.pasa_out.3utr.sort2.no_source.no_dup.gff3
# Affected genes: Gene.7652, Gene.23443, Gene.45811

# After this the submission was succesfully validated and we submit the data
#java -jar ../../../../pn29fi-dss-0003/software/bin_manish/webin/webin-cli-7.2.0.jar -username Webin-59528 -password Mpipzmpipz2021! -context genome -manifest cur_manifest2.txt -submit
# Ran on local PC as webin submit didnt work from the cluster
java -jar /home/ra98jam/software/webin/webin-cli-7.2.0.jar -username Webin-59528 -password Mpipzmpipz2021! -context genome -manifest cur_manifest2.txt -submit

# </editor-fold>



# <editor-fold desc="For Orange Red genome and GFF3">

## Remove the source information from the GFF3
awk '{$2="."; print $0}' FS="\t" OFS="\t" ora.pasa_out.3utr.sort.gff3 > ora.pasa_out.3utr.sort.no_source.gff3
## Remove duplicate entries
agat_sp_fix_features_locations_duplicated.pl --gff ora.pasa_out.3utr.sort.no_source.gff3 -o ora.pasa_out.3utr.sort.no_source.no_dup.gff3
## Flag short introns
agat_sp_flag_short_introns.pl --gff ora.pasa_out.3utr.sort.no_source.no_dup.gff3 --out tmp.gff3
mv tmp.gff3 ora.pasa_out.3utr.sort.no_source.no_dup.gff3

# Shorten the 3'UTR for Gene.43189 by 1bp

# Convert to .embl format
EMBLmyGFF3 ora.pasa_out.3utr.sort.no_source.no_dup.gff3 ora.genome.v1.fasta \
        --data_class STD \
        --topology linear \
        --molecule_type "genomic DNA" \
        --transl_table 1  \
        --species 36596 \
        --locus_tag ARP2 \
        --project_id PRJEB71142 \
        --rg schneeberger \
        --author 'Manish Goel' \
        -o ora.genome.v1.embl
# Gzip the file
pigz -p 4 ora.genome.v1.embl
# Validate using webin
java -jar ../../../../pn29fi-dss-0003/software/bin_manish/webin/webin-cli-7.2.0.jar -username Webin-59528 -password Mpipzmpipz2021! -context genome -manifest ora_manifest.txt -validate
# Use the output of validation to remove any remaining duplicated entries using the GFF_rename_for_ENA.py

# Fix the ERRORs identified in ora.pasa_out.3utr.sort.no_source.no_dup.gff3 using the code in GFF_rename_for_ENA.py

# Second round of cleaning
## Remove duplicate entries
agat_sp_fix_features_locations_duplicated.pl --gff ora.pasa_out.3utr.sort.no_source.no_dup.clean.gff3 -o ora.pasa_out.3utr.sort2.no_source.no_dup.gff3
# Convert to .embl format
EMBLmyGFF3 ora.pasa_out.3utr.sort2.no_source.no_dup.gff3 ora.genome.v1.fasta \
        --data_class STD \
        --topology linear \
        --molecule_type "genomic DNA" \
        --transl_table 1  \
        --species 36596 \
        --locus_tag ARP2 \
        --project_id PRJEB71142 \
        --rg schneeberger \
        --author 'Manish Goel' \
        -o ora.genome2.v1.embl
# Gzip the file
pigz -p 4 ora.genome2.v1.embl
# Validate using webin
java -jar ../../../../pn29fi-dss-0003/software/bin_manish/webin/webin-cli-7.2.0.jar -username Webin-59528 -password Mpipzmpipz2021! -context genome -manifest ora_manifest2.txt -validate
# Use the output of validation to remove any remaining duplicated entries using the GFF_rename_for_ENA.py

# Few positions were still problematic. Fixed them manually in ora.pasa_out.3utr.sort2.no_source.no_dup.gff3
# Affected genes: Gene.19416, Gene.25621, Gene.32832, Gene.43189, Gene.49266
java -jar /home/ra98jam/software/webin/webin-cli-7.2.0.jar -username Webin-59528 -password Mpipzmpipz2021! -context genome -manifest ora_manifest2.txt -submit


# </editor-fold>