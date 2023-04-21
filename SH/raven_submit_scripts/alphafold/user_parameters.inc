## Set user parameters for alphafold, common to MSA and prediction jobs

# Comma-separated list of protein data, passed to alphafold via --fasta_paths
FASTA_PATHS=mRNA.9970.1.fa

# Output directory to write the results to, passed via --output
OUTPUT_DIR=/ptmp/${USER}/mRNA.9970.1

# Preset for alphafold, passed via --db_preset
PRESET="full_dbs"



# --- More flags and paths, picked from <https://github.com/deepmind/alphafold/blob/main/docker/run_docker.py> ---

# Path to the Uniref90 database for use by JackHMMER.
#uniref90_database_path = os.path.join(FLAGS.data_dir, 'uniref90', 'uniref90.fasta')
uniref90_database_path=${ALPHAFOLD_DATA}/uniref90/uniref90.fasta

# Path to the MGnify database for use by JackHMMER.
#mgnify_database_path = os.path.join( FLAGS.data_dir, 'mgnify', 'mgy_clusters_2018_12.fa')
mgnify_database_path=${ALPHAFOLD_DATA}/mgnify/mgy_clusters_2022_05.fa

# Path to the Uniprot database for use by JackHMMER.
#uniprot_database_path = os.path.join(FLAGS.data_dir, 'uniprot', 'uniprot.fasta')
uniprot_database_path=${ALPHAFOLD_DATA}/uniprot/uniprot.fasta

# Path to the BFD database for use by HHblits.
#bfd_database_path = os.path.join(FLAGS.data_dir, 'bfd', 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt')
bfd_database_path=${ALPHAFOLD_DATA}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt

# Path to the Small BFD database for use by JackHMMER.
# small_bfd_database_path = os.path.join(FLAGS.data_dir, 'small_bfd', 'bfd-first_non_consensus_sequences.fasta')
small_bfd_database_path=${ALPHAFOLD_DATA}/small_bfd/bfd-first_non_consensus_sequences.fasta

# Path to the PDB seqres database for use by hmmsearch.
#pdb_seqres_database_path = os.path.join(FLAGS.data_dir, 'pdb_seqres', 'pdb_seqres.txt')
pdb_seqres_database_path=${ALPHAFOLD_DATA}/pdb_seqres/pdb_seqres.txt

# Path to the Uniref30 database for use by HHblits.
# uniref30_database_path = os.path.join( FLAGS.data_dir, 'uniref30', 'UniRef30_2021_03')
uniref30_database_path=${ALPHAFOLD_DATA}/uniref30/UniRef30_2021_03

# Path to the PDB70 database for use by HHsearch.
# pdb70_database_path = os.path.join(FLAGS.data_dir, 'pdb70', 'pdb70')
pdb70_database_path=${ALPHAFOLD_DATA}/pdb70/pdb70

# Path to a directory with template mmCIF structures, each named <pdb_id>.cif.
#template_mmcif_dir = os.path.join(FLAGS.data_dir, 'pdb_mmcif', 'mmcif_files')
template_mmcif_dir=${ALPHAFOLD_DATA}/pdb_mmcif/mmcif_files

# Path to a file mapping obsolete PDB IDs to their replacements.
#obsolete_pdbs_path = os.path.join(FLAGS.data_dir, 'pdb_mmcif', 'obsolete.dat')
obsolete_pdbs_path=${ALPHAFOLD_DATA}/pdb_mmcif/obsolete.dat
