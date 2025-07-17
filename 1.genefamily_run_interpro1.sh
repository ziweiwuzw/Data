#!/bin/bash
#SBATCH -c 32
#SBATCH --job-name=InterPro_only
#SBATCH --time=14-00:00:00
###SBATCH --partition=fos
###SBATCH --qos=fos

# ==== 🔧 可配置部分 ====
IMMUNITY_GENE_DIR="protein_renamed/protein_pro_inflamatory_download/"
PPN=32
RUN_ID=pro_inflamatory_only
INTERPRO_OUT_DIR="interpro_results_immunity${RUN_ID}"
# =======================

# Load module
module load InterProScan

# Make output directory
mkdir -p "$INTERPRO_OUT_DIR"

# Run InterProScan for each immunity gene fasta file
for fasta_file in ${IMMUNITY_GENE_DIR}/*.faa; do
  base_name=$(basename "$fasta_file" .faa)
  interproscan.sh -i "$fasta_file" -b ${INTERPRO_OUT_DIR}/${base_name}_output -goterms -iprlookup -pa -cpu ${PPN}
done

echo "✅ InterProScan completed for immunity_gene files. Results in ${INTERPRO_OUT_DIR}/"

