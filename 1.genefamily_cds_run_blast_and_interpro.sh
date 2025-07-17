#!/bin/bash
#SBATCH -c 32
#SBATCH --job-name=CDS_BLAST_InterPro
#SBATCH --time=14-00:00:00
###SBATCH --partition=gpu-l40s-preempt  ##gpu-a100-mig #fos
####SBATCH --qos=gpu-a100-mig   #fos
###SBATCH --gres=gpu:1

# ==== 🔧 配置 ====
PPN=32
RUN_ID="rest3_cds-0.85"
PRIMATES_CDS_DIR="/data/scratch/projects/punim2285/ziwei_data/150_genome_cds_with_protein"
IMMUNITY_GENE_DIR="protein_renamed_rest3"

BLAST_OUT_DIR="blastn_results_${RUN_ID}"
BLAST_DB_DIR="blast_dbs_${RUN_ID}"
EXTRACTED_FASTA_DIR="extracted_cds_${RUN_ID}"
INTERPRO_OUT_DIR="interpro_results_${RUN_ID}"
PY_SCRIPT="$HOME/shell/extract_cds_and_translate_dual.py"

module load InterProScan
module load BEDTools
module load Biopython
module load BLAST

mkdir -p "$BLAST_OUT_DIR"
mkdir -p "$BLAST_DB_DIR"
mkdir -p "$EXTRACTED_FASTA_DIR"
mkdir -p "$INTERPRO_OUT_DIR"

echo "🚀 Step 1: Creating BLAST databases and running TBLASTN..."
for cds_file in ${PRIMATES_CDS_DIR}/*.fasta; do
  species_name=$(basename "$cds_file" .fasta)
  db_path="${BLAST_DB_DIR}/${species_name}"

  # ✅ 跳过已存在的 BLAST 数据库
  if [[ ! -f "${db_path}.nin" ]]; then
    echo "🔨 Building BLAST DB: $species_name"
    makeblastdb -in "$cds_file" -dbtype nucl -out "$db_path"
  else
    echo "⏩ Skipping DB build: $species_name (already exists)"
  fi

  for immune_protein in ${IMMUNITY_GENE_DIR}/*.faa; do
    gene_name=$(basename "$immune_protein" .faa)
    sample_name="${gene_name}_vs_${species_name}"

    tblastn -query "$immune_protein" -db "$db_path" -outfmt 6 \
      -evalue 1e-5 -max_target_seqs 500 -num_threads ${PPN} \
      -out "${BLAST_OUT_DIR}/${sample_name}.txt"
  done
done

echo "✅ Step 2: Filtering BLAST results with identity > 85..."
for blast_file in ${BLAST_OUT_DIR}/*.txt; do
  awk '{
    if ($3 > 85 && $11 < 2e-28) print $0;
  }' "$blast_file" > "${blast_file%.txt}_filtered.txt"
done

echo "🔄 Step 3–5: Extracting CDS + Translate + InterProScan..."
for filtered_file in ${BLAST_OUT_DIR}/*_filtered.txt; do
  sample_name=$(basename "$filtered_file" _filtered.txt)
  species_name=$(echo "$sample_name" | sed 's/.*_vs_//')
  cds_path="${PRIMATES_CDS_DIR}/${species_name}.fasta"

  if [[ ! -f "$cds_path" ]]; then
    echo "❌ Missing CDS file: $cds_path"
    continue
  fi

  bed_file="tmp_${sample_name}.bed"
  awk '{
    chr=$2;
    start=($9<$10 ? $9 : $10) - 1;
    end=($9>$10 ? $9 : $10);
    print chr "\t" start "\t" end;
  }' "$filtered_file" | sort -k1,1 -k2,2n > "$bed_file"

  merged_bed="tmp_${sample_name}_merged.bed"
  bedtools merge -i "$bed_file" > "$merged_bed"

  output_prefix="${EXTRACTED_FASTA_DIR}/${sample_name}_merged_hits"
  python3 "$PY_SCRIPT" "$merged_bed" "$cds_path" "$output_prefix"

  protein_file="${output_prefix}.faa"
  if [[ -s "$protein_file" ]]; then
    interproscan.sh -i "$protein_file" -b "${INTERPRO_OUT_DIR}/${sample_name}_output" \
      -goterms -iprlookup -pa -cpu ${PPN}
    echo "✅ InterProScan done: $sample_name"
  else
    echo "⚠️ No protein sequences extracted for $sample_name"
  fi

  rm -f "$bed_file" "$merged_bed"
done

echo "🎉 All tasks completed:"
echo "  ➔ TBLASTN: $BLAST_OUT_DIR"
echo "  ➔ Extracted CDS: $EXTRACTED_FASTA_DIR"
echo "  ➔ InterProScan results: $INTERPRO_OUT_DIR"

