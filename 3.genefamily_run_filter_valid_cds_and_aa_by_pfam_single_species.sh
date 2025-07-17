#!/bin/bash
#SBATCH -c 1
#SBATCH --job-name=filter_pfam_cds_aa
#SBATCH --time=2-00:00:00
#SBATCH --partition=fos
#SBATCH --qos=fos

# ==== 配置参数 ====
MIN_HITS=1
PFAM_SUMMARY="/data/scratch/projects/punim2285/ziwei_data/ISG_primates_ncbi/pfam_summary_ISG_inflamatory_edit.tsv" #"/data/scratch/projects/punim2285/ziwei_data/anti_pro-inflamatory/genefamily_data/pfam_summary_pro_inflamatory_edit.tsv"
PFAM_TXT_DIR="interpro_results_all_rest_cds-0.85" #"interpro_results_rest_cds-0.85_pro"
CDS_DIR="extracted_cds_all_rest_cds-0.85" #"extracted_cds_rest_cds-0.85_pro"
OUTPUT_DIR="filtered_cds_by_pfam_min${MIN_HITS}_rest_ISG" #"filtered_cds_by_pfam_min${MIN_HITS}_rest_pro_inflamatory"

module load Biopython

echo "🚀 Running filter_valid_cds_and_aa_by_pfam..."
python3 ~/shell/genefamily_filter_valid_cds_and_aa_by_pfam_single_species.py \
  --pfam_summary "$PFAM_SUMMARY" \
  --pfam_txt_dir "$PFAM_TXT_DIR" \
  --cds_dir "$CDS_DIR" \
  --output_dir "$OUTPUT_DIR" \
  --min_hits "$MIN_HITS"

echo "✅ Done. Filtered CDS and AA written to: $OUTPUT_DIR"

