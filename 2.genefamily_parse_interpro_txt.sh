#!/bin/bash

# ==== 🔧 配置 ====
INPUT_DIR="interpro_results_immunityonly" #interpro_results_all_rest_cds-0.85 #"interpro_results_first_cds-0.85" #"interpro_results_pro_inflamatory_first_cds-0.85" #"interpro_results_immunitypro_inflamatory_only" #"interpro_results_rest_cds-0.85_anti" #"interpro_results_anti_inflamatory_first_cds-0.85/"
OUTPUT_PREFIX="interpro"
# =================

for tsv_file in ${INPUT_DIR}/*_output.tsv; do
  basename=$(basename "$tsv_file" _output.tsv)
  output_file=${INPUT_DIR}/"${OUTPUT_PREFIX}_${basename}.txt"

  echo -e "Protein_ID\tDatabase\tDomain_ID\tDomain_Name\tStart\tEnd\tGO_Terms" > "$output_file"

  while IFS=$'\t' read -r protein md5 length db sig_id sig_desc start end score status date ipr_id ipr_desc go_terms pathways; do
    [[ "$protein" == "#"* || -z "$protein" ]] && continue

    echo -e "${protein}\t${db}\t${sig_id}\t${sig_desc}\t${start}\t${end}\t${go_terms}" >> "$output_file"
  done < "$tsv_file"

  echo "✅ Parsed $tsv_file → $output_file"
done

