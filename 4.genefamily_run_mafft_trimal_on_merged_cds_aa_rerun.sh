#!/bin/bash
#SBATCH -c 8
#SBATCH --job-name=cds_direct_pipeline
#SBATCH --time=7-00:00:00
#SBATCH --partition=fos
#SBATCH --qos=fos

# ==== 参数 ====
PPN=8
MAX_GAP=5000
MAX_OVERLAP=5000
RUN_ID="rerun2"
MERGE_SCRIPT=~/shell/smart_merge_multiple_versions.py
#CODON_CHECK_SCRIPT=~/shell/check_codon_integrity.py

# ==== 输入输出路径 ====
CDS_DIR="merged_fragment_cds_rerun1_edit/" #"merged_fragment_codon_edit"
ALIGN_DIR="align_cds_${RUN_ID}"
MERGED_DIR="merged_fragment_cds_${RUN_ID}"
TRIM_DIR="trimal_cds_${RUN_ID}"
#VALID_DIR="codon_valid_${RUN_ID}"
#REJECT_DIR="codon_reject_${RUN_ID}"

# ==== 创建目录 ====
mkdir -p "$ALIGN_DIR" "$MERGED_DIR" "$TRIM_DIR"
#mkdir -p "$VALID_DIR" "$REJECT_DIR"

# ==== 加载模块 ====
module load MAFFT
module load trimAl

# ==== 主流程 ====
for cds_file in ${CDS_DIR}/*.fasta; do
    gene=$(basename "$cds_file" .fasta)
    echo "🔬 Processing $gene"

    TMP_CLEAN="${ALIGN_DIR}/tmp_${gene}_cleaned.fasta"
    ALIGN_FILE="${ALIGN_DIR}/${gene}_mafft_align.fas"
    MERGED_FILE="${MERGED_DIR}/${gene}_merged.fasta"
    TRIMMED_FILE="${TRIM_DIR}/${gene}_trimmed.fasta"
    VALID_FILE="${VALID_DIR}/${gene}_valid.fasta"
    REJECT_FILE="${REJECT_DIR}/${gene}_invalid.fasta"

    # 0️⃣ 安全清洗 header，保留 species_name 和 merge/1 等标记
    sed '/^>/ {
        s/__/_/g;
        s/mergedd/merged/g;
        s/_merged\([0-9]*\)/|merged\1/g;
        s/\([^m]\)_\([0-9]\+\)/\1|\2/g;
        s/_[A-Z]/|/g;
    }' "$cds_file" > "$TMP_CLEAN"

    # 1️⃣ MAFFT 对 CDS 直接比对
    mafft --auto --maxiterate 1000 --thread ${PPN} "$TMP_CLEAN" > "$ALIGN_FILE"
    rm "$TMP_CLEAN"

    # 2️⃣ 贪婪拼接
    python3 "$MERGE_SCRIPT" "$ALIGN_FILE" "$MERGED_FILE" $MAX_GAP $MAX_OVERLAP

    # 3️⃣ TrimAl 修剪
    trimal -in "$MERGED_FILE" -out "$TRIMMED_FILE" -automated1

    # 4️⃣ 密码子完整性检查（如需启用请取消注释）
    #python3 "$CODON_CHECK_SCRIPT" "$TRIMMED_FILE" "$VALID_FILE" "$REJECT_FILE"

    echo "✅ $gene 完成"
done

echo "🎉 所有 CDS 已完成：header清洗 ➜ MAFFT ➜ 贪婪拼接 ➜ TrimAl"

