#!/bin/bash
#SBATCH -c 4
#SBATCH --job-name=cds_codon_pipeline
#SBATCH --time=7-00:00:00
#SBATCH --partition=fos
#SBATCH --qos=fos

# ==== 参数 ====
PPN=4
MAX_GAP=500
MAX_OVERLAP=500
RUN_ID="pro_run"
PAL2NAL=~/bin/pal2nal.pl  # 请确认 pal2nal.pl 路径正确

# ==== 输入输出路径 ====
INPUT_DIR="filtered_cds_by_pfam_min1_all_pro_inflamatory" #"filtered_cds_by_pfam_min1_all"
ALIGN_AA_DIR="align_aa_${RUN_ID}"
CODON_ALIGN_DIR="align_codon_${RUN_ID}"
MERGE_DIR="merged_fragment_codon_${RUN_ID}"
TRIM_DIR="trimal_codon_${RUN_ID}"
#VALID_DIR="codon_valid_${RUN_ID}"
#REJECT_DIR="codon_reject_${RUN_ID}"

# ==== 脚本路径 ====
MERGE_SCRIPT=~/shell/smart_merge_multiple_versions.py
#CODON_CHECK_SCRIPT=~/shell/check_codon_integrity.py

# ==== 创建目录 ====
mkdir -p "$ALIGN_AA_DIR" "$CODON_ALIGN_DIR" "$MERGE_DIR" "$TRIM_DIR"
#mkdir -p "$VALID_DIR" "$REJECT_DIR"

# ==== 加载模块 ====
module load MAFFT
module load trimAl

# ==== 主流程 ====
for faa_file in ${INPUT_DIR}/*.faa; do
    gene=$(basename "$faa_file" _filtered_min1.faa)
    cds_file="${INPUT_DIR}/${gene}_filtered_min1.cds.fasta"

    echo "🔬 Processing $gene"

    ALIGN_AA="${ALIGN_AA_DIR}/${gene}_aligned.faa"
    CODON_ALIGN="${CODON_ALIGN_DIR}/${gene}_codon_aligned.fasta"
    MERGED="${MERGE_DIR}/${gene}_merged_codon.fasta"
    TRIMMED="${TRIM_DIR}/${gene}_trimmed.fasta"
    VALID="${VALID_DIR}/${gene}_valid.fasta"
    REJECT="${REJECT_DIR}/${gene}_invalid.fasta"

    # 0️⃣ 安全清洗 AA header
    TMP_AA="${ALIGN_AA_DIR}/tmp_${gene}_cleaned.faa"
    sed '/^>/ {
        s/__/_/g;
        s/mergedd/merged/g;
        s/_merged\([0-9]\+\)/|merged\1/g;
        s/\([^m]\)_\([0-9]\+\)/\1|\2/g;
        s/_[A-Z]/|/g;
    }' "$faa_file" > "$TMP_AA"

    # 1️⃣ MAFFT 蛋白比对
    mafft --auto --maxiterate 1000 --thread ${PPN} "$TMP_AA" > "$ALIGN_AA"
    rm "$TMP_AA"

    # 2️⃣ pal2nal 密码子对齐（必须步骤）
    "$PAL2NAL" "$ALIGN_AA" "$cds_file" -output fasta > "$CODON_ALIGN"

    # 3️⃣ 贪婪拼接
    python3 "$MERGE_SCRIPT" "$CODON_ALIGN" "$MERGED" $MAX_GAP $MAX_OVERLAP

    # 4️⃣ TrimAl 修剪
    trimal -in "$MERGED" -out "$TRIMMED" -automated1

    # 5️⃣ 密码子完整性检查（暂不启用）
    #python3 "$CODON_CHECK_SCRIPT" "$TRIMMED" "$VALID" "$REJECT"

    echo "✅ $gene 完成"
done

echo "🎉 所有基因已完成：蛋白对齐 ➜ pal2nal ➜ 拼接 ➜ TrimAl"
#echo "🎉 有效序列保存在 $VALID_DIR，无效序列保存在 $REJECT_DIR"

