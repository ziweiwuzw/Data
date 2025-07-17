#!/usr/bin/env python3
import os
import glob
import pandas as pd
from collections import defaultdict

# ===== 🔧 配置区域 =====
input_dir = "interpro_results_immunitypro_inflamatory_only"   # 👈 你只需要改这里
output_file = "pfam_summary_pro_inflamatory.tsv"
# ======================

input_files = glob.glob(os.path.join(input_dir, "*.txt"))

pfam_count = defaultdict(int)
pfam_proteins = defaultdict(set)
pfam_names = {}

for file in input_files:
    try:
        df = pd.read_csv(file, sep="\t", dtype=str)
    except Exception as e:
        print(f"⚠️ Skipping {file}: {e}")
        continue

    if "Database" not in df.columns or "Domain_ID" not in df.columns:
        print(f"⚠️ Skipping {file}: missing required columns")
        continue

    pfam_df = df[df["Database"] == "Pfam"]

    for _, row in pfam_df.iterrows():
        pfam_id = row["Domain_ID"]
        name = row["Domain_Name"]
        protein = row["Protein_ID"]

        pfam_count[pfam_id] += 1
        pfam_proteins[pfam_id].add(protein)
        pfam_names[pfam_id] = name

# 输出汇总结果
with open(output_file, "w") as out:
    out.write("Pfam_ID\tDomain_Name\tTotal_Count\tUnique_Proteins\n")
    for pfam_id in sorted(pfam_count, key=pfam_count.get, reverse=True):
        name = pfam_names.get(pfam_id, "-")
        count = pfam_count[pfam_id]
        unique_proteins = len(pfam_proteins[pfam_id])
        out.write(f"{pfam_id}\t{name}\t{count}\t{unique_proteins}\n")

print(f"✅ Pfam summary saved to {output_file}")

