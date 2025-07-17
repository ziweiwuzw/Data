#!/usr/bin/env python3
import sys
import re
import csv
from Bio import SeqIO
from collections import defaultdict

def get_alignment_range(seq_str):
    start = next((i for i, c in enumerate(seq_str) if c != '-'), None)
    end = len(seq_str) - next((i for i, c in enumerate(reversed(seq_str)) if c != '-'), None)
    return start, end

def greedy_merge(seq_info, max_gap=55, max_overlap=20, ids=None):
    used = [False] * len(seq_info)
    merged_results = []
    merge_flags = []
    merged_id_map = []

    for i in range(len(seq_info)):
        if used[i]:
            continue

        curr_seq, start_i, end_i = seq_info[i]
        merged = list(curr_seq)
        used[i] = True
        merged_flag = False
        absorb_ids = []

        for j in range(i + 1, len(seq_info)):
            if used[j]:
                continue
            next_seq, start_j, end_j = seq_info[j]

            if start_j is None or end_j is None:
                continue

            gap = start_j - end_i
            overlap = end_i - start_j

            if (0 <= gap <= max_gap) or (0 <= overlap < max_overlap):
                for k in range(start_j, end_j):
                    if next_seq[k] != '-':
                        merged[k] = next_seq[k]
                end_i = max(end_i, end_j)
                used[j] = True
                merged_flag = True
                absorb_ids.append(ids[j])

        merged_str = ''.join(merged)
        merged_results.append(merged_str)
        merge_flags.append(merged_flag)
        merged_id_map.append((ids[i], absorb_ids))

    return merged_results, merge_flags, merged_id_map

def decide_final_id(original_id, merged_flag, species_name, nonmerge_counter):
    """
    输出 ID 命名逻辑：
    - 若原始 ID 中包含 merge：如果合并了则追加 _merge，否则保留原样
    - 若未包含 merge：合并则命名为 Species_merge，否则编号
    """
    already_merged = bool(re.search(r'merge', original_id, re.IGNORECASE))

    if already_merged:
        if merged_flag:
            return f"{original_id}_merge"
        else:
            return original_id.split()[0]
    else:
        if merged_flag:
            return f"{species_name}_merge"
        else:
            return f"{species_name}_{nonmerge_counter}"

def main(input_fasta, output_fasta, max_gap=55, max_overlap=20, output_csv="id_mapping.csv", missing_output="missing_ids.txt"):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    species_seq_dict = defaultdict(list)

    for seq in sequences:
        species = seq.id.split('|')[0]
        full_id = seq.id
        species_seq_dict[species].append((full_id, str(seq.seq)))

    id_mapping = []
    all_original_ids = set()

    with open(output_fasta, "w") as f:
        for species, id_seq_list in species_seq_dict.items():
            if not id_seq_list:
                continue

            seqs = [seq for _, seq in id_seq_list]
            ids = [full_id for full_id, _ in id_seq_list]
            all_original_ids.update(ids)

            seq_info = [(s, *get_alignment_range(s)) for s in seqs]
            seq_info.sort(key=lambda x: (x[1] if x[1] is not None else float('inf')))
            merged_versions, merge_flags, merge_id_map = greedy_merge(seq_info, max_gap, max_overlap, ids)

            nonmerged_counter = 1
            for i, (merged_seq, merged_flag) in enumerate(zip(merged_versions, merge_flags)):
                species_name = ids[i].split('|')[0]
                new_id = decide_final_id(ids[i], merged_flag, species_name, nonmerged_counter)

                if not merged_flag and "merge" not in new_id:
                    nonmerged_counter += 1

                f.write(f">{new_id}\n{merged_seq}\n")
                id_mapping.append([ids[i], new_id])

    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Original ID", "New ID"])
        writer.writerows(id_mapping)

    written_ids = set()
    with open(output_fasta) as f:
        for line in f:
            if line.startswith(">"):
                written_ids.add(line[1:].strip().split()[0])

    missing_ids = all_original_ids - written_ids
    if missing_ids:
        print(f"\n⚠️ Missing {len(missing_ids)} sequence(s) from output:")
        for mid in sorted(missing_ids):
            print(f"  - {mid}")
        with open(missing_output, "w") as fout:
            fout.writelines([mid + "\n" for mid in sorted(missing_ids)])
    else:
        print("\n✅ All original sequences were retained or merged properly.")

if __name__ == "__main__":
    if len(sys.argv) not in [3, 5]:
        print("Usage: smart_merge_multiple_versions.py <input_fasta> <output_fasta> [max_gap max_overlap]")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    max_gap = int(sys.argv[3]) if len(sys.argv) == 5 else 55
    max_overlap = int(sys.argv[4]) if len(sys.argv) == 5 else 20
    output_csv = "id_mapping.csv"
    missing_output = "missing_ids.txt"

    main(input_fasta, output_fasta, max_gap=max_gap, max_overlap=max_overlap,
         output_csv=output_csv, missing_output=missing_output)

