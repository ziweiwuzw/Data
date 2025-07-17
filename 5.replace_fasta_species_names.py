#!/usr/bin/env python3
import argparse
from pathlib import Path

def parse_match_log(log_path):
    """解析 simplified 或 fuzzy log，返回 dict[filename][old_name] = new_name"""
    mapping = {}
    current_file = None
    with open(log_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.endswith(".fas - Simplified matches:") or line.endswith(".fas - Fuzzy matches:"):
                current_file = line.split(" - ")[0]
                mapping[current_file] = {}
            elif "→" in line and current_file:
                old, new = map(str.strip, line.split("→"))
                mapping[current_file][old] = new
    return mapping

def merge_mappings(*dicts):
    """合并多个 dict[filename][old] = new"""
    merged = {}
    for d in dicts:
        for fname, submap in d.items():
            if fname not in merged:
                merged[fname] = {}
            merged[fname].update(submap)
    return merged

def replace_fasta_headers(fasta_path, species_map, output_path):
    """替换 fasta 中 header 的物种名，并保留原始名为注释"""
    with open(fasta_path) as f_in, open(output_path, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                header = line[1:].strip()
                species = header.split("|")[0].split("_1")[0]
                rest = header[len(species):]
                if species in species_map:
                    new_species = species_map[species]
                    new_header = f">{new_species}{rest}  # original: {species}\n"
                    f_out.write(new_header)
                else:
                    f_out.write(line)
            else:
                f_out.write(line)

def main():
    parser = argparse.ArgumentParser(description="Replace species names in fasta files using simplified and fuzzy match logs")
    parser.add_argument("--simplified_log", required=True, help="Path to simplified_matches.log")
    parser.add_argument("--fuzzy_log", required=True, help="Path to fuzzy_matches.log")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing original fasta files")
    parser.add_argument("--output_dir", required=True, help="Directory to write updated fasta files")
    args = parser.parse_args()

    simplified_matches = parse_match_log(args.simplified_log)
    fuzzy_matches = parse_match_log(args.fuzzy_log)
    merged_matches = merge_mappings(simplified_matches, fuzzy_matches)

    fasta_dir = Path(args.fasta_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    for fasta_file, species_map in merged_matches.items():
        input_path = fasta_dir / fasta_file
        output_path = output_dir / fasta_file
        if input_path.exists():
            print(f"Processing: {fasta_file}")
            replace_fasta_headers(input_path, species_map, output_path)
        else:
            print(f"[WARNING] File not found: {input_path}")

if __name__ == "__main__":
    main()

