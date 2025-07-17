#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import Phylo
import copy
import difflib

def extract_species(fasta_path):
    """提取 fasta 文件中的物种名（保留原始名）"""
    species = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                species_name = header.split("|")[0].split("_1")[0]
                species.add(species_name)
    return species

def prune_and_save(tree, species_set, output_path, unmatched_log, simplified_log, fuzzy_log, fasta_filename):
    """修剪树并记录未匹配、简化匹配和模糊匹配物种"""
    tree_copy = copy.deepcopy(tree)
    tree_species = set(term.name for term in tree_copy.get_terminals())

    matched_species = set()
    unmatched_species = set()
    simplified_matches = {}
    fuzzy_matches = {}

    for sp in species_set:
        # 1. 精确匹配
        if sp in tree_species:
            matched_species.add(sp)
            continue

        # 2. 简化前两段尝试
        simplified = "_".join(sp.split("_")[:2])
        if simplified in tree_species:
            matched_species.add(simplified)
            simplified_matches[sp] = simplified
            continue

        # 3. 模糊匹配
        close = difflib.get_close_matches(sp, tree_species, n=1, cutoff=0.8)
        if close:
            matched_species.add(close[0])
            fuzzy_matches[sp] = close[0]
        else:
            unmatched_species.add(sp)

    # unmatched
    with unmatched_log.open("a") as log:
        if not matched_species:
            print(f"[ERROR] All species in {fasta_filename} are missing. Skipping.")
            log.write(f"[ALL MISSING] {fasta_filename} - No species matched:\n")
            for sp in sorted(unmatched_species):
                log.write(f"  - {sp}\n")
            log.write("\n")
            return
        elif unmatched_species:
            print(f"[WARNING] Unmatched species in {fasta_filename}:")
            log.write(f"{fasta_filename} - Unmatched species:\n")
            for sp in sorted(unmatched_species):
                print(f"  - {sp}")
                log.write(f"  - {sp}\n")
            log.write("\n")

    # simplified matches
    with simplified_log.open("a") as log:
        if simplified_matches:
            log.write(f"{fasta_filename} - Simplified matches:\n")
            for orig, simp in sorted(simplified_matches.items()):
                log.write(f"  {orig}  →  {simp}\n")
            log.write("\n")

    # fuzzy matches
    with fuzzy_log.open("a") as log:
        if fuzzy_matches:
            log.write(f"{fasta_filename} - Fuzzy matches:\n")
            for orig, matched in sorted(fuzzy_matches.items()):
                log.write(f"  {orig}  →  {matched}\n")
            log.write("\n")

    # 修剪树
    for term in tree_copy.get_terminals():
        if term.name not in matched_species:
            tree_copy.prune(term)

    Phylo.write(tree_copy, output_path, "newick")

def main():
    parser = argparse.ArgumentParser(description="Prune a species tree based on taxa found in .fas files.")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing .fas files")
    parser.add_argument("--tree_file", required=True, help="Newick tree file (full species tree)")
    parser.add_argument("--output_dir", required=True, help="Directory to save pruned trees and logs")
    args = parser.parse_args()

    fasta_dir = Path(args.fasta_dir)
    tree_file = Path(args.tree_file)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    unmatched_log = output_dir / "unmatched_species.log"
    simplified_log = output_dir / "simplified_matches.log"
    fuzzy_log = output_dir / "fuzzy_matches.log"

    unmatched_log.write_text("")
    simplified_log.write_text("")
    fuzzy_log.write_text("")

    fasta_files = list(fasta_dir.glob("*.fas"))
    if not fasta_files:
        print("[ERROR] No .fas files found in the specified directory.")
        return

    print(f"Loading species tree: {tree_file}")
    full_tree = Phylo.read(tree_file, "newick")
    print(f"Tree contains {len(full_tree.get_terminals())} taxa")

    for fasta_path in fasta_files:
        print(f"\nProcessing: {fasta_path.name}")
        species = extract_species(fasta_path)
        output_tree_path = output_dir / (f"{fasta_path.stem}_pruned.nwk")
        prune_and_save(full_tree, species, output_tree_path, unmatched_log, simplified_log, fuzzy_log, fasta_path.name)
        print(f"Finished: {fasta_path.name}")

if __name__ == "__main__":
    main()

