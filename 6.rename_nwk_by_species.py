#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import Phylo

def read_species_list(file_path):
    with open(file_path, 'r') as f:
        return {line.strip() for line in f if line.strip()}

def mark_tree(tree_file, output_file, target_species):
    tree = Phylo.read(tree_file, "newick")

    for clade in tree.find_clades():
        if clade.name:
            if clade.name in target_species:
                clade.name += "#1"
            else:
                clade.name += "#0"

    Phylo.write(tree, output_file, "newick")

def main():
    parser = argparse.ArgumentParser(description="Label tree leaves with #1 (target species) or #0 (others).")
    parser.add_argument("--input_dir", type=Path, required=True, help="Input directory containing .nwk tree files")
    parser.add_argument("--output_dir", type=Path, required=True, help="Directory to save labeled trees")
    parser.add_argument("--species_file", type=Path, required=True, help="File containing target species list (one per line)")

    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    species_file = args.species_file

    output_dir.mkdir(parents=True, exist_ok=True)
    target_species = read_species_list(species_file)

    for tree_file in input_dir.glob("*_pruned.nwk"):
        output_file = output_dir / tree_file.name
        print(f"Labeling: {tree_file.name}")
        mark_tree(tree_file, output_file, target_species)

    print(f"\n✅ All labelled trees saved to: {output_dir.resolve()}")

if __name__ == "__main__":
    main()

