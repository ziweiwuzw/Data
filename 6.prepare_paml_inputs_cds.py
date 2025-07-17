#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

def fasta_to_phylip(input_fasta, output_phylip, min_spaces=4):
    """从FASTA生成 codon-level sequential PHYLIP 文件"""
    content = Path(input_fasta).read_text()
    entries = re.findall(r'(?m)^>([^\s|#]+)[^\n]*\n([^>]*)', content)

    sequences = []
    for identifier, seq_block in entries:
        clean_seq = ''.join(seq_block.split()).upper()
        genus_species = "_".join(identifier.split("_")[:2])
        if len(clean_seq) % 3 != 0:
            print(f"[WARNING] {genus_species} sequence length not a multiple of 3. Skipping.")
            continue
        sequences.append((genus_species, clean_seq))

    if len(sequences) < 2:
        print(f"[SKIP] Too few sequences in {input_fasta}")
        return

    aln_len = len(sequences[0][1])
    with open(output_phylip, "w") as fout:
        fout.write(f" {len(sequences)} {aln_len}\n")
        for name, seq in sequences:
            fout.write(f"{name}{' ' * min_spaces}{seq}\n")

def generate_ctl(seqfile, treefile, outfile, ctlfile, model):
    ctl_template = f"""
      seqfile = {seqfile}
      treefile = {treefile}
      outfile = {outfile}

      noisy = 9
      verbose = 1
      runmode = 0

      seqtype = 1               * 1 = codons
      CodonFreq = 2             * F3x4
      model = {model}           * 0 = one-ratio, 1 = free-ratio, 2 = branch model
      NSsites = 0
      Mgene = 0
      clock = 0
      icode = 0
      fix_kappa = 0
      kappa = 2
      fix_omega = 0
      omega = 1
      fix_alpha = 1
      alpha = 0.0
      ncatG = 8
      cleandata = 1
      getSE = 0
      RateAncestor = 0
      Small_Diff = .5e-6
    """
    Path(ctlfile).write_text(ctl_template.strip() + "\n")

def fix_phy_spacing(phy_path, min_spaces=4):
    lines = Path(phy_path).read_text().splitlines()
    if not lines or not lines[0].strip():
        return False

    try:
        num_seq, length = map(int, lines[0].strip().split())
    except ValueError:
        print(f"[ERROR] Invalid header in {phy_path.name}")
        return False

    fixed_lines = [f"{num_seq} {length}"]
    modified = False

    for line in lines[1:]:
        match = re.match(r"^(\S+)\s+([A-Z\-]+)$", line.strip())
        if not match:
            name_part = line[:25].strip()
            seq_part = line[25:].strip()
            fixed_lines.append(f"{name_part}{' ' * min_spaces}{seq_part}")
            modified = True
        else:
            name, seq = match.groups()
            fixed_lines.append(f"{name}{' ' * min_spaces}{seq}")

    if modified:
        Path(phy_path).write_text("\n".join(fixed_lines) + "\n")
        print(f"[FIXED] {phy_path.name}")
    else:
        print(f"[OK]    {phy_path.name}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Prepare CDS PHYLIP & codeml control files for PAML branch model")
    parser.add_argument("--fasta_dir", default="updated_fasta_for_tree", help="Directory with *.fas/*.aln")
    parser.add_argument("--tree_dir", default="trimmed_trees", help="Unlabelled tree directory (for model0 and model1)")
    parser.add_argument("--labelled_tree_dir", default="labelled_trees", help="Labelled tree directory (for model2)")
    parser.add_argument("--output_dir", default="paml_inputs", help="Where to save .phy/.ctl files")
    args = parser.parse_args()

    fasta_dir = Path(args.fasta_dir)
    tree_dir = Path(args.tree_dir)
    labelled_tree_dir = Path(args.labelled_tree_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = list(fasta_dir.glob("*.fas")) + list(fasta_dir.glob("*.aln")) + list(fasta_dir.glob("*.fasta"))

    for fasta_file in fasta_files:
        gene = fasta_file.stem
        phylip_file = output_dir / f"{gene}.phy"
        fasta_to_phylip(fasta_file, phylip_file)

        tree0 = tree_dir / f"{gene}_pruned.nwk"
        tree1 = tree_dir / f"{gene}_pruned.nwk"
        tree2 = labelled_tree_dir / f"{gene}_pruned.nwk"

        if not tree0.exists() or not tree2.exists():
            print(f"[WARNING] Missing tree file(s) for {gene}, skipping")
            continue

        ctl0 = output_dir / f"{gene}_model0.ctl"
        out0 = output_dir / f"{gene}_model0.out"
        generate_ctl(phylip_file.name, tree0.resolve(), out0.name, ctl0, model=0)

        ctl1 = output_dir / f"{gene}_model1.ctl"
        out1 = output_dir / f"{gene}_model1.out"
        generate_ctl(phylip_file.name, tree1.resolve(), out1.name, ctl1, model=1)

        ctl2 = output_dir / f"{gene}_model2.ctl"
        out2 = output_dir / f"{gene}_model2.out"
        generate_ctl(phylip_file.name, tree2.resolve(), out2.name, ctl2, model=2)

    print("\n🔍 Checking and fixing .phy spacing...")
    for phy_file in output_dir.glob("*.phy"):
        fix_phy_spacing(phy_file)

    print(f"\n✅ All files written to: {output_dir.resolve()}")

if __name__ == "__main__":
    main()

