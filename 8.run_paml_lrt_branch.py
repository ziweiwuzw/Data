#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from scipy.stats import chi2

def extract_lnL(paml_out_path):
    """提取 codeml 输出中的 log likelihood 值"""
    text = Path(paml_out_path).read_text()
    for line in text.splitlines():
        if line.startswith("lnL(ntime:"):
            match = re.search(r"lnL.*?:\s*(-?\d+\.\d+)", line)
            if match:
                return float(match.group(1))
    return None

def extract_treefile_from_ctl(ctl_path):
    """从 .ctl 文件中提取 treefile 路径"""
    text = Path(ctl_path).read_text()
    for line in text.splitlines():
        if line.strip().startswith("treefile"):
            parts = line.strip().split("=")
            if len(parts) == 2:
                path = parts[1].strip()
                tree_path = (ctl_path.parent / path).resolve()
                if tree_path.exists():
                    return tree_path
    return None

def count_omega_groups(tree_path):
    """从树文件中读取有多少个不同的 ω 分组（标签）"""
    if not tree_path or not tree_path.exists():
        return None
    text = Path(tree_path).read_text()
    labels = re.findall(r"#(\d+)", text)
    return len(set(labels)) + 1  # +1 for background ω0

def main():
    parser = argparse.ArgumentParser(description="Compare model0 and model2 codeml results for branch model LRT (auto df from ctl tree)")
    parser.add_argument("--input_dir", default="paml_inputs", help="Directory with codeml output and ctl files")
    parser.add_argument("--output_tsv", type=str, help="Optional output TSV file")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    gene_names = sorted({f.stem.rsplit("_model", 1)[0] for f in input_dir.glob("*_model0.out")})

    results = []
    print("Gene\tlnL_model0\tlnL_model2\tLRT\tp-value\tdf\tSignificant")

    for gene in gene_names:
        f0 = input_dir / f"{gene}_model0.out"
        f2 = input_dir / f"{gene}_model2.out"
        ctl2 = input_dir / f"{gene}_model2.ctl"

        if not f0.exists() or not f2.exists() or not ctl2.exists():
            print(f"{gene}\t[Missing file(s)]")
            continue

        lnL0 = extract_lnL(f0)
        lnL2 = extract_lnL(f2)
        if lnL0 is None or lnL2 is None:
            print(f"{gene}\t[Error extracting lnL]")
            continue

        tree = extract_treefile_from_ctl(ctl2)
        df = count_omega_groups(tree)
        if df is None:
            print(f"{gene}\t[Missing or unreadable tree]")
            continue
        df -= 1  # degree of freedom = num_groups - 1

        LRT = 2 * (lnL2 - lnL0)
        pval = chi2.sf(LRT, df)
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""

        print(f"{gene}\t{lnL0:.2f}\t{lnL2:.2f}\t{LRT:.2f}\t{pval:.4g}\t{df}\t{sig}")
        results.append((gene, lnL0, lnL2, LRT, pval, df, sig))

    if args.output_tsv:
        with open(args.output_tsv, "w") as out:
            out.write("Gene\tlnL_model0\tlnL_model2\tLRT\tp-value\tdf\tSignificant\n")
            for r in results:
                out.write(f"{r[0]}\t{r[1]:.2f}\t{r[2]:.2f}\t{r[3]:.2f}\t{r[4]:.4g}\t{r[5]}\t{r[6]}\n")
        print(f"\n✅ Results saved to: {args.output_tsv}")

if __name__ == "__main__":
    main()

