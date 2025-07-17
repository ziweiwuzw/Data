#!/usr/bin/env python3
import argparse
import math
from pathlib import Path
import subprocess

BATCH_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=codeml_batch{batch_id}
#SBATCH --partition=fos
#SBATCH --qos=fos
#SBATCH --time=14-00:00:00
#SBATCH --mem=4G
#SBATCH -c 1

module load Biopython
module load PAML

PAML_INPUTS="{ctl_dir}"
TMP_BASE="$PAML_INPUTS/tmp_batch_{batch_id}"
mkdir -p "$TMP_BASE"

{codeml_blocks}

wait
"""

CODEML_BLOCK = """
# === {job_name} ===
mkdir -p "$TMP_BASE/{job_name}"
cd "$TMP_BASE/{job_name}"

cp "$PAML_INPUTS/{ctl_file}" ./run.ctl
cp "$PAML_INPUTS/{phy_file}" ./

codeml run.ctl > {job_name}.log 2>&1

mv {job_name}.out "$PAML_INPUTS/"
mv {job_name}.log "$PAML_INPUTS/"

cd "$PAML_INPUTS"
rm -rf "$TMP_BASE/{job_name}"
"""

def main():
    parser = argparse.ArgumentParser(description="Submit codeml batches with isolated temp folders")
    parser.add_argument("--ctl_dir", required=True, help="Directory with .ctl and .phy files")
    parser.add_argument("--cpu", type=int, default=4, help="Number of concurrent codeml jobs per batch")
    parser.add_argument("--slurm_dir", default="slurm_scripts", help="Where to store generated SLURM scripts")
    args = parser.parse_args()

    ctl_dir = Path(args.ctl_dir).resolve()
    ctl_files = sorted(ctl_dir.glob("*.ctl"))
    if not ctl_files:
        print("❌ No .ctl files found in the directory.")
        return

    slurm_dir = Path(args.slurm_dir)
    slurm_dir.mkdir(exist_ok=True)

    total = len(ctl_files)
    cpu = args.cpu
    num_batches = math.ceil(total / cpu)

    print(f"🔍 Found {total} .ctl files, preparing {num_batches} SLURM batches with {cpu} jobs each...")

    for batch_id in range(num_batches):
        batch_files = ctl_files[batch_id * cpu : (batch_id + 1) * cpu]
        block_texts = []

        for ctl_file in batch_files:
            job_name = ctl_file.stem
            phy_file = job_name.replace("_model0", "").replace("_model1", "").replace("_model2", "") + ".phy"
            block = CODEML_BLOCK.format(
                job_name=job_name,
                ctl_file=ctl_file.name,
                phy_file=phy_file
            )
            block_texts.append(block.strip())

        batch_script = BATCH_TEMPLATE.format(
            batch_id=batch_id,
            ctl_dir=ctl_dir,
            codeml_blocks="\n\n".join(block_texts)
        )

        batch_script_path = slurm_dir / f"run_batch_{batch_id}.sh"
        batch_script_path.write_text(batch_script)
        subprocess.run(["sbatch", str(batch_script_path)])
        print(f"✅ Submitted batch {batch_id} with {len(batch_files)} job(s)")

if __name__ == "__main__":
    main()

