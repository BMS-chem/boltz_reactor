"""
cache_msas.py

Copies new MSAs into the shared cache so that later runs with the same sequence don't have to fetch them again.

Reads:  experiments/<experiment>/msa_outputs/.../msa/*.csv   (written by the first boltz run)
        experiments/<experiment>/inputs/proteins.csv         (to check protein sequence hashes)
Writes: cache_msa/<complex_hash>/<chain_hash>.csv

Each downloaded file is named <complex_hash>_<number>.csv, where the number is
Boltz's own chain numbering. This script reads the first sequence in the file
and uses that to verify which chain the file belongs to.

Example:
    python Scripts/cache_msas.py --experiment <experiment>
"""

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Copy downloaded MSAs into cache_msa/.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Work out the folder paths
# ---------------------------------------------------------------------------
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
inputs_dir = experiment_dir / "inputs"
msa_outputs_dir = experiment_dir / "msa_outputs"
cache_dir = project_dir / "cache_msa"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")


# ---------------------------------------------------------------------------
# 3. Find the downloaded MSA files
# ---------------------------------------------------------------------------
# Boltz puts them in a folder called 'msa' somewhere inside the output folder
# (currently msa_outputs/boltz_results_msa_jobs/msa/). rglob searches all the way
# down, so this keeps working if Boltz changes the folder names around it.
downloaded = []
if msa_outputs_dir.exists():
    for msa_folder in msa_outputs_dir.rglob("msa"):
        if msa_folder.is_dir():
            downloaded.extend(sorted(msa_folder.glob("*.csv")))

if not downloaded:
    print("No downloaded MSA files found - nothing to add to the cache.")
    sys.exit(0)

print(f"Found {len(downloaded)} downloaded MSA file(s)")


# ---------------------------------------------------------------------------
# 4. Work out which chain sequence each hash belongs to
# ---------------------------------------------------------------------------
# Build a lookup of every chain sequence in this experiment, so that each
# downloaded file can be matched to the chain it actually contains.
# Key: the chain sequence. Value: its hash.
proteins = pd.read_csv(inputs_dir / "proteins.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)
proteins.columns = proteins.columns.str.strip()
proteins = proteins.apply(lambda column: column.str.strip())

chain_hash_of = {}
for _, protein in proteins.iterrows():
    for chain in protein["SEQUENCE"].split(":"):
        chain = chain.upper()
        chain_hash_of[chain] = hashlib.sha256(chain.encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# 5. Copy each file into the cache
# ---------------------------------------------------------------------------
n_copied = 0
problems = []

for msa_file in downloaded:

    # Everything before the "_" in the filename is the complex hash.
    complex_hash = msa_file.stem.rsplit("_", 1)[0]

    # Read the first sequence in the file. The columns are 'key' and 'sequence',
    # and the first row is the protein we asked about.
    table = pd.read_csv(msa_file, dtype=str, keep_default_na=False)
    if "sequence" not in table.columns or len(table) == 0:
        problems.append(f"{msa_file.name}: not in the expected 'key,sequence' format")
        continue

    # Gaps and lower-case letters mark insertions in an alignment; remove them
    # to get back the plain protein sequence.
    query = table["sequence"].iloc[0].upper().replace("-", "")

    if query not in chain_hash_of:
        problems.append(
            f"{msa_file.name}: its first sequence does not match any chain in proteins.csv "
            f"(starts with {query[:30]}...)"
        )
        continue

    destination = cache_dir / complex_hash / f"{chain_hash_of[query]}.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(msa_file, destination)
    print(f"  cached: {msa_file.name}  ->  cache_msa/{complex_hash}/{destination.name}")
    n_copied += 1


# ---------------------------------------------------------------------------
# 6. Report
# ---------------------------------------------------------------------------
print(f"Added {n_copied} MSA file(s) to the cache")

if problems:
    # Stop here: carrying on would run predictions with missing MSAs.
    sys.exit("ERROR: some downloaded MSAs could not be matched to a protein:\n  " + "\n  ".join(problems))
