"""
make_msa_jobs.py

Works out which proteins do not have a cached MSA yet, and writes one protein-only
Boltz job for each of them. Those jobs are then run with --use_msa_server, which
downloads the MSAs (cache_msas.py puts them into the cache afterwards).

Reads:  experiments/<experiment>/inputs/proteins.csv
Writes: experiments/<experiment>/msa_jobs/<complex_hash>.yaml   (one per uncached protein)

If every protein is already cached, it writes nothing and the MSA fetch is skipped.

Example:
    python Scripts/make_msa_jobs.py --experiment <experiment>
"""

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Write a protein-only Boltz job for each uncached MSA.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Work out the folder paths
# ---------------------------------------------------------------------------
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
inputs_dir = experiment_dir / "inputs"
msa_jobs_dir = experiment_dir / "msa_jobs"
cache_dir = project_dir / "cache_msa"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")


# ---------------------------------------------------------------------------
# 3. Read the protein list
# ---------------------------------------------------------------------------
proteins = pd.read_csv(inputs_dir / "proteins.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)
proteins.columns = proteins.columns.str.strip()
proteins = proteins.apply(lambda column: column.str.strip())


# ---------------------------------------------------------------------------
# 4. Empty the msa_jobs folder
# ---------------------------------------------------------------------------
# Boltz runs every .yaml in this folder, so leftovers from a previous run
# (already cached by now) must not be fetched again.
msa_jobs_dir.mkdir(exist_ok=True)
for old_file in msa_jobs_dir.glob("*.yaml"):
    old_file.unlink()


# ---------------------------------------------------------------------------
# 5. Decide which proteins still need an MSA
# ---------------------------------------------------------------------------
# One MSA job per distinct SEQUENCE text, not per row: two rows with the same
# sequence share one MSA. seen_complexes stops us writing the same job twice.
seen_complexes = set()
n_written = 0
n_cached = 0

for _, protein in proteins.iterrows():

    sequence_text = protein["SEQUENCE"]
    if sequence_text in seen_complexes:
        continue
    seen_complexes.add(sequence_text)

    # Work out what the cache folder for this complex would be
    complex_hash = hashlib.sha256(sequence_text.encode()).hexdigest()[:16]
    chains = [chain.upper() for chain in sequence_text.split(":")]

    # Work out whether every chain already has a file
    distinct_chains = list(dict.fromkeys(chains))
    all_present = True
    for chain in distinct_chains:
        chain_hash = hashlib.sha256(chain.encode()).hexdigest()[:16]
        if not (cache_dir / complex_hash / f"{chain_hash}.csv").exists():
            all_present = False

    if all_present:
        print(f"  cached:   {protein['ID']}")
        n_cached += 1
        continue

    # If not already cached, write a protein-only job. There is no 'msa' key and no ligand:
    # Boltz will fetch the MSA from the server, and we throw the structure away.
    # The file is named after the complex hash so cache_msas.py can tell which
    # complex each downloaded MSA belongs to.
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    sequences = []
    for i, chain in enumerate(chains):
        sequences.append({"protein": {"id": chain_letters[i], "sequence": chain}})

    job = {"version": 1, "sequences": sequences}
    with open(msa_jobs_dir / f"{complex_hash}.yaml", "w") as f:
        yaml.safe_dump(job, f, sort_keys=False, width=100000)

    print(f"  to fetch: {protein['ID']}  ({len(chains)} chain(s) -> {complex_hash}.yaml)")
    n_written += 1


# ---------------------------------------------------------------------------
# 6. Report
# ---------------------------------------------------------------------------
print(f"{n_cached} protein(s) already cached, {n_written} to fetch")
if n_written == 0:
    print("Nothing to fetch.")
