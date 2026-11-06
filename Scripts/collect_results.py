"""
collect_results.py

Collects the Boltz result files for one experiment into a single table.

Reads:  experiments/<experiment>/outputs/boltz_results_jobs/predictions/<id>/
            affinity_<id>.json              (only if affinity was on)
            confidence_<id>_model_0.json
Writes: experiments/<experiment>/results.csv

The <id> is the job YAML's file name without '.yaml', i.e. <proteinID>__<ligandID>,
so the two IDs are split back out into their own columns.

Example:
    python Scripts/collect_results.py --experiment <experiment>
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Collect the Boltz results into one table.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
outputs_dir = experiment_dir / "outputs"
# The table sits at the top of the experiment folder, beside reference/ and
# enrichment/, not inside outputs/ where Boltz writes its own files.
results_file = experiment_dir / "results.csv"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")


# ---------------------------------------------------------------------------
# 3. Find the predictions folder
# ---------------------------------------------------------------------------
# Boltz writes to outputs/boltz_results_<name of the jobs folder>/predictions/,
# where <name> is normally 'jobs'. rglob searches for it rather than assuming the
# name, so this still works if the folder is ever renamed.
predictions_dirs = sorted(outputs_dir.rglob("predictions")) if outputs_dir.exists() else []

if not predictions_dirs:
    sys.exit(f"ERROR: no predictions folder found under {outputs_dir}\n"
             "       Has the Boltz prediction step run yet?")
if len(predictions_dirs) > 1:
    found = "\n  ".join(str(d) for d in predictions_dirs)
    sys.exit(f"ERROR: more than one predictions folder found under {outputs_dir}:\n  {found}\n"
             "       Delete the outputs of the run you do not want and try again.")

predictions_dir = predictions_dirs[0]
print(f"Reading results from {predictions_dir}")


# ---------------------------------------------------------------------------
# 4. Read one row per prediction folder
# ---------------------------------------------------------------------------
# Each pair has its own folder, named after the job YAML (<proteinID>__<ligandID>).
rows = []
n_no_confidence = 0
n_no_affinity = 0

for pair_dir in sorted(predictions_dir.iterdir()):

    if not pair_dir.is_dir():
        continue

    pair_id = pair_dir.name

    # Split the ID back into the protein and the ligand. Both make_jobs.py and the
    # input checks forbid '__' inside an ID, so there is exactly one separator.
    if "__" in pair_id:
        protein_id, ligand_id = pair_id.split("__", 1)
    else:
        protein_id, ligand_id = pair_id, ""

    row = {"PROTEIN_ID": protein_id, "LIGAND_ID": ligand_id, "PAIR_ID": pair_id}

    # 4a. Affinity. The file is missing when the run used --no_affinity, which is fine.
    affinity_file = pair_dir / f"affinity_{pair_id}.json"
    if affinity_file.exists():
        with open(affinity_file) as f:
            affinity = json.load(f)
        # affinity_pred_value is log10(IC50 in uM): lower means a tighter binder.
        row["AFFINITY_PRED_VALUE"] = affinity.get("affinity_pred_value")
        row["AFFINITY_PROBABILITY_BINARY"] = affinity.get("affinity_probability_binary")
    else:
        n_no_affinity += 1

    # 4b. Confidence, for model 0 (the best-ranked pose).
    confidence_file = pair_dir / f"confidence_{pair_id}_model_0.json"
    if confidence_file.exists():
        with open(confidence_file) as f:
            confidence = json.load(f)
        # Only the overall scores are taken. The per-chain scores (chains_ptm and
        # pair_chains_iptm) are nested, so they are left in the JSON files.
        for key in ["confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm",
                    "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"]:
            row[key.upper()] = confidence.get(key)
    else:
        # No confidence file means this pair did not finish, so it is worth knowing about.
        print(f"  WARNING: no confidence file for '{pair_id}'")
        n_no_confidence += 1

    # 4c. Where the structure is, so the row can be traced back to a file.
    structure_file = pair_dir / f"{pair_id}_model_0.cif"
    if not structure_file.exists():
        structure_file = pair_dir / f"{pair_id}_model_0.pdb"
    row["STRUCTURE"] = str(structure_file.relative_to(project_dir)) if structure_file.exists() else ""

    rows.append(row)

if not rows:
    sys.exit(f"ERROR: {predictions_dir} has no prediction folders in it")


# ---------------------------------------------------------------------------
# 5. Write the table
# ---------------------------------------------------------------------------
results = pd.DataFrame(rows)

# Sort the tightest predicted binders first (a lower log10 IC50 is tighter).
# na_position="last" keeps any pair without an affinity value at the bottom.
if "AFFINITY_PRED_VALUE" in results.columns:
    results = results.sort_values("AFFINITY_PRED_VALUE", na_position="last")
else:
    results = results.sort_values(["PROTEIN_ID", "LIGAND_ID"])

results.to_csv(results_file, index=False, encoding="utf-8")


# ---------------------------------------------------------------------------
# 6. Report
# ---------------------------------------------------------------------------
print(f"Wrote {len(results)} rows to {results_file}")
if n_no_affinity:
    print(f"Note: {n_no_affinity} pair(s) had no affinity file (expected if the run used --no_affinity)")
if n_no_confidence:
    print(f"WARNING: {n_no_confidence} pair(s) had no confidence file: these predictions did not finish")
