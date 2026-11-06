"""
new_experiment.py

Sets up a new experiment folder with template input CSVs.

Makes:  experiments/<experiment>/inputs/proteins.csv   (header row only)
        experiments/<experiment>/inputs/ligands.csv    (header row only)
        experiments/<experiment>/jobs/  outputs/  logs/

Existing files are never overwritten, so running this twice is harmless.

Example:
    python Scripts/new_experiment.py --experiment my_screen
"""

import argparse
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Set up a new experiment folder.")
parser.add_argument("--experiment", required=True, help="Name for the new folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Check the name is usable as a folder name
# ---------------------------------------------------------------------------
# The name becomes a folder name and part of the Slurm job name, so spaces and
# slashes are rejected here rather than causing a confusing error later on.
name = args.experiment.strip()

if name == "":
    sys.exit("ERROR: the experiment name is empty")
for bad in [" ", "/", "\\"]:
    if bad in name:
        sys.exit(f"ERROR: the experiment name must not contain '{bad}': {name}")


# ---------------------------------------------------------------------------
# 3. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / name
inputs_dir = experiment_dir / "inputs"

if experiment_dir.exists():
    print(f"Note: {experiment_dir} already exists. Only the missing parts will be made.")


# ---------------------------------------------------------------------------
# 4. Make the folders
# ---------------------------------------------------------------------------
# msa_jobs/ and msa_outputs/ are left out on purpose: run_pipeline.sh makes them
# when it needs them and deletes them again at step 3, so an experiment that has
# not been run should not have them. exist_ok=True leaves an existing folder alone.
for folder in ["inputs", "jobs", "outputs", "logs"]:
    (experiment_dir / folder).mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# 5. Write the two template CSVs
# ---------------------------------------------------------------------------
# Just the header row. A file that already exists is left untouched, so this
# cannot wipe out input you have already typed in.
templates = {
    "proteins.csv": "ID,SEQUENCE,COFACTORS\n",
    "ligands.csv": "ID,SMILES\n",
}

for file_name, header in templates.items():
    file_path = inputs_dir / file_name
    if file_path.exists():
        print(f"  kept:    {file_name} (already there, left as it is)")
    else:
        with open(file_path, "w", encoding="utf-8", newline="") as f:
            f.write(header)
        print(f"  created: {file_name}")


# ---------------------------------------------------------------------------
# 6. Report, and say what to do next
# ---------------------------------------------------------------------------
print(f"\nExperiment folder ready: {experiment_dir}")
print("\nNext:")
print(f"  1. Fill in {inputs_dir / 'proteins.csv'}")
print("     ID,SEQUENCE,COFACTORS  - one row per protein.")
print("     SEQUENCE separates chains with ':'. COFACTORS holds CCD codes")
print("     separated by ';' (e.g. ZN;MG) and may be left empty.")
print(f"  2. Fill in {inputs_dir / 'ligands.csv'}")
print("     ID,SMILES  - one row per ligand.")
print("     IDs must be unique and must not contain spaces, slashes or '__'.")
print(f"  3. On the cluster, submit it:  ./Scripts/submit_boltz.sh {name}")
