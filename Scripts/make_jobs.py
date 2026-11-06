"""
make_jobs.py

Makes one Boltz YAML file for every protein x ligand pair in an experiment.

Reads:  experiments/<experiment>/inputs/proteins.csv   (ID, SEQUENCE, COFACTORS[, COVALENT_SITE])
        experiments/<experiment>/inputs/ligands.csv    (ID, SMILES[, WARHEAD])
        experiments/<experiment>/inputs/covalent_prepared.csv   (if this is a covalent experiment)
Writes: experiments/<experiment>/jobs/<proteinID>__<ligandID>.yaml

If covalent_prepared.csv is present, every job gets a bond constraint joining that atom to the protein's COVALENT_SITE. 

Example:
    python Scripts/make_jobs.py --experiment <experiment>
    python Scripts/make_jobs.py --experiment <experiment> --no_affinity
"""

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Make one Boltz YAML per protein x ligand pair.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
parser.add_argument("--no_affinity", action="store_true", help="Leave out the affinity prediction")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
# Using this (rather than the current folder) means the script works wherever you run it from.
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
inputs_dir = experiment_dir / "inputs"
jobs_dir = experiment_dir / "jobs"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")


# ---------------------------------------------------------------------------
# 3. Read the two CSV files
# ---------------------------------------------------------------------------
# encoding="utf-8-sig" removes the invisible marker some text editors put at the start of a file.
# dtype=str and keep_default_na=False keep every cell as text, with empty cells as "" (not NaN).
proteins = pd.read_csv(inputs_dir / "proteins.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)
ligands = pd.read_csv(inputs_dir / "ligands.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)

# Remove stray spaces around column names and values
proteins.columns = proteins.columns.str.strip()
ligands.columns = ligands.columns.str.strip()
proteins = proteins.apply(lambda column: column.str.strip())
ligands = ligands.apply(lambda column: column.str.strip())

print(f"Read {len(proteins)} proteins and {len(ligands)} ligands")

# If prepare_covalents.py has run and found warheads, it leaves this file behind.
# Its presence is what makes this a covalent experiment.
prepared_file = inputs_dir / "covalent_prepared.csv"
covalent = prepared_file.exists()
prepared_by_id = {}

if covalent:
    prepared = pd.read_csv(prepared_file, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    prepared.columns = prepared.columns.str.strip()
    prepared = prepared.apply(lambda column: column.str.strip())
    prepared_by_id = {row["ID"]: row for _, row in prepared.iterrows()}
    print(f"Covalent experiment: read {len(prepared_by_id)} prepared ligand(s)")


# ---------------------------------------------------------------------------
# 4. Check the inputs
# ---------------------------------------------------------------------------
# Every problem found is added to this list, and the script stops at the end
# of this section if there are any, so you can fix them all in one go.
problems = []

# 4a. The expected columns are there
for column in ["ID", "SEQUENCE", "COFACTORS"]:
    if column not in proteins.columns:
        problems.append(f"proteins.csv has no '{column}' column")
for column in ["ID", "SMILES"]:
    if column not in ligands.columns:
        problems.append(f"ligands.csv has no '{column}' column")
if problems:
    sys.exit("ERROR:\n  " + "\n  ".join(problems))

# 4b. IDs are filled in, unique, and safe to use in a file name
# ("__" is not allowed because it separates the protein and ligand IDs in the job name)
for file_name, table in [("proteins.csv", proteins), ("ligands.csv", ligands)]:
    for id_ in table["ID"]:
        if id_ == "":
            problems.append(f"{file_name}: a row has an empty ID")
        elif " " in id_ or "/" in id_ or "\\" in id_ or "__" in id_:
            problems.append(f"{file_name}: ID '{id_}' contains a space, a slash or '__'")
    duplicated = table["ID"][table["ID"].duplicated()].unique()
    for id_ in duplicated:
        problems.append(f"{file_name}: ID '{id_}' appears more than once")

# 4c. Protein sequences only contain the 20 standard amino acid letters
amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
for _, row in proteins.iterrows():
    if row["SEQUENCE"] == "":
        problems.append(f"proteins.csv: '{row['ID']}' has an empty SEQUENCE")
        continue
    for chain in row["SEQUENCE"].split(":"):
        bad_letters = set(chain.upper()) - amino_acids
        if chain == "":
            problems.append(f"proteins.csv: '{row['ID']}' has an empty chain (check the ':' separators)")
        elif bad_letters:
            problems.append(f"proteins.csv: '{row['ID']}' has letters that are not amino acids: {sorted(bad_letters)}")

# 4d. SMILES can be read by RDKit
for _, row in ligands.iterrows():
    if row["SMILES"] == "":
        problems.append(f"ligands.csv: '{row['ID']}' has an empty SMILES")
    elif Chem.MolFromSmiles(row["SMILES"]) is None:
        problems.append(f"ligands.csv: '{row['ID']}' has a SMILES that RDKit cannot read: {row['SMILES']}")

# 4e. Prevent make_jobs.py from being run on a stale or incomplete covalent_prepared.csv.
if covalent:
    for ligand_id in ligands["ID"]:
        if ligand_id not in prepared_by_id:
            problems.append(f"'{ligand_id}' is in ligands.csv but not in covalent_prepared.csv. "
                            "Re-run prepare_covalents.py.")
    if "COVALENT_SITE" not in proteins.columns:
        problems.append("this is a covalent experiment, but proteins.csv has no "
                        "'COVALENT_SITE' column")
    else:
        for _, row in proteins.iterrows():
            if row["COVALENT_SITE"] == "":
                problems.append(f"proteins.csv: '{row['ID']}' has an empty COVALENT_SITE")

if problems:
    sys.exit("ERROR: please fix these problems in the input files:\n  " + "\n  ".join(problems))

print("Input checks passed")


# ---------------------------------------------------------------------------
# 5. Empty the jobs folder
# ---------------------------------------------------------------------------
# Boltz runs every .yaml file in this folder, so old files from an earlier run
# (e.g. a ligand you have since removed from the CSV) are deleted first.
jobs_dir.mkdir(exist_ok=True)
for old_file in jobs_dir.glob("*.yaml"):
    old_file.unlink()


# ---------------------------------------------------------------------------
# 6. Write one YAML per protein x ligand pair
# ---------------------------------------------------------------------------
chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
missing_msas = set()

for _, protein in proteins.iterrows():

    # The MSA cache folder for this protein is named after a hash of the full SEQUENCE text,
    # and each chain's MSA file is named after a hash of that chain's sequence.
    complex_hash = hashlib.sha256(protein["SEQUENCE"].encode()).hexdigest()[:16]
    chains = [chain.upper() for chain in protein["SEQUENCE"].split(":")]

    cofactors = []
    if protein["COFACTORS"] != "":
        cofactors = [code.strip().upper() for code in protein["COFACTORS"].split(";") if code.strip() != ""]

    # Where a covalent ligand attaches, e.g. "A:C16" -> chain A, residue 16.
    # prepare_covalents.py has already checked the chain exists, the position is
    # inside it, and the residue there really is the letter written down, so this
    # only has to pull the two numbers back out.
    site_chain = None
    site_position = None
    if covalent:
        chain_part, residue_part = protein["COVALENT_SITE"].split(":", 1)
        site_chain = chain_part.strip().upper()
        site_position = int(residue_part.strip()[1:])   # drop the residue letter

    for _, ligand in ligands.iterrows():

        sequences = []      # the list of chains that goes into the YAML
        letter_index = 0    # which chain letter to use next

        # Protein chains first: A, B, ...
        for chain in chains:
            chain_hash = hashlib.sha256(chain.encode()).hexdigest()[:16]
            msa_path = f"cache_msa/{complex_hash}/{chain_hash}.csv"   # relative to the project folder
            if not (project_dir / msa_path).exists():
                missing_msas.add(msa_path)
            sequences.append({"protein": {
                "id": chain_letters[letter_index],
                "sequence": chain,
                "msa": msa_path,
            }})
            letter_index += 1

        # Then the cofactors, given as CCD codes (e.g. ZN)
        for code in cofactors:
            sequences.append({"ligand": {
                "id": chain_letters[letter_index],
                "ccd": code,
            }})
            letter_index += 1

        # The screened ligand goes last. In a covalent experiment it goes in as its
        # bound (open-valence) form, which is what the bond constraint attaches to.
        ligand_letter = chain_letters[letter_index]
        if covalent:
            ligand_smiles = prepared_by_id[ligand["ID"]]["BOUND_SMILES"]
        else:
            ligand_smiles = ligand["SMILES"]
        sequences.append({"ligand": {
            "id": ligand_letter,
            "smiles": ligand_smiles,
        }})

        job = {"version": 1, "sequences": sequences}

        # The covalent bond itself. A SMILES ligand is one residue, numbered 1, and
        # its atoms are named by Boltz's own rule - see PLAN.md for both facts.
        if covalent:
            prepared_row = prepared_by_id[ligand["ID"]]
            job["constraints"] = [{"bond": {
                "atom1": [site_chain, site_position, prepared_row["TARGET_ATOM"]],
                "atom2": [ligand_letter, 1, prepared_row["ATTACHMENT_ATOM"]],
            }}]

        if not args.no_affinity:
            job["properties"] = [{"affinity": {"binder": ligand_letter}}]

        # width=100000 stops long sequences being split across lines
        job_file = jobs_dir / f"{protein['ID']}__{ligand['ID']}.yaml"
        with open(job_file, "w") as f:
            yaml.safe_dump(job, f, sort_keys=False, width=100000)


# ---------------------------------------------------------------------------
# 7. Report
# ---------------------------------------------------------------------------
n_jobs = len(proteins) * len(ligands)
print(f"Wrote {n_jobs} job files to {jobs_dir}")
if covalent:
    print("These are covalent jobs: each one bonds the ligand to its protein's COVALENT_SITE")
if args.no_affinity:
    print("Affinity prediction is OFF for these jobs")

if missing_msas:
    # This is expected on the laptop or before the MSA step has run, so it is a note, not an error.
    print(f"Note: {len(missing_msas)} MSA file(s) are not in cache_msa/ yet. "
          "Boltz will fail on these jobs unless the MSA step runs first.")
