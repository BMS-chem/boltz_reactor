"""
calculate_enrichment.py

Scores one experiment against a list of known actives and inactives - useful
for benchmarking or validating a model against experimental data.

Reads:  experiments/<experiment>/results.csv
        experiments/<experiment>/reference/actives_inactives.txt
            ID,LABEL          (header row; LABEL is 'active' or 'inactive')
            The ID must match the LIGAND_ID column of results.csv.

Writes: experiments/<experiment>/enrichment/enrichment.txt
        experiments/<experiment>/enrichment/roc_affinity_<protein>.png
        experiments/<experiment>/enrichment/results_annotated.csv
            (results.csv annotated with active/inactive labels)

The score ranked is AFFINITY_PRED_VALUE, which is log10(IC50 uM).

The AUC is a plain ROC AUC over predicted/labelled compounds with no weighting.

Example:
    python Scripts/calculate_enrichment.py --experiment 001_eLife_Sigma2_test
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

# Fixing matplotlib import - probably something wrong with my environment? need to fix
_dll_dir = os.path.join(sys.prefix, "Library", "bin")
if os.name == "nt" and os.path.isdir(_dll_dir):
    os.environ["PATH"] = _dll_dir + os.pathsep + os.environ.get("PATH", "")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# The two calculations, kept as functions to keep the main flow readable
# ---------------------------------------------------------------------------
def roc_auc(scores, is_active):
    """Plain ROC AUC. Higher score = predicted active. Ties get the average rank."""
    ranks = scores.rank(method="average")
    n_active = int(is_active.sum())
    n_inactive = int((~is_active).sum())
    # Mann-Whitney U: the chance a random active outranks a random inactive.
    sum_active_ranks = ranks[is_active].sum()
    return (sum_active_ranks - n_active * (n_active + 1) / 2) / (n_active * n_inactive)


def roc_curve(scores, is_active):
    """False positive rate and true positive rate, for plotting."""
    # Walk down the compounds from the best score to the worst. Every compound
    # with the same score must be counted at once, or ties make a jagged curve.
    table = pd.DataFrame({"score": scores, "active": is_active})
    table = table.sort_values("score", ascending=False)

    counts = table.groupby("score", sort=False)["active"].agg(["sum", "count"])
    true_positives = counts["sum"].cumsum()
    false_positives = (counts["count"] - counts["sum"]).cumsum()

    n_active = int(is_active.sum())
    n_inactive = int((~is_active).sum())

    # Start the curve at the origin (nothing called active yet).
    fpr = [0.0] + list(false_positives / n_inactive)
    tpr = [0.0] + list(true_positives / n_active)
    return fpr, tpr


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Calculate ROC AUC for known actives against inactives.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Work out the folder paths
# ---------------------------------------------------------------------------
# project folder is one level up
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
results_file = experiment_dir / "results.csv"
labels_file = experiment_dir / "reference" / "actives_inactives.txt" # The label file is probably going to be a hand-edited CSV-style list
enrichment_dir = experiment_dir / "enrichment"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")
if not results_file.exists():
    sys.exit(f"ERROR: results table not found: {results_file}\n"
             "       Run collect_results.py first.")
if not labels_file.exists():
    sys.exit(f"ERROR: label file not found: {labels_file}\n"
             "       It needs a header row 'ID,LABEL' and one row per compound,\n"
             "       where LABEL is 'active' or 'inactive' and ID matches LIGAND_ID.")

enrichment_dir.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# 3. Read the results and the labels
# ---------------------------------------------------------------------------
results = pd.read_csv(results_file, encoding="utf-8-sig")

for column in ["PROTEIN_ID", "LIGAND_ID", "AFFINITY_PRED_VALUE"]:
    if column not in results.columns:
        sys.exit(f"ERROR: {results_file} has no '{column}' column.\n"
                 "       Was this experiment run with affinity prediction on?")

results["LIGAND_ID"] = results["LIGAND_ID"].astype(str).str.strip()

# The label file is probably going to be a hand-edited CSV-style list
labels = pd.read_csv(labels_file, encoding="utf-8-sig")

if list(labels.columns[:2]) != ["ID", "LABEL"]:
    sys.exit(f"ERROR: {labels_file} must start with the columns 'ID,LABEL', "
             f"not {list(labels.columns)}")

labels["ID"] = labels["ID"].astype(str).str.strip()
labels["LABEL"] = labels["LABEL"].astype(str).str.strip().str.lower()

unexpected = sorted(set(labels["LABEL"]) - {"active", "inactive"})
if unexpected:
    sys.exit(f"ERROR: {labels_file} has LABEL values other than 'active' or 'inactive': "
             f"{unexpected}")

if labels["ID"].duplicated().any():
    repeated = sorted(labels.loc[labels["ID"].duplicated(), "ID"].unique())
    sys.exit(f"ERROR: {labels_file} lists the same ID more than once: {repeated}")

print(f"Labels:  {len(labels)} compounds "
      f"({(labels['LABEL'] == 'active').sum()} active, "
      f"{(labels['LABEL'] == 'inactive').sum()} inactive)")
print(f"Results: {len(results)} rows, {results['PROTEIN_ID'].nunique()} protein(s)")


# ---------------------------------------------------------------------------
# 4. Report what the two files do not have in common
# ---------------------------------------------------------------------------
# A run that was cut short or resumed is a normal state here, so this is a
# report rather than an error. The AUCs below use the overlap.
predicted_ids = set(results["LIGAND_ID"])
labelled_ids = set(labels["ID"])

missing_predictions = sorted(labelled_ids - predicted_ids)
unlabelled_ligands = sorted(predicted_ids - labelled_ids)

if missing_predictions:
    print(f"Note: {len(missing_predictions)} labelled compound(s) have no prediction, "
          f"e.g. {missing_predictions[:3]}")
if unlabelled_ligands:
    print(f"Note: {len(unlabelled_ligands)} predicted ligand(s) are not in the label file, "
          f"e.g. {unlabelled_ligands[:3]}")
if not labelled_ids & predicted_ids:
    sys.exit("ERROR: no labelled compound has a prediction. Do the IDs in the label "
             "file match the LIGAND_ID column of results.csv?")


# ---------------------------------------------------------------------------
# 5. Join, then score each protein in turn
# ---------------------------------------------------------------------------
# The pipeline pairs every protein with every ligand, so an experiment can hold
# more than one protein. Each gets its own AUC and its own plot.
scored = results.merge(labels, left_on="LIGAND_ID", right_on="ID", how="inner")
scored["IS_ACTIVE"] = scored["LABEL"] == "active"


# ---------------------------------------------------------------------------
# 5a. The annotated results table, for looking at individual compounds
# ---------------------------------------------------------------------------
# Every prediction is kept, including any ligand the label file does not cover,
# which gets an empty LABEL. Sorted the way results.csv is: tightest first.
annotated = results.merge(labels, left_on="LIGAND_ID", right_on="ID", how="left")
annotated = annotated.drop(columns=["ID"])
annotated["LABEL"] = annotated["LABEL"].fillna("")
annotated = annotated.sort_values("AFFINITY_PRED_VALUE", na_position="last")

annotated_file = enrichment_dir / "results_annotated.csv"
annotated.to_csv(annotated_file, index=False, encoding="utf-8")
print(f"Wrote {len(annotated)} rows to {annotated_file.name}")

report_lines = []
report_lines.append(f"Enrichment for experiment: {args.experiment}")
report_lines.append(f"Results file: {results_file.relative_to(project_dir)}")
report_lines.append(f"Label file:   {labels_file.relative_to(project_dir)}")
report_lines.append("Score ranked: AFFINITY_PRED_VALUE, negated (lower is a tighter binder)")
report_lines.append("")
report_lines.append(f"Labelled compounds:            {len(labels)}")
report_lines.append(f"Labelled without a prediction: {len(missing_predictions)}")
report_lines.append(f"Predicted without a label:     {len(unlabelled_ligands)}")
report_lines.append("")

for protein_id, protein_rows in scored.groupby("PROTEIN_ID"):

    # Drop any pair with no affinity value (a missing affinity file).
    usable = protein_rows[protein_rows["AFFINITY_PRED_VALUE"].notna()]
    n_dropped = len(protein_rows) - len(usable)

    is_active = usable["IS_ACTIVE"]
    n_active = int(is_active.sum())
    n_inactive = int((~is_active).sum())

    if n_active == 0 or n_inactive == 0:
        message = (f"{protein_id}: skipped, needs both actives and inactives "
                   f"(have {n_active} active, {n_inactive} inactive)")
        print(message)
        report_lines.append(message)
        continue

    # Negate, so that a higher score means a tighter predicted binder.
    ranking_score = -usable["AFFINITY_PRED_VALUE"]
    auc = roc_auc(ranking_score, is_active)
    fpr, tpr = roc_curve(ranking_score, is_active)

    line = (f"{protein_id}: AUC {auc:.3f}  "
            f"(n = {len(usable)}: {n_active} active, {n_inactive} inactive"
            + (f", {n_dropped} dropped with no affinity value" if n_dropped else "") + ")")
    print(line)
    report_lines.append(line)

    # 5a. The plot for this protein.
    figure, axes = plt.subplots(figsize=(5, 5))
    axes.plot(fpr, tpr, linewidth=2, label=f"AUC = {auc:.3f}")
    axes.plot([0, 1], [0, 1], linestyle="--", linewidth=1, color="grey", label="random")
    axes.set_xlabel("False positive rate")
    axes.set_ylabel("True positive rate")
    axes.set_title(f"{protein_id}\nPredicted affinity", fontsize=10)
    axes.set_xlim(0, 1)
    axes.set_ylim(0, 1)
    axes.legend(loc="lower right")
    figure.tight_layout()

    # Protein IDs are checked for spaces and slashes upstream, so they are safe
    # to put in a file name.
    plot_file = enrichment_dir / f"roc_affinity_{protein_id}.png"
    figure.savefig(plot_file, dpi=150)
    plt.close(figure)
    print(f"  wrote {plot_file.name}")


# ---------------------------------------------------------------------------
# 6. Write the numbers alongside the plots
# ---------------------------------------------------------------------------
report_file = enrichment_dir / "enrichment.txt"
report_file.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

print(f"\nWrote {report_file}")
