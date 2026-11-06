# Boltz batch pipeline

Takes a list of proteins and a list of ligands, makes one [Boltz-2](https://github.com/jwohlwend/boltz)
job for every protein x ligand pair, runs them on a Slurm cluster, and collects
the results into one table. It also supports covalent ligands, and can score an
experiment against known actives and inactives.

It is a set of small, plain Python and shell scripts rather than a framework.
There is no package to install and nothing imports anything else.

## Requirements

- Boltz-2 (developed against v2.2.0) and a CUDA GPU, on a Slurm cluster
- Python 3.12 with `pandas`, `pyyaml` and `rdkit` (plus `matplotlib`, for the
  enrichment plots only)

```
conda create -n boltz_pipeline_env -c conda-forge python=3.12 pandas pyyaml rdkit matplotlib
```

The Python scripts run anywhere, including a laptop with no GPU, so inputs can be
checked and results analysed locally. Only Boltz itself needs the cluster.

## How it fits together

```
submit_boltz.sh          login node. Checks the inputs, then calls sbatch.
  └── run_pipeline.sh    GPU node. Runs the steps below in order.
        1.  make_msa_jobs.py     for protein sequences without a cached MSA
        2.  boltz predict        fetches any missing MSAs
        3.  cache_msas.py        files them in cache_msa/
        3b. prepare_covalents.py applies the warhead reaction (covalent only)
        4.  make_jobs.py         checks the CSVs, writes one YAML per pair
        5.  boltz predict        runs the predictions
        6.  collect_results.py   gathers the JSON output into results.csv
```

Boltz is always started from the project root, because the job YAMLs refer to
cached MSAs by relative path.

MSAs are cached per complex in `cache_msa/`, so adding ligands to an existing
screen costs nothing extra: cached MSAs are reused and Boltz skips pairs that
already have output.

## Folders

```
Scripts/          the pipeline
docs/             design notes and the benchmark write-up
cache_msa/        protein MSAs, named by sequence hash, reused across experiments
cache_weights/    Boltz model weights and CCD data
experiments/<name>/
  inputs/         proteins.csv and ligands.csv, written by hand
  jobs/           one YAML per protein x ligand pair
  outputs/        structures and scores
  logs/           Slurm logs for this experiment
  reference/      source material for this experiment, not read by the pipeline
  enrichment/     benchmarking output, if the ligands have known labels
  results.csv     the summary table
```

`cache_weights/`, `cache_msa/` and the per-experiment `jobs/`, `outputs/` and
`logs/` are not in version control.

## Input format

`inputs/proteins.csv`

| Column | |
|---|---|
| `ID` | used in output file names; no spaces, slashes or `__` |
| `SEQUENCE` | one-letter amino acids; separate chains of a complex with `:` |
| `COFACTORS` | CCD codes separated by `;` (e.g. `ZN;MG`), or empty |
| `COVALENT_SITE` | covalent experiments only: where the ligand attaches, as `<chain>:<residue><position>`, e.g. `A:C123` |

`inputs/ligands.csv`

| Column | |
|---|---|
| `ID` | as above |
| `SMILES` | the ligand, written unreacted |
| `WARHEAD` | covalent experiments only: `acrylamide`, `chloroacetamide` or `CF2_acrylamide` |

`make_jobs.py` checks these before anything is submitted, and stops on duplicate
or malformed IDs, empty cells, non-amino-acid letters, or SMILES RDKit cannot
parse.

## Running it

Make a new experiment folder with empty input templates:

```
python Scripts/new_experiment.py --experiment my_screen
```

Fill in the two CSVs, then check them locally by writing the YAMLs without
running anything:

```
python Scripts/make_jobs.py --experiment my_screen
```

Then, on the cluster:

```
./Scripts/submit_boltz.sh my_screen
```

Options: `--no_affinity` (structure only), `--override` (redo pairs that already
have output), `--time 1-00:00:00` (Slurm time limit).

The Slurm partition, resources and notification address are set at the top of
`Scripts/submit_boltz.sh` and will need changing for your cluster. The scripts
need `chmod +x Scripts/*.sh` once.

## Covalent ligands

Give the ligand as its unreacted SMILES and name a `WARHEAD`.
`prepare_covalents.py` applies that warhead's reaction SMARTS, opens one valence
on the attachment atom, names the atom the way Boltz names SMILES ligand atoms,
and writes `inputs/covalent_prepared.csv`. `make_jobs.py` then writes a `bond`
constraint from the protein's `COVALENT_SITE` to that atom.

`COVALENT_SITE` is checked against the sequence, so a wrong residue or an
off-by-one fails before any GPU time is spent. Add a warhead by adding an entry
to the `WARHEADS` dictionary at the top of `prepare_covalents.py`.

An experiment is wholly covalent or wholly not. Boltz crashes if one job folder
mixes jobs that ask for affinity with jobs that do not; `docs/design.md` has the
details.

## Benchmarking against known actives

Where an experiment's ligands have measured labels, `calculate_enrichment.py`
reports how well the predicted affinity separates them. It runs after
`collect_results.py`, and needs `reference/actives_inactives.txt`:

```
ID,LABEL
ZINC000657934399,active
ZINC000466933503,inactive
```

`ID` must match the `LIGAND_ID` column of `results.csv`.

```
python Scripts/calculate_enrichment.py --experiment 001_eLife_Sigma2_test
```

This writes `enrichment/enrichment.txt`, a ROC plot per protein, and
`enrichment/results_annotated.csv` (the results table with the labels joined on,
tightest predicted binder first).

The pipeline was checked this way against a published virtual screen, giving a
ROC AUC of 0.848 where the paper reports 0.838 for the same compounds and model:
see [docs/sigma2-benchmark.md](docs/sigma2-benchmark.md).

## What isn't established

- **Affinity on covalent complexes.** Imposing a bond changed the predicted
  binder probability far more than the predicted potency, and the two ensemble
  members disagreed across the decision boundary. Do not rank a covalent screen
  on affinity until it has been tested against known covalent binders and
  non-binders, as the sigma-2 set was used above.
- **Bond geometry.** The covalent bond forms at the correct atom, but its length
  comes out compressed and variable (1.27-1.67 A observed, against about 1.81 A
  for C-S). Use it only as a yes/no check that the constraint applied, with a
  cutoff near 2.5 A. Never rank on it.
- The benchmark above covers one target, which Boltz-2 had seen in training.

`docs/design.md` records the design decisions and the Boltz behaviours they rely
on, each with how it was checked.

## Licence

MIT, see [LICENSE](LICENSE).
