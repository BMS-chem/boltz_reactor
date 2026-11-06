# Boltz batch pipeline

Takes a list of proteins and a list of ligands, makes one [Boltz-2](https://github.com/jwohlwend/boltz)
job for every protein x ligand pair, runs them on a Slurm cluster, and collects
the results into one table. 
Covalent ligands are written as regular SMILES then reacted via SMARTS strings.
Enrichment of actives/inactives during benchmarking/validation can also be calculated.


## Requirements

- Boltz-2 (v2.2.0) and a CUDA GPU on a Slurm cluster
- Python 3.12 with `pandas`, `pyyaml`, `rdkit`, `matplotlib`

```
conda create -n boltz_pipeline_env -c conda-forge python=3.12 pandas pyyaml rdkit matplotlib
```

The Python scripts run anywhere (including a laptop with no GPU) so inputs can be
checked, and results can be analysed locally; the cluster is only needed for the actual Boltz predictions.

## Workflow

```
submit_boltz.sh                  login node. Checks the inputs, then calls sbatch.
  └── run_pipeline.sh            GPU node. Runs the steps below in order.
        1.  make_msa_jobs.py     writes jobs for protein sequences without a cached MSA
        2.  boltz predict        fetches any missing MSAs
        3.  cache_msas.py        files them in cache_msa/
        3b. prepare_covalents.py applies the warhead reaction (for covalent experiments)
        4.  make_jobs.py         writes jobs for protein/ligand pairs, using cached MSAs
        5.  boltz predict        runs the predictions
        6.  collect_results.py   collates the prediction outputs into results.csv
```

Scripts are run from the project root because the job YAMLs refer to cached MSAs by relative path.
MSAs are cached per complex in `cache_msa/`.

## Layout

```
Scripts/          the pipeline
docs/             notes and write-ups
cache_msa/        protein MSAs, named by sequence hash, reused across experiments
cache_weights/    Boltz model weights and CCD data
experiments/<name>/
  inputs/         proteins.csv and ligands.csv, written by hand
  jobs/           one YAML per protein x ligand pair
  outputs/        structures and scores
  logs/           Slurm logs for this experiment
  reference/      source material for the experiment, not read by the pipeline
  enrichment/     benchmarking output (if relevant)
  results.csv     the summary table from Boltz predictions
```

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
| `WARHEAD` | covalent experiments only. So far includes `acrylamide`, `chloroacetamide` or `CF2_acrylamide` |

`make_jobs.py` checks these before anything is submitted and stops if any input string isn't formatted correctly.

## Running it

Make a new experiment folder with empty input templates:

```
python Scripts/new_experiment.py --experiment <experiment>
```

Fill in the protein and ligand CSVs. Optionally, check them locally by running the job writing script:

```
python Scripts/make_jobs.py --experiment <experiment>
```

Then, log in to the cluster, navigate to the project directory, and run:

```
./Scripts/submit_boltz.sh <experiment>
```

Options: `--no_affinity` (structure only), `--override` (redo pairs that already
have output), `--time 1-00:00:00` (Slurm time limit).

The Slurm partition, resources and notification address are set at the top of
`Scripts/submit_boltz.sh` and will need changing for different clusters. Scripts may
need `chmod +x Scripts/*.sh` once to run.

## Covalent ligands

The ligands are given as unreacted SMILES and categorised in the `WARHEAD` column.
New warheads can be written into the dictionary at the top of `prepare_covalents.py`.
This script applies the reaction SMARTS, opens one valence on the attachment atom, 
names the atom the way Boltz names SMILES ligand atoms, and writes 
`inputs/covalent_prepared.csv`. `make_jobs.py` will then also write a `bond`
constraint between that atom and the protein's `COVALENT_SITE`, which itself is 
checked against the sequence.

Experiments can't mix covalent and non-covalent jobs.

## Benchmarking against known actives

Where an experiment consists of ligands that can be categorised into known
active and inactives based on experimental data, `calculate_enrichment.py`
reports how well the predicted affinity separates them. It runs after
`collect_results.py`, and you need to have provided a list in 
`reference/actives_inactives.txt`:

```
ID,LABEL
Ligand1,active
Ligand2,inactive
```

`ID` must match the `LIGAND_ID` column of `results.csv`.

```
python Scripts/calculate_enrichment.py --experiment <virtual_vs_experimental_screen>
```

This writes `enrichment/enrichment.txt`, a ROC plot per protein, and
`enrichment/results_annotated.csv`.

The pipeline was checked this way against a published virtual screen, giving a
ROC AUC of 0.848 where the paper reports 0.838 for the same compounds and model:
see [docs/sigma2-benchmark.md](docs/sigma2-benchmark.md).

## Still figuring out:

- **Affinity on covalent complexes.** Jack told me that Boltz can predict covalent
binding affinity after all (not just pose prediction) - analogous to Ki; it doesn't
factor in Kinact - is this just "affinity_pred_value" on a covalent job, or am I missing
something?

## Licence

MIT, see [LICENSE](LICENSE).
