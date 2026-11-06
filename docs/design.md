# Pipeline plan (v1)

Take lists of proteins and ligands, makes one Boltz-2 job per protein x ligand pair, runs them on HPC (NEMO), and collects the results into one table.

## Decisions

- **Pairs:** every protein x every ligand.
- **Inputs:** `SEQUENCE` separates protein chains with `:` (`SEQA:SEQB`, or `SEQA:SEQA` for a homodimer). `COFACTORS` holds CCD codes separated by `;` (`ZN;MG`) and may be empty.
- **Input checks:** `make_jobs.py` stops with a clear message, before anything is submitted, on duplicate IDs, IDs containing spaces or `/`, empty cells, non-amino-acid letters, or SMILES that RDKit can't parse.
- **File names:** the YAMLs are named `<proteinID>__<ligandID>.yaml`. Boltz names its outputs after the YAML.
- **Chain letters:** protein chains first (A, B, ...), then cofactors, then the screened ligand. `affinity: binder:` is the ligand's letter.
- **Affinity:** on by default, turned off with `--no_affinity`. Covalent jobs get it too.
- **Covalent experiments:** an experiment is wholly covalent or wholly not, because
  Boltz cannot mix affinity and non-affinity jobs in one folder (see below). A
  ligand names a `WARHEAD` and is written **unreacted**; a protein names where the
  ligand attaches in `COVALENT_SITE`, as `<chain>:<residue><position>` (`A:C16`).
  - `prepare_covalents.py` (step 3b) applies the warhead's reaction SMARTS, opens one
    valence on the attachment atom, names it, and writes `inputs/covalent_prepared.csv`.
    That file's presence is what makes the experiment covalent for `make_jobs.py`.
  - Warheads live in a `WARHEADS` dictionary at the top of that script, each with a
    reaction SMARTS (attachment atom mapped `:1` on both sides), a target residue and
    a target atom. `acrylamide` and `chloroacetamide` to begin with.
  - **More than one warhead match is an error**, even when the products are identical,
    because which end reacts is not something to guess at.
  - `COVALENT_SITE` is checked against the sequence - the residue there must really be
    the letter written down, and must be one the warheads react with. This catches an
    off-by-one or the wrong construct before any GPU time is spent.
  - Tested end to end on the laptop with `experiments/002_covalent_test/`, and the
    generated YAML for N-methylacrylamide matches the hand-written one that was
    verified on NEMO. Twelve bad-input cases were each confirmed to be caught.
  - **Run end to end on NEMO, 2026-11-05** (job 63198532, 3 pairs). All steps including
    3b ran under `module load boltz/2.2.0`; no errors; `results.csv` written.

- **The covalent bond forms at the right atom, but its length is not physical.**
  Measured from the output structures, the constrained atom was the ligand atom
  nearest the cysteine SG in every case, but the distances were:

  | pair | SG-atom | distance |
  | --- | --- | --- |
  | `testpep__ACRYL-ANILIDE` | SG-C12 | 1.27 A |
  | `testpep__NMA` | SG-C10 | 1.43 A |
  | `testpep__CHLOROACETAMIDE` | SG-C8 | 1.63 A |
  | the bonded verification job | SG-C10 | 1.67 A |

  A C-S single bond is about 1.81 A, so all four are compressed, one of them below a
  C-C single bond length. Boltz satisfies the constraint topologically while
  distorting the local geometry.

  **Do not rank or filter on bond length or bond geometry.** It says nothing about
  whether the pose is good. Use it only as a yes/no check that the constraint was
  applied at all, with a generous cutoff (under about 2.5 A).
- **MSA cache:** one folder per complex, `cache_msa/<hash of the SEQUENCE string>/`, holding one CSV per distinct chain, `<hash of the chain sequence>.csv`.
  - The MSAs are fetched by a protein-only Boltz run with `--use_msa_server`. Boltz pairs the MSAs for multi-chain complexes, and the pairing is valid only within that complex, which is why the cache is per complex.
  - `cache_msas.py` matches each fetched CSV to its chain by reading the first sequence row, which is the query sequence.
- **MSA paths:** the YAMLs use relative paths (`cache_msa/...`), and Boltz is always run from the project root. Confirmed working on NEMO, 2026-11-05.
- **Reruns:** the YAMLs are regenerated, cached MSAs are reused, and Boltz skips pairs that already have output. `--override` forces a full redo.
- **Intermediates:** `run_pipeline.sh` deletes `msa_jobs/` and `msa_outputs/` after step 3, so each run cleans up after itself. Everything worth keeping is in `cache_msa/` by then, and `make_msa_jobs.py` decides what to fetch from the cache alone, never from `msa_outputs/`. `set -e` means the deletion is skipped if step 3 fails, leaving the downloads to inspect. Without this, `msa_outputs/` accumulates across runs and step 3 re-scans and re-copies all of it every time.
- **Slurm:** one job runs every step in order. It keeps the settings from the old `run_boltz_job.sh` (ga100, 1 GPU, 32 CPUs, 100G, 7 days). Logs go to `experiments/X/logs/`.
- **Boltz settings:** v1 uses Boltz's defaults, set as variables at the top of `run_pipeline.sh`.
- **Python:** pandas, pyyaml and rdkit only. On the laptop they come from a conda env. On NEMO they come from the Python in `module load boltz/2.2.0`, which has all three (confirmed by the first full run, 2026-11-05).
- **Testing:** there is no pytest suite. A small experiment, `experiments/000_test_experiment/`, is run end to end instead. First full run on NEMO: 2026-11-05, all six steps, 4 pairs, `results.csv` written.
  - The GPU node can reach the ColabFold MSA server, so step 2 works inside the Slurm job.
  - `000_test_experiment` fetches ubiquitin's MSA twice, once alone and once paired in `two_chain_test`. The two files differ (16384 vs 16383 rows), which is why the cache is per complex and not per sequence.


## Boltz facts this relies on (v2.2.0 source)

- `boltz predict <folder>` runs every `.yaml` file in the folder. The folder must contain no other files or subfolders.
- `msa:` accepts `.a3m` or `.csv` (`key,sequence`). A relative path resolves from the folder Boltz is started in (confirmed on NEMO, 2026-11-05).
- The MSA server writes `<out>/boltz_results_<folder>/msa/<yamlname>_<entityid>.csv`. Paired rows only appear when the YAML has more than one distinct protein sequence.
- One YAML can't mix cached and server-fetched MSAs.
- The outputs go in `<out>/boltz_results_<folder>/predictions/<id>/`: `<id>_model_0.cif`, `confidence_<id>_model_0.json` and `affinity_<id>.json` (with `affinity_pred_value` = log10 IC50 in uM, and `affinity_probability_binary`).
- `msa: empty` is accepted and runs that chain in single sequence mode, with a warning.

### Covalent bonds (checked on NEMO `boltz/2.2.0`, 2026-11-05)

The Boltz **documentation** says `bond` constraints are only for CCD ligands. The
code does not enforce it: a constraint pointing at a **SMILES** ligand works. So a
covalent ligand needs no custom CCD entry and can be given as SMILES like any other.

- Constraint syntax, with a 1-based residue index. A SMILES ligand's residue index is `1`:

      constraints:
        - bond:
            atom1: [A, 16, SG]
            atom2: [B, 1, C10]

- **Atom naming for a SMILES ligand:** element symbol + RDKit canonical rank + 1,
  where the rank is computed on the molecule **after `AddHs`**. The hydrogens are
  then dropped from the model, so heavy-atom names start part-way up the numbering.
  Reproduce it with `Chem.CanonicalRankAtoms(Chem.AddHs(mol))`.
  For `[CH2]CC(=O)NC` this gives `C10 C11 C13 C14 N12 O9`, radical carbon `C10` -
  predicted with rdkit on the laptop and then confirmed in the run.
- **Evidence:** two 33-residue jobs with one cysteine, identical but for the
  constraint. Bonded: SG-C10 **1.67 A**, C10 the nearest ligand atom,
  confidence 0.89. No-constraint control: SG-C10 **14.15 A**, confidence 0.85.
  No warning or error about the SMILES ligand.
  (These jobs lived in a `verify_covalent/` folder that is no longer in the repo;
  the numbers above are the record of them. Worth rebuilding as the re-check to
  run after a Boltz upgrade.)
- **Affinity does run on a covalent complex.** A YAML can carry both a `bond`
  constraint and a `properties: affinity` block; Boltz writes the usual
  `affinity_<id>.json`. Confirmed on NEMO `boltz/2.2.0`, 2026-11-05.

  What the number means there is **not** established, and two things in that run
  argue for checking before ranking anything by it. Same ligand, same protein,
  bond constraint the only difference:

  | | covalent | non-covalent |
  | --- | --- | --- |
  | `affinity_pred_value` | 1.32 | 1.15 |
  | `affinity_probability_binary` | 0.39 | 0.10 |
  | the two models' `probability_binary` | 0.046 / 0.733 | 0.008 / 0.187 |

  - The imposed bond barely moved the potency but quadrupled the binder
    probability, so that score is partly reading the constraint, not the ligand.
  - The two ensemble members disagreed across the decision boundary on the
    covalent job, which makes the reported mean of 0.39 hard to trust.

  Caveat: one ligand, a designed sequence, `msa: empty`, 1 recycling step and
  1 diffusion sample. These are degraded settings and could explain the spread.
  Before a covalent screen is ranked on affinity, test it the way the sigma-2
  benchmark was tested: known covalent binders and non-binders, and check that
  `affinity_probability_binary` separates them.

### A jobs folder is all-affinity or no-affinity, never mixed

**This constrains the pipeline design, so it is not optional.** If *any* YAML in a
folder asks for affinity, the affinity stage runs over **every** record in that
folder and crashes on the ones that never asked, because they have no
`pre_affinity_<id>.npz`:

    FileNotFoundError: .../predictions/covalent_bonded/pre_affinity_covalent_bonded.npz

The cause is in `main.py`, `filter_inputs_affinity`. `r.affinity` is only used to
build the set of records that already *have* an affinity result; the returned list
does not filter on it:

    return Manifest([r for r in manifest.records if r.id not in existing])

Observed on NEMO `boltz/2.2.0`, 2026-11-05, and read from the v2.2.0 source.

Consequences:
- Covalent pairs (no affinity) and ordinary pairs (affinity) **cannot share one
  `jobs/` folder**. They need separate folders and separate `boltz predict` calls.
- The record list comes from `processed/records/`, not from the YAMLs present, so
  deleting a YAML is not enough to take a job out of a folder - its processed
  record has to go too.
