#!/bin/env bash
#
# run_pipeline.sh - this runs on a GPU NODE. Do not run it by hand on the login node.
#
# Submit it with:  ./Scripts/submit_boltz.sh <experiment>
#
# It runs the whole pipeline in order, in one Slurm job:
#   1. work out which proteins still need an MSA
#   2. fetch those MSAs from the ColabFold server (a protein-only Boltz run)
#   3. copy the fetched MSAs into the shared cache
#   3b. prepare the covalent ligands, if this experiment has any
#   4. make one job YAML per protein x ligand pair
#   5. run the Boltz predictions
#   6. collect the results into one table
#
# The #SBATCH settings (time, memory, partition, log paths) are all given by
# submit_boltz.sh on the sbatch command line, so there are no #SBATCH lines here.

set -e            # stop as soon as any command fails
set -u            # stop if an undefined variable is used

# ---------------------------------------------------------------------------
# 1. Boltz settings (currently Boltz' defaults)
# ---------------------------------------------------------------------------
DIFFUSION_SAMPLES=1       # poses generated per job
RECYCLING_STEPS=3
SAMPLING_STEPS=200
OUTPUT_FORMAT=mmcif       # mmcif or pdb

# ---------------------------------------------------------------------------
# 2. Set up
# ---------------------------------------------------------------------------
# EXPERIMENT, PROJECT_DIR, NO_AFFINITY and OVERRIDE are passed in by submit_boltz.sh.
EXPERIMENT_DIR="$PROJECT_DIR/experiments/$EXPERIMENT"

# Boltz is always run from the project folder, because the job YAMLs give their
# MSA paths relative to it (cache_msa/...).
cd "$PROJECT_DIR"

mkdir -p "$EXPERIMENT_DIR/msa_jobs" "$EXPERIMENT_DIR/msa_outputs" \
         "$EXPERIMENT_DIR/jobs" "$EXPERIMENT_DIR/outputs" cache_msa

echo "=========================================="
echo "Experiment:  $EXPERIMENT"
echo "Project:     $PROJECT_DIR"
echo "Affinity:    $( [ "$NO_AFFINITY" = "yes" ] && echo OFF || echo ON )"
echo "Override:    $OVERRIDE"
echo "Started:     $(date)"
echo "=========================================="

module purge
module load boltz/2.2.0

# Turn yes/no settings into the flags the commands actually take.
AFFINITY_FLAG=""
if [ "$NO_AFFINITY" = "yes" ]; then AFFINITY_FLAG="--no_affinity"; fi

OVERRIDE_FLAG=""
if [ "$OVERRIDE" = "yes" ]; then OVERRIDE_FLAG="--override"; fi


# ---------------------------------------------------------------------------
# 3. Step 1: find the proteins that have no cached MSA
# ---------------------------------------------------------------------------
# Writes one protein-only YAML into msa_jobs/ for each complex not already in
# cache_msa/, and empties msa_jobs/ first. If every protein is cached it writes
# nothing, and the fetch below is skipped.
echo ""
echo "--- Step 1: checking the MSA cache ---"
python "$PROJECT_DIR/Scripts/make_msa_jobs.py" --experiment "$EXPERIMENT"


# ---------------------------------------------------------------------------
# 4. Step 2: fetch the missing MSAs
# ---------------------------------------------------------------------------
# A protein-only Boltz run with --use_msa_server. The predicted structures are
# thrown away; we only want the MSA files it downloads on the way.
echo ""
echo "--- Step 2: fetching MSAs ---"

# Count the YAML files that step 1 wrote. 2>/dev/null keeps it quiet if there are none.
N_MSA_JOBS=$(ls -1 "$EXPERIMENT_DIR"/msa_jobs/*.yaml 2>/dev/null | wc -l)

if [ "$N_MSA_JOBS" -eq 0 ]; then
    echo "Every protein already has a cached MSA, nothing to fetch."
else
    echo "Fetching MSAs for $N_MSA_JOBS protein(s)..."
    boltz predict "$EXPERIMENT_DIR/msa_jobs" \
        --out_dir "$EXPERIMENT_DIR/msa_outputs" \
        --cache "$PROJECT_DIR/cache_weights" \
        --use_msa_server \
        --output_format "$OUTPUT_FORMAT"
fi


# ---------------------------------------------------------------------------
# 5. Step 3: copy the fetched MSAs into the shared cache
# ---------------------------------------------------------------------------
# Copies the CSVs Boltz downloaded into cache_msa/<complex_hash>/<chain_hash>.csv,
# so later experiments with the same protein do not fetch them again.
echo ""
echo "--- Step 3: adding MSAs to the cache ---"
python "$PROJECT_DIR/Scripts/cache_msas.py" --experiment "$EXPERIMENT"


# ---------------------------------------------------------------------------
# 5b. Tidy up the MSA working folders
# ---------------------------------------------------------------------------
# Everything worth keeping is now in cache_msa/, so the fetch intermediates go:
# the raw ColabFold downloads, the throwaway protein-only structures, and the
# job YAMLs (which are rewritten from the CSVs on every run anyway).
#
# This runs only if step 3 succeeded, because set -e stops the script on failure.
# That is deliberate: if the caching fails, the downloads stay put to look at.
#
# The 'case' guards the rm: $EXPERIMENT_DIR must start with $PROJECT_DIR/experiments/
# and must not end there. rm -rf on an empty or wrong variable would be expensive.
case "$EXPERIMENT_DIR" in
    "$PROJECT_DIR/experiments/"?*)
        rm -rf "$EXPERIMENT_DIR/msa_outputs" "$EXPERIMENT_DIR/msa_jobs"
        echo "Removed msa_outputs/ and msa_jobs/ (the MSAs are in cache_msa/)"
        ;;
    *)
        echo "WARNING: '$EXPERIMENT_DIR' is not where it should be, so nothing was deleted."
        ;;
esac


# ---------------------------------------------------------------------------
# 5c. Step 3b: prepare the covalent ligands
# ---------------------------------------------------------------------------
# If any ligand names a WARHEAD, this works out its bound (open-valence) form and
# which atom attaches to the protein, and writes inputs/covalent_prepared.csv.
# If none do, it says so and does nothing, so this runs for every experiment.
echo ""
echo "--- Step 3b: preparing covalent ligands ---"
python "$PROJECT_DIR/Scripts/prepare_covalents.py" --experiment "$EXPERIMENT"


# ---------------------------------------------------------------------------
# 6. Step 4: make the job files
# ---------------------------------------------------------------------------
# Checks the input CSVs and writes one YAML per protein x ligand pair into jobs/.
# It stops the whole pipeline here if anything in the CSVs is wrong.
echo ""
echo "--- Step 4: making the job files ---"
python "$PROJECT_DIR/Scripts/make_jobs.py" --experiment "$EXPERIMENT" $AFFINITY_FLAG


# ---------------------------------------------------------------------------
# 7. Step 5: run the predictions
# ---------------------------------------------------------------------------
# Boltz runs every YAML in jobs/ in one go, and skips any pair that already has
# output, unless --override is given.
echo ""
echo "--- Step 5: running the predictions ---"
boltz predict "$EXPERIMENT_DIR/jobs" \
    --out_dir "$EXPERIMENT_DIR/outputs" \
    --cache "$PROJECT_DIR/cache_weights" \
    --diffusion_samples "$DIFFUSION_SAMPLES" \
    --recycling_steps "$RECYCLING_STEPS" \
    --sampling_steps "$SAMPLING_STEPS" \
    --output_format "$OUTPUT_FORMAT" \
    $OVERRIDE_FLAG


# ---------------------------------------------------------------------------
# 8. Step 6: collect the results
# ---------------------------------------------------------------------------
# Reads the affinity and confidence JSON files Boltz wrote for each pair and puts
# them into one table: experiments/<experiment>/results.csv
echo ""
echo "--- Step 6: collecting the results ---"
python "$PROJECT_DIR/Scripts/collect_results.py" --experiment "$EXPERIMENT"

echo ""
echo "=========================================="
echo "Finished:    $(date)"
echo "Results:     $EXPERIMENT_DIR/results.csv"
echo "Structures:  $EXPERIMENT_DIR/outputs/"
echo "=========================================="
