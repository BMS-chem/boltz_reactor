#!/bin/env bash
#
# submit_boltz.sh - run this on the cluster LOGIN NODE.
#
# It does no computing itself. It checks the experiment folder exists, then hands
# the real work to run_pipeline.sh on a GPU node via sbatch.
#
# Usage:
#   ./Scripts/submit_boltz.sh <experiment>
#   ./Scripts/submit_boltz.sh <experiment> --no_affinity
#   ./Scripts/submit_boltz.sh <experiment> --override --time 1-00:00:00
#
# Options:
#   --no_affinity   skip the affinity prediction (structure only)
#   --override      redo predictions that already have output (Boltz skips them by default)
#   --time          Slurm time limit (default below)

set -e            # stop as soon as any command fails
set -u            # stop if an undefined variable is used


# ---------------------------------------------------------------------------
# 1. Settings for SBATCH setup
# ---------------------------------------------------------------------------
# These are site-specific: set them for your own cluster.
TIME_LIMIT="7-00:00:00"       # days-hours:minutes:seconds
PARTITION="ga100"             # a GPU partition
CPUS=32
MEMORY="100G"
EMAIL=""                      # leave empty for no email notification


# ---------------------------------------------------------------------------
# 2. Read the command-line arguments
# ---------------------------------------------------------------------------
# The first argument is the experiment name. Anything after it is an option.
if [ $# -lt 1 ]; then
    echo "ERROR: please give an experiment name, e.g. ./Scripts/submit_boltz.sh boltz_run_001"
    exit 1
fi

EXPERIMENT="$1"
shift                 # drop the experiment name, leaving only the options

NO_AFFINITY="no"      # these are passed on to run_pipeline.sh
OVERRIDE="no"

while [ $# -gt 0 ]; do
    case "$1" in
        --no_affinity) NO_AFFINITY="yes"; shift ;;
        --override)    OVERRIDE="yes"; shift ;;
        --time)        TIME_LIMIT="$2"; shift 2 ;;
        *) echo "ERROR: unknown option '$1'"; exit 1 ;;
    esac
done


# ---------------------------------------------------------------------------
# 3. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
# Doing it this way means you can run the script from anywhere.
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
PROJECT_DIR="$( dirname "$SCRIPT_DIR" )"
EXPERIMENT_DIR="$PROJECT_DIR/experiments/$EXPERIMENT"

if [ ! -d "$EXPERIMENT_DIR" ]; then
    echo "ERROR: experiment folder not found: $EXPERIMENT_DIR"
    exit 1
fi
if [ ! -f "$EXPERIMENT_DIR/inputs/proteins.csv" ] || [ ! -f "$EXPERIMENT_DIR/inputs/ligands.csv" ]; then
    echo "ERROR: $EXPERIMENT_DIR/inputs/ needs both proteins.csv and ligands.csv"
    exit 1
fi

# The logs folder has to exist before sbatch starts, or Slurm has nowhere to write.
mkdir -p "$EXPERIMENT_DIR/logs"


# ---------------------------------------------------------------------------
# 4. Submit the job
# ---------------------------------------------------------------------------
# The log paths are given here, on the sbatch command line, rather than in the
# #SBATCH lines inside run_pipeline.sh. That is the only way to get them into each
# experiment's own folder, because #SBATCH lines cannot contain variables.
# %j is replaced by the Slurm job number, so reruns do not overwrite old logs.
#
# --export passes our settings into the job as environment variables.
# ALL keeps the rest of the environment as well.

# Only ask for email if an address was set above. An address has no spaces, so
# leaving $MAIL_OPTIONS unquoted below splits it into the two options correctly.
MAIL_OPTIONS=""
if [ -n "$EMAIL" ]; then
    MAIL_OPTIONS="--mail-type=END,FAIL --mail-user=$EMAIL"
fi

sbatch \
    --job-name="boltz_$EXPERIMENT" \
    --output="$EXPERIMENT_DIR/logs/boltz_%j.out" \
    --error="$EXPERIMENT_DIR/logs/boltz_%j.err" \
    --time="$TIME_LIMIT" \
    --partition="$PARTITION" \
    --cpus-per-task="$CPUS" \
    --mem="$MEMORY" \
    --gres=gpu:1 \
    $MAIL_OPTIONS \
    --export=ALL,EXPERIMENT="$EXPERIMENT",PROJECT_DIR="$PROJECT_DIR",NO_AFFINITY="$NO_AFFINITY",OVERRIDE="$OVERRIDE" \
    "$SCRIPT_DIR/run_pipeline.sh"

echo ""
echo "Submitted experiment '$EXPERIMENT'"
echo "  logs:    $EXPERIMENT_DIR/logs/"
echo "  watch:   squeue -u $USER"
echo "  follow:  tail -f $EXPERIMENT_DIR/logs/boltz_<jobid>.out"
