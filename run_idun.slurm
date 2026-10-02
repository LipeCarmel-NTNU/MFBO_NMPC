#!/bin/bash
# Run one case on an Idun compute node.
#
#   sbatch run_idun.slurm baseline        # SF arm, both phases
#   sbatch run_idun.slurm case2 bo        # MF arm, optimisation phase only
#   sbatch run_idun.slurm                 # case1, both phases
#
# The run writes to results/<arm>_<case>_s<sobol_seed>, which run_config.py
# derives from the case and the Sobol seed of this checkout, or to
# MFBO_RESULTS_DIR when it is exported before sbatch. One job runs one case. Two
# campaigns started from one checkout would share inbox/theta.txt and
# matlab.lock, so each campaign that runs at the same time needs its own clone.
#
# The job runs run_supervised.py, which starts MATLAB itself, watches it, and
# relaunches it for as long as the run needs. Neither half runs on a login node,
# and MATLAB runs headless, so no X11 connection exists that a dropped ssh
# session can take down.
#
# MATLAB has no command window in this mode. Its output goes to
# <results folder>/logs/matlab_console.log. To watch it while the job runs:
#
#   ssh <the node in squeue>
#   cd <this directory> && python watch_matlab_log.py --case <case>
#
# The account below is the one association this user holds, at normal priority.
# `sacctmgr show assoc format=Account%15,User,QOS | grep $USER` lists them, and
# `idun-slurm-quota` shows the CPU hours left this month.
#
#SBATCH --job-name=mfbo-nmpc
#SBATCH --account=nv-ikp-ps
#SBATCH --partition=CPUQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm_mfbo_%j.out
#SBATCH --mail-type=FAIL,TIME_LIMIT
##SBATCH --mail-user=your.name@ntnu.no

set -uo pipefail

# Settings
CASE="${1:-case1}"
PHASE="${2:-both}"
MATLAB_MODULE="MATLAB/2025a"
CONDA_MODULE="Anaconda3/2024.02-1"   # `module spider anaconda` lists what is installed
CONDA_ENV="mfbo-nmpc"                # the name declared in environment.yml
MATLAB_WORKERS=8                     # NumWorkers in main_initialization.m and main_BO.m

# Argument and resource checks
# These run before MATLAB starts, so a mistake costs a few seconds of queue time
# rather than a day of a reservation that cannot produce a usable result.
case "$PHASE" in
    both|init|bo) ;;
    *) echo "[wrapper] phase must be both, init or bo, got '$PHASE'" >&2
       exit 2 ;;
esac

CORES="${SLURM_CPUS_PER_TASK:-0}"
if [ "$CORES" -le "$MATLAB_WORKERS" ]; then
    echo "[wrapper] the reservation has $CORES core(s) and the MATLAB scripts open" >&2
    echo "[wrapper] $MATLAB_WORKERS workers plus a client. Raise --cpus-per-task to at" >&2
    echo "[wrapper] least $((MATLAB_WORKERS + 1)), or lower NumWorkers in the .m files." >&2
    exit 2
fi

cd "${SLURM_SUBMIT_DIR:-$(pwd)}" || exit 2

echo "[wrapper] job ${SLURM_JOB_ID:-none} on $(hostname), $CORES cores, ${SLURM_MEM_PER_NODE:-?} MB"
echo "[wrapper] case $CASE, phase $PHASE"

# Job-private MATLAB state
# parpool('Processes') keeps its job files under the MATLAB preference
# directory, which defaults to a fixed path in $HOME. Two jobs of this kind
# running at once would share it and interfere. A per-job prefdir and TMPDIR
# keep them apart, and both live on the node rather than on the shared
# filesystem.
JOB_TMP="${TMPDIR:-/tmp}/mfbo_${SLURM_JOB_ID:-$$}"
mkdir -p "$JOB_TMP/prefs" "$JOB_TMP/tmp" || exit 2
export MATLAB_PREFDIR="$JOB_TMP/prefs"
export TMPDIR="$JOB_TMP/tmp"
trap 'rm -rf "$JOB_TMP"' EXIT

# Environment
module purge
module load "$MATLAB_MODULE" || { echo "[wrapper] could not load $MATLAB_MODULE" >&2; exit 2; }

if module load "$CONDA_MODULE" 2>/dev/null; then
    CONDA_BASE="$(conda info --base 2>/dev/null)"
    # shellcheck disable=SC1091
    [ -n "$CONDA_BASE" ] && source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV" || {
        echo "[wrapper] could not activate the conda environment '$CONDA_ENV'." >&2
        echo "[wrapper] Create it once with: conda env create -f environment.yml" >&2
        exit 2
    }
else
    echo "[wrapper] $CONDA_MODULE is not available. Falling back to the python on PATH."
fi

python -c "import torch, botorch, gpytorch, numpy, scipy" || {
    echo "[wrapper] the driver's packages are not importable in this environment." >&2
    exit 2
}

command -v matlab >/dev/null || {
    echo "[wrapper] matlab is not on PATH after loading $MATLAB_MODULE." >&2
    exit 2
}

# Results folder
# The same derivation the driver makes, so the checks and the log paths below
# look at the tree the run writes. Exported, so the driver and MATLAB use it.
if [ -z "${MFBO_RESULTS_DIR:-}" ]; then
    MFBO_RESULTS_DIR="$(python -c 'import sys; from run_config import results_dir_for; print(results_dir_for(sys.argv[1]))' "$CASE")" || {
        echo "[wrapper] could not derive the results folder for case '$CASE'." >&2
        exit 2
    }
fi
export MFBO_RESULTS_DIR
echo "[wrapper] results folder $MFBO_RESULTS_DIR"

# The optimisation phase reads the design rows, which an earlier init phase, or
# a copy from the SF clone, puts there. It fits or publishes phi itself.
if [ "$PHASE" = "bo" ] && [ ! -f "$MFBO_RESULTS_DIR/init/results.csv" ]; then
    echo "[wrapper] the bo phase needs $MFBO_RESULTS_DIR/init/results.csv, which" >&2
    echo "[wrapper] the init phase writes. Run it first, or copy the design in." >&2
    exit 2
fi

# The run
# run_supervised.py owns the MATLAB process: it waits for the serve loop before
# the first request is published, relaunches MATLAB whenever it exits or stops
# holding matlab.lock with a request pending, and stops it on the way out. Its
# SIGTERM handler runs that teardown when Slurm ends the step, so `python -u`
# keeps the wrapper log in order up to that point.
LOGS="$MFBO_RESULTS_DIR/logs"
echo "[wrapper] MATLAB console: $LOGS/matlab_console.log"
python -u run_supervised.py --case "$CASE" --phase "$PHASE"
STATUS=$?

if [ "$STATUS" -ne 0 ]; then
    echo "[wrapper] the run exited with status $STATUS. Last lines of the MATLAB console:" >&2
    tail -n 40 "$LOGS/matlab_console.log" >&2
else
    echo "[wrapper] the run finished."
fi

if [ -f "$LOGS/restarts.csv" ]; then
    echo "[wrapper] MATLAB restarts recorded this run and earlier ones:"
    tail -n 20 "$LOGS/restarts.csv"
fi

# Peak resident memory of this job step, for sizing --mem on the next run.
sstat -j "${SLURM_JOB_ID:-0}.batch" --format=JobID,MaxRSS,MaxVMSize 2>/dev/null

exit "$STATUS"
