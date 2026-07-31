#!/bin/bash
# Run one case on an Idun compute node.
#
#   sbatch run_idun.slurm                 # case1, both phases
#   sbatch run_idun.slurm case2           # case2, both phases
#   sbatch run_idun.slurm all             # case1, then case2, unattended
#   sbatch run_idun.slurm case1 bo        # resume the optimisation phase only
#
# With `all`, the job runs case1 to the end, moves results/ to
# results_archive/case1_<timestamp>/, and starts case2 in a fresh results/. Size
# --time for both cases in that mode.
#
# The job runs run_supervised.py, which starts MATLAB itself, watches it, and
# relaunches it for as long as the run needs. Neither half runs on a login node,
# and MATLAB runs headless, so no X11 connection exists that a dropped ssh
# session can take down.
#
# MATLAB has no command window in this mode. Its output goes to
# results/logs/matlab_console.log. To watch it while the job runs:
#
#   ssh <the node in squeue>
#   cd <this directory> && python watch_matlab_log.py
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
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm_mfbo_%j.out
#SBATCH --mail-type=FAIL,TIME_LIMIT
##SBATCH --mail-user=your.name@ntnu.no

set -uo pipefail

# Settings
CASE="${1:-case1}"
PHASE="${2:-both}"
[ "$CASE" = "all" ] && CASE="case1 case2"
MATLAB_MODULE="MATLAB/2025a"
CONDA_MODULE="Anaconda3/2024.02-1"   # `module spider anaconda` lists what is installed
CONDA_ENV="mfbo-nmpc"                # the name declared in environment.yml
MATLAB_WORKERS=31                    # NumWorkers in main_initialization.m and main_BO.m

# Argument and resource checks
# These run before MATLAB starts, so a mistake costs a few seconds of queue time
# rather than a day of a reservation that cannot produce a usable result.
case "$PHASE" in
    both|init|bo) ;;
    *) echo "[wrapper] phase must be both, init or bo, got '$PHASE'" >&2
       exit 2 ;;
esac

if [ "$PHASE" = "bo" ] && [ ! -f results/surrogate/phi_coeffs.mat ]; then
    echo "[wrapper] the bo phase needs results/surrogate/phi_coeffs.mat, which the" >&2
    echo "[wrapper] design phase publishes. Run the init phase first." >&2
    exit 2
fi

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

# The run
# run_supervised.py owns the MATLAB process: it waits for the serve loop before
# the first request is published, relaunches MATLAB whenever it exits or stops
# holding matlab.lock with a request pending, and stops it on the way out. Its
# SIGTERM handler runs that teardown when Slurm ends the step, so `python -u`
# keeps the wrapper log in order up to that point.
echo "[wrapper] MATLAB console: results/logs/matlab_console.log"
# CASE may name more than one case, so the split is intended here.
# shellcheck disable=SC2086
python -u run_supervised.py --case $CASE --phase "$PHASE"
STATUS=$?

if [ "$STATUS" -ne 0 ]; then
    echo "[wrapper] the run exited with status $STATUS. Last lines of the MATLAB console:" >&2
    tail -n 40 results/logs/matlab_console.log >&2
else
    echo "[wrapper] the run finished."
fi

if [ -f results/logs/restarts.csv ]; then
    echo "[wrapper] MATLAB restarts recorded this run and earlier ones:"
    tail -n 20 results/logs/restarts.csv
fi

# Peak resident memory of this job step, for sizing --mem on the next run.
sstat -j "${SLURM_JOB_ID:-0}.batch" --format=JobID,MaxRSS,MaxVMSize 2>/dev/null

exit "$STATUS"
