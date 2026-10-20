#!/bin/bash
# =====================================================================
# submit_all.sh - submit every seed folder that is waiting to run.
#
# It scans runs/*/seed_b*_s*/ (or only the run folders you name) for
# seeds that were never submitted (no 'submitted' file) and never started
# (no log.lammps), splits them into at most -j jobs, and submits ONE array
# job.  Each array task runs its seeds one after the other.  A summary job
# follows when the whole array has ended.
#
#   ./submit_all.sh                    one array task per seed (if your queue allows it)
#   ./submit_all.sh -j 1               ALL waiting seeds in a single job, one after another
#   ./submit_all.sh -j 5               5 jobs, seeds shared out evenly
#   ./submit_all.sh -j 5 -r 2          5 jobs, at most 2 running at the same time
#   ./submit_all.sh -t                 test: only the first waiting seed, no summary
#   ./submit_all.sh -j 3 runs/WC_basal_C_100eV runs/WC_prism_100eV
#
# Options
#   -j N        number of jobs (array tasks) to use; default = one per seed
#   -r N        at most N of them running at once (sbatch --array=...%N)
#   -T HH:MM:SS wall time per job.  A job running K seeds needs K x (time of
#               one seed); without -T the #SBATCH -t of run_replicas.slurm is used.
#   -n          no summary job (saves one submission; run
#               python3 summarize_replicas.py yourself later)
#   -t          test: submit only the first waiting seed, without a summary
#   -d          dry run: print what would be submitted, submit nothing
#
# Your queue's submission limit counts every array task as a job, so use
# -j to stay under it (the summary job counts as one more unless -n).
# Check the limit with:  sacctmgr show assoc user=$USER format=user,maxsubmit
#                        sacctmgr show qos format=name,maxsubmitpu
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(pwd)

NJ=0; MAXRUN=""; WALL=""; SUMMARY=1; TEST=0; DRY=0
while getopts "j:r:T:ntdh" opt; do
  case $opt in
    j) NJ=$OPTARG ;; r) MAXRUN=$OPTARG ;; T) WALL=$OPTARG ;;
    n) SUMMARY=0 ;; t) TEST=1; SUMMARY=0 ;; d) DRY=1 ;;
    *) sed -n '2,32p' "$0"; exit 0 ;;
  esac
done
shift $((OPTIND - 1))
if [ $# -gt 0 ]; then RUNS=("$@"); else RUNS=(runs/*/); fi

# ---- seeds waiting to run ------------------------------------------------
SEEDS=()
for RUN in "${RUNS[@]}"; do
  for D in "${RUN%/}"/seed_b*_s*/; do
    D=${D%/}
    [ -f "$D/args.txt" ] && [ ! -f "$D/log.lammps" ] && [ ! -f "$D/submitted" ] && SEEDS+=("$ROOT/$D")
  done
done
NS=${#SEEDS[@]}
[ "$NS" -gt 0 ] || { echo "no seed folders waiting to run (run ./setup_replicas.sh first?)"; exit 0; }
[ "$TEST" = 1 ] && { SEEDS=("${SEEDS[0]}"); NS=1; NJ=1; }
[ "$NJ" -le 0 ] || [ "$NJ" -gt "$NS" ] && NJ=$NS

# ---- share the seeds out over NJ tasks (round robin) ---------------------
mkdir -p logs
STAMP=$(date +%Y%m%d_%H%M%S)
TASKS="$ROOT/logs/tasks_$STAMP.txt"           # line i = seed folders of array task i
declare -a LINE
for ((k = 0; k < NS; k++)); do
  t=$(( k % NJ ))
  LINE[$t]="${LINE[$t]:-} ${SEEDS[$k]}"
done
: > "$TASKS"
for ((t = 0; t < NJ; t++)); do echo "${LINE[$t]# }" >> "$TASKS"; done
PER=$(( (NS + NJ - 1) / NJ ))
echo "$NS seed(s) waiting -> $NJ job(s), up to $PER seed(s) each (run one after another)"
for RUN in "${RUNS[@]}"; do
  n=$(grep -o "${RUN%/}/seed_" "$TASKS" | wc -l || true)
  [ "$n" -gt 0 ] && echo "  ${RUN%/}: $n"
done

ARR="1-$NJ"; [ -n "$MAXRUN" ] && ARR="$ARR%$MAXRUN"
OPTS=(--parsable --array="$ARR" --export=ALL,TASKS="$TASKS")
[ -n "$WALL" ] && OPTS+=(--time="$WALL")
if [ "$PER" -gt 1 ] && [ -z "$WALL" ]; then
  echo "NOTE: each job runs up to $PER seeds back to back - make sure the wall time"
  echo "      (#SBATCH -t in run_replicas.slurm, or -T here) covers $PER seeds."
fi
if [ "$DRY" = 1 ]; then
  echo "dry run: sbatch ${OPTS[*]} run_replicas.slurm"
  [ "$SUMMARY" = 1 ] && echo "dry run: sbatch --dependency=afterany:<array id> summarize.slurm"
  echo "task list written to $TASKS (nothing submitted, nothing marked)"
  exit 0
fi

JID=$(sbatch "${OPTS[@]}" run_replicas.slurm)
for D in "${SEEDS[@]}"; do echo "$JID" > "$D/submitted"; done   # never submitted twice
echo "submitted array job $JID ($NJ task(s)); task list $TASKS"
if [ "$SUMMARY" = 1 ]; then
  SID=$(sbatch --parsable --dependency=afterany:"$JID" summarize.slurm)
  echo "summary job $SID will run when $JID has ended"
fi
