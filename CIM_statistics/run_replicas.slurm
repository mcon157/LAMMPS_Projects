#!/bin/bash
# Runs the seed folders listed on line $SLURM_ARRAY_TASK_ID of the task
# file, one after another.  Submitted by submit_all.sh, which writes the
# task file (logs/tasks_<date>.txt) and passes it as TASKS.  Don't sbatch
# this file directly - use ./submit_all.sh (-t for a one-seed test).
#
# EDIT the four lines marked EDIT (same as your main CIM jobs).  If a job
# runs several seeds, the wall time must cover all of them (or pass
# -T HH:MM:SS to submit_all.sh).
#SBATCH -J cimrep
#SBATCH -A ACF-UTK0016
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=2-00:00:00 
#SBATCH -o logs/cimrep_%A_%a.out
#SBATCH -e logs/cimrep_%A_%a.err

module load lammps/29Sep2021_intel
LMP=lmp_mpi

: "${TASKS:?no task file - submit with ./submit_all.sh}"
TID=${SLURM_ARRAY_TASK_ID:-1}
LIST=$(sed -n "${TID}p" "$TASKS")
[ -n "$LIST" ] || { echo "line $TID of $TASKS is empty" >&2; exit 1; }

FAILED=0
for D in $LIST; do
  echo "=== $(date '+%F %T')  seed folder: $D"
  if [ -f "$D/log.lammps" ]; then
    echo "    already started before (log.lammps exists) - skipped"; continue
  fi
  cd "$D" || { echo "    cannot cd - skipped"; FAILED=$((FAILED + 1)); continue; }
  echo "    lammps args: $(cat args.txt)"
  if srun -n "${SLURM_NTASKS:-1}" $LMP $(cat args.txt) -log log.lammps; then
    echo "    finished $(date '+%F %T')"
  else
    echo "    LAMMPS FAILED (see $D/log.lammps)"; FAILED=$((FAILED + 1))
  fi
done
echo "=== done: $FAILED seed(s) failed in this job"
exit "$FAILED"
