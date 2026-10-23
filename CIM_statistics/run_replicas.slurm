#!/bin/bash
# Runs the seed folders on line $SLURM_ARRAY_TASK_ID of the task file, one
# after another.  Submitted by ./submit_all.sh (which passes TASKS); do not
# sbatch it directly.  EDIT account, partition/qos, cores and wall time to
# match your main CIM jobs; the wall time must cover all seeds of a task.
#SBATCH -J cimrep
#SBATCH -A YOUR_ACCOUNT            # EDIT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=long
#SBATCH --qos=long
#SBATCH --time=2-00:00:00          # EDIT
#SBATCH -o logs/cimrep_%A_%a.out
#SBATCH -e logs/cimrep_%A_%a.err

module load lammps/29Sep2021_intel  # EDIT: your LAMMPS module
LMP=lmp_mpi                         # EDIT

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
