#!/bin/bash
# One CIM replica per array task.  The job list comes from make_cim_replicas.sh.
#   sbatch --array=1-<N> run_cim_replicas.slurm <MAIN>/reps/jobs.txt
# (make_cim_replicas.sh prints the exact line, with N filled in.)
# Works for WC and W replicas alike: each directory's args.txt names its deck.
#
# EDIT the four lines marked EDIT for your cluster (same module / binary /
# core count / wall time as your main CIM jobs).
#SBATCH -J cimrep
#SBATCH -N 1
#SBATCH --ntasks=12                 # EDIT
#SBATCH -t 24:00:00                 # EDIT: ~100 cycles; time one replica first
#SBATCH -o cimrep_%A_%a.out

module load lammps                  # EDIT: the same LAMMPS build as the main runs
LMP=lmp_mpi                         # EDIT

JOBS=${1:?usage: sbatch --array=1-N run_cim_replicas.slurm <MAIN>/reps/jobs.txt}
if ! [[ "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-9]+$ ]]; then
  echo "submit with --array=1-N (N = number of lines in $JOBS)" >&2
  exit 1
fi
D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBS")
[ -n "$D" ] && [ -f "$D/args.txt" ] || { echo "no replica directory on line $SLURM_ARRAY_TASK_ID of $JOBS" >&2; exit 1; }
cd "$D"
echo "replica directory: $D"
echo "lammps args: $(cat args.txt)"
srun -n "${SLURM_NTASKS:-1}" $LMP $(cat args.txt) -log log.lammps
