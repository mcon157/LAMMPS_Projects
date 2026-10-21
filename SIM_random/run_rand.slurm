#!/bin/bash
# One random surface per array task (task ID = surface ID s): build it
# (step 1), then NIMP He impacts at each energy (step 2).
#   sbatch --array=1-50 run_rand.slurm WC                     # 50 random WC surfaces
#   sbatch --array=1-50 run_rand.slurm W                      # 50 random W surfaces
#   sbatch --array=1-30 run_rand.slurm WC bas "-var beta 0 -var gamma 0"         # basal refs
#   sbatch --array=1-50 run_rand.slurm W r "-var T 300 -var R 35 -var a 3.152"   # Jussila
# The 3rd argument goes to both steps.  Run in a folder with the decks and
# potentials (W: w_eam4.fs WHe.table HeHe.table; WC: WCH.tersoff WHe.table
# CHe.table HeHe.table).  EDIT the marked lines for your cluster.
#SBATCH -J randsim
#SBATCH -N 1
#SBATCH --ntasks=8                  # EDIT
#SBATCH -t 48:00:00                 # EDIT: time a 20-impact pilot first
#SBATCH -o rand_%A_%a.out

module load lammps                  # EDIT: same build as the other SIM runs
LMP=lmp_mpi                         # EDIT

MAT=${1:?W or WC}
TAG=${2:-r}
EXTRA=${3:-}
ENERGIES="50 100 200"               # EDIT: W validation adds 300 500 (see README)
NIMP=1000                           # impacts per surface per energy
S=$SLURM_ARRAY_TASK_ID

srun -n "$SLURM_NTASKS" $LMP -in in_rand_${MAT}.step1 -var s $S -var tag $TAG $EXTRA -log log.${MAT}_${TAG}${S}.step1
for E in $ENERGIES; do
  srun -n "$SLURM_NTASKS" $LMP -in in_rand_${MAT}.step2 -var s $S -var tag $TAG -var e $E -var nimp $NIMP $EXTRA \
       -log log.${MAT}_${TAG}${S}_${E}eV.step2
done
