#!/bin/bash
# Writes replica_yields.csv / replica_summary.csv / replica_yields.png into
# every run folder, and runs/all_replica_summary.csv.  submit_all.sh queues
# this to start when the seed array has finished; you can also run it any
# time by hand:  sbatch summarize.slurm   (or just: python3 summarize_replicas.py)
#SBATCH -J cimrep_sum
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 00:30:00
#SBATCH -o logs/cimrep_summary_%j.out

module load python                  # EDIT: any python3 with numpy (matplotlib optional)
cd "$SLURM_SUBMIT_DIR"
python3 summarize_replicas.py runs
