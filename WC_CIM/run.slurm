#!/bin/bash
#SBATCH -J wccim
#SBATCH -A YOUR_ACCOUNT            # EDIT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH -e wche.e%j
#SBATCH -o wche.o%j

# ============================================================
# Self-chaining job for a CIM campaign on a short (3 h) partition.
# Each job reads last_cycle.txt, queues its successor (afterany, so the
# chain survives the wall-time kill), then runs LAMMPS from the next cycle
# to NIMP, reading he_cycle.restart.  Seeds depend on the cycle number, so
# resuming continues the random stream.
#   sbatch run.slurm          (or: sbatch run.slurm 19  = at most 19 more links)
# Fresh campaign: first
#   rm -f last_cycle.txt he_cycle.restart sputtering.dat sputter_summary.dat cycle_time.dat
# (the deck starts these files only at cycle 1).
# ============================================================

# ---------------- settings ----------------
NIMP=1000                          # total implant cycles wanted (-> -var nimplants)
INPUT=in.step2                   # LAMMPS input file (in_CIM_092926.step2 copied as in.step2)
LAST_FILE=last_cycle.txt         # resume anchor written by the deck
STATE_RESTART=he_cycle.restart   # carry-over state (cycles > 1)
NTASKS=24                        # MPI ranks; = --ntasks-per-node
MAX_LINKS=${1:-19}               # this job + this many more (safety cap)
# ------------------------------------------------

cd "$SLURM_SUBMIT_DIR" || exit 1

# 1) next cycle = last completed + 1
if [[ -f "$LAST_FILE" ]]; then
    LAST=$(awk 'NF{v=$1+0} END{print v+0}' "$LAST_FILE")
else
    LAST=0
fi
ISTART=$(( LAST + 1 ))

echo "[chain] job=$SLURM_JOB_ID  last_cycle=$LAST  istart=$ISTART  nimplants=$NIMP  max_links=$MAX_LINKS"

# 2) all done: stop the chain
if (( ISTART > NIMP )); then
    echo "[chain] all $NIMP cycles complete -- chain finished."
    exit 0
fi

# 3) resuming needs the restart; stop without queueing a successor if it is missing
if (( ISTART > 1 )) && [[ ! -f "$STATE_RESTART" ]]; then
    echo "[chain] ERROR: $LAST_FILE says cycle $LAST done but $STATE_RESTART is missing."
    echo "[chain] Not resuming. Restore the restart, or start fresh (see header)."
    exit 1
fi

# 4) queue the next link before running (afterany: runs however this job ends)
if (( MAX_LINKS > 0 )); then
    sbatch --dependency=afterany:$SLURM_JOB_ID run.slurm $(( MAX_LINKS - 1 ))
fi

# 5) run (with -in, so the jumps resolve)
module load lammps                 # EDIT: your LAMMPS module (29Oct2020 or newer)
LMP=lmp_mpi                        # EDIT
srun -n $NTASKS $LMP -var istart $ISTART -var nimplants $NIMP -in $INPUT
