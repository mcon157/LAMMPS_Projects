#!/bin/bash
# =====================================================================
# make_cim_replicas.sh - set up seed-varied CIM replicas that branch from
# the checkpoints of ONE existing CIM campaign (one surface, one energy).
#
# Each replica gets its own directory with the potentials, the replica
# deck, the right restart and a one-line LAMMPS argument list (args.txt):
#
#   <MAIN>/reps/b0000_r01/   cycles    1-100   from the cycle-0 restart
#   <MAIN>/reps/b0500_r01/   cycles  501-600   from he_cycle_500.restart
#   <MAIN>/reps/b1000_r01/   cycles 1001-1100  from he_cycle_1000.restart
#
# Usage:
#   ./make_cim_replicas.sh <WC|W> <MAIN_RUN_DIR> <E_eV> [NREP] [NCYC] [BRANCHES...]
#
#   ./make_cim_replicas.sh WC ~/runs/WC_basal_C_100eV 100              # 10 reps x 100 cycles at 0/500/1000
#   ./make_cim_replicas.sh WC ~/runs/WC_prism_50eV   50 10 100 0 500   # only branches 0 and 500
#   ./make_cim_replicas.sh W  ~/runs/W001_100eV      100 10 100 0      # pure W, branch 0 only
#
# Then submit what was created (the script prints this line too):
#   sbatch --array=1-<N> run_cim_replicas.slurm <MAIN>/reps/jobs.txt
#
# Optional environment overrides:
#   CYCLE0=<file>     name of the cycle-0 restart in MAIN
#                     (default WC: wche0.restart, W: w0.restart)
#   DUMP_EVERY=<n>    keep dump.bulk/full every n cycles (default 50)
#
# Nothing here touches the main run's own files.  Re-running the script
# (e.g. with a larger NREP) leaves existing replica directories alone and
# lists ONLY the new ones in reps/jobs.txt, so you can always submit the
# whole list.  FORCE=1 rebuilds existing directories too (this overwrites
# their restart - only do it before they have run).
# =====================================================================
set -euo pipefail

MAT=${1:?material: WC or W}
MAIN=${2:?main run directory}
E=${3:?He energy in eV}
NREP=${4:-10}            # keep < 100
NCYC=${5:-100}
shift $(( $# < 5 ? $# : 5 ))
BRANCHES=${@:-0 500 1000}
DUMP_EVERY=${DUMP_EVERY:-50}
HERE=$(cd "$(dirname "$0")" && pwd)

case "$MAT" in
  WC) DECK=in_CIM_rep.step2;   POTS="WCH.tersoff WHe.table CHe.table HeHe.table"
      CYCLE0=${CYCLE0:-wche0.restart}; EVAR=e;    NVAR=nimplants ;;
  W)  DECK=in_W_CIM_rep.step2; POTS="w_eam4.fs WHe.table HeHe.table"
      CYCLE0=${CYCLE0:-w0.restart};    EVAR=Eion; NVAR=nend ;;
  *)  echo "first argument must be WC or W" >&2; exit 1 ;;
esac
[ -f "$HERE/$DECK" ] || { echo "$DECK not found next to this script ($HERE)" >&2; exit 1; }
[ "$NREP" -lt 100 ] || { echo "NREP must be < 100 (replica numbers are branch + r)" >&2; exit 1; }
for P in $POTS; do
  [ -f "$MAIN/$P" ] || { echo "potential file $P not found in $MAIN" >&2; exit 1; }
done

mkdir -p "$MAIN/reps"
JOBS="$MAIN/reps/jobs.txt"
: > "$JOBS"

for B in $BRANCHES; do
  if [ "$B" -eq 0 ]; then SRC="$MAIN/$CYCLE0"; else SRC="$MAIN/he_cycle_${B}.restart"; fi
  if [ ! -f "$SRC" ]; then
    echo "skip branch $B: $SRC not found (the main run has not reached it)" >&2
    continue
  fi
  ISTART=$(( B + 1 ))
  IEND=$(( B + NCYC ))
  NEW=0
  for R in $(seq 1 "$NREP"); do
    D=$(printf "%s/reps/b%04d_r%02d" "$MAIN" "$B" "$R")
    if [ -f "$D/args.txt" ] && [ "${FORCE:-0}" != 1 ]; then
      continue                      # already set up (and maybe run): leave it alone
    fi
    NEW=$(( NEW + 1 ))
    mkdir -p "$D"
    for P in $POTS; do cp "$MAIN/$P" "$D/"; done
    cp "$HERE/$DECK" "$D/$DECK"
    if [ "$B" -eq 0 ]; then
      cp "$SRC" "$D/$CYCLE0"
      cp "$MAIN"/dump.bulk-0.lammps "$MAIN"/dump.full-0.lammps "$D/" 2>/dev/null || true
    else
      cp "$SRC" "$D/he_cycle.restart"
    fi
    # rep = branch + r is unique across the campaign (branches are multiples
    # of 500 and r < 100), so no two replicas share a random stream; the main
    # run is rep 0.
    REP=$(( B + R ))
    echo "-in $DECK -var $EVAR $E -var rep $REP -var istart $ISTART -var $NVAR $IEND -var dump_every $DUMP_EVERY" > "$D/args.txt"
    echo "$D" >> "$JOBS"
  done
  echo "branch $B: $NEW new replica(s) of $NREP, cycles $ISTART-$IEND (rep $((B+1))-$((B+NREP)))"
done
N=$(wc -l < "$JOBS")
echo "jobs listed in $JOBS ($N directories)"
[ "$N" -gt 0 ] && echo "submit with:  sbatch --array=1-$N run_cim_replicas.slurm $JOBS"
exit 0
