#!/bin/bash
# =====================================================================
# setup_replicas.sh - create one folder per seed inside each run folder.
#
# Expected layout (run everything from the CIM_replicas folder):
#
#   CIM_replicas/
#     potentials/            WCH.tersoff WHe.table CHe.table HeHe.table w_eam4.fs
#     runs/
#       WC_basal_C_50eV/     wche0.restart  he_cycle_500.restart  he_cycle_1000.restart
#       W001_100eV/          w0.restart     he_cycle_500.restart ...
#     in_CIM_rep.step2  in_W_CIM_rep.step2
#     setup_replicas.sh  run_replicas.slurm  summarize.slurm  submit_all.sh
#     summarize_replicas.py
#
# The material comes from the cycle-0 restart in the run folder
# (wche0.restart = WC, w0.restart = W), the energy from the folder name
# (..._<E>eV).  For every save point present (0 = the cycle-0 restart,
# 500 = he_cycle_500.restart, ...) it makes NSEED folders:
#
#   runs/WC_basal_C_50eV/seed_b0000_s01/   cycles    1-100
#   runs/WC_basal_C_50eV/seed_b0500_s01/   cycles  501-600
#   runs/WC_basal_C_50eV/seed_b1000_s01/   cycles 1001-1100
#
# each holding the potentials, the input file, its restart and args.txt.
# The seed number entering LAMMPS is rep = branch + s, unique in a run
# folder (the main run is rep 0).
#
# Usage:
#   ./setup_replicas.sh                       # every folder in runs/
#   ./setup_replicas.sh runs/WC_basal_C_50eV  # only these folders
#   ./setup_replicas.sh -n 20 -c 100 -b "0 500" runs/WC_basal_W_100eV
#
# Options:  -n seeds per save point (default 10, max 99)
#           -c cycles per seed       (default 100)
#           -b save points to use    (default "0 500 1000"; missing ones are skipped)
#           -d dump every N cycles   (default 50)
#           -x "extra LAMMPS args" appended to every args.txt
#              (e.g. -x "-var Tsub 739"; must match the main run)
#           -f rebuild seed folders that already exist (overwrites their
#              restart - only before they have run)
#
# Existing seed folders are never touched (unless -f), so re-running with
# a larger -n only adds seeds.  It ends by listing, for every run folder,
# the seeds that have not been submitted or started yet; submit_all.sh
# submits exactly those.
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(pwd)

NSEED=10; NCYC=100; BRANCHES="0 500 1000"; DUMP_EVERY=50; FORCE=0; EXTRA=""
while getopts "n:c:b:d:x:f" opt; do
  case $opt in
    n) NSEED=$OPTARG ;; c) NCYC=$OPTARG ;; b) BRANCHES=$OPTARG ;;
    d) DUMP_EVERY=$OPTARG ;; x) EXTRA=$OPTARG ;; f) FORCE=1 ;;
    *) sed -n '2,45p' "$0"; exit 1 ;;
  esac
done
shift $((OPTIND - 1))
[ "$NSEED" -ge 1 ] && [ "$NSEED" -le 99 ] || { echo "-n must be 1-99" >&2; exit 1; }

if [ $# -gt 0 ]; then RUNS=("$@"); else RUNS=(runs/*/); fi
[ -d potentials ] || { echo "no potentials/ folder in $ROOT" >&2; exit 1; }

TOTAL=0
for RUN in "${RUNS[@]}"; do
  RUN=${RUN%/}
  [ -d "$RUN" ] || { echo "skip $RUN: not a folder" >&2; continue; }
  NAME=$(basename "$RUN")
  if [ -f "$RUN/wche0.restart" ]; then
    MAT=WC; CYCLE0=wche0.restart; DECK=in_CIM_rep.step2;   EVAR=e;    NVAR=nimplants
    POTS="WCH.tersoff WHe.table CHe.table HeHe.table"
  elif [ -f "$RUN/w0.restart" ]; then
    MAT=W;  CYCLE0=w0.restart;    DECK=in_W_CIM_rep.step2; EVAR=Eion; NVAR=nend
    POTS="w_eam4.fs WHe.table HeHe.table"
  else
    echo "skip $NAME: no wche0.restart or w0.restart" >&2; continue
  fi
  if [[ "$NAME" =~ _([0-9]+)eV$ ]]; then E=${BASH_REMATCH[1]}; else
    echo "skip $NAME: folder name must end in _<E>eV" >&2; continue
  fi
  for P in $POTS; do
    [ -f "potentials/$P" ] || { echo "missing potentials/$P (needed for $NAME)" >&2; exit 1; }
  done

  for B in $BRANCHES; do
    if [ "$B" -eq 0 ]; then SRC="$RUN/$CYCLE0"; else SRC="$RUN/he_cycle_${B}.restart"; fi
    if [ ! -f "$SRC" ]; then
      [ "$B" -ne 0 ] && echo "  $NAME: no he_cycle_${B}.restart - save point $B skipped"
      continue
    fi
    ISTART=$(( B + 1 )); IEND=$(( B + NCYC )); NEW=0
    for S in $(seq 1 "$NSEED"); do
      D=$(printf "%s/seed_b%04d_s%02d" "$RUN" "$B" "$S")
      if [ -f "$D/args.txt" ] && [ "$FORCE" != 1 ]; then continue; fi
      mkdir -p "$D"
      for P in $POTS; do cp "potentials/$P" "$D/"; done
      cp "$DECK" "$D/"
      if [ "$B" -eq 0 ]; then cp "$SRC" "$D/$CYCLE0"; else cp "$SRC" "$D/he_cycle.restart"; fi
      REP=$(( B + S ))
      echo "-in $DECK -var $EVAR $E -var rep $REP -var istart $ISTART -var $NVAR $IEND -var dump_every $DUMP_EVERY $EXTRA" > "$D/args.txt"
      NEW=$(( NEW + 1 ))
    done
    echo "  $NAME ($MAT, $E eV): save point $B -> $NEW new seed folder(s), cycles $ISTART-$IEND"
    TOTAL=$(( TOTAL + NEW ))
  done
done

# Report what is waiting to run, over ALL run folders (not only the ones
# named in this call).  submit_all.sh finds the same seeds by itself.
echo "$TOTAL seed folder(s) created now.  Waiting to run (never submitted, never started):"
NWAIT=0
for RUN in runs/*/; do
  n=0
  for D in "${RUN}"seed_b*_s*/; do
    [ -f "$D/args.txt" ] && [ ! -f "$D/log.lammps" ] && [ ! -f "$D/submitted" ] && n=$((n + 1))
  done
  [ "$n" -gt 0 ] && echo "  ${RUN%/}: $n"
  NWAIT=$((NWAIT + n))
done
echo "  total: $NWAIT"
[ "$NWAIT" -gt 0 ] && echo "submit with ./submit_all.sh  (see ./submit_all.sh -h for packing several seeds per job)"
exit 0
