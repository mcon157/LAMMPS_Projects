#!/bin/bash
# =====================================================================
# setup_replicas.sh - one folder per seed inside each run folder.
#
#   CIM_replicas/
#     potentials/         WCH.tersoff WHe.table CHe.table HeHe.table w_eam4.fs
#     runs/WC_basal_C_50eV/   wche0.restart he_cycle_500.restart ...  (WC)
#     runs/W001_100eV/        w0.restart    he_cycle_500.restart ...  (W)
#     in_CIM_rep.step2  in_W_CIM_rep.step2  and the scripts
#
# Material from the cycle-0 restart, energy from the folder name (_<E>eV).
# For each save point present (0, 500 = he_cycle_500.restart, ...) it makes
# NSEED folders  seed_b<branch>_s<seed>/  (potentials, deck, restart,
# args.txt); LAMMPS gets rep = branch + seed.  Existing seed folders are
# left alone (re-run with a larger -n to add seeds).
#
#   ./setup_replicas.sh [options] [run folders]      default: all of runs/
#     -n seeds per save point (10, max 99)   -c cycles per seed (100)
#     -b save points ("0 500 1000")           -d dump every N cycles (50)
#     -x "extra LAMMPS args" (must match the main run)
#     -f rebuild existing seed folders (only before they have run)
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$(pwd)

NSEED=10; NCYC=100; BRANCHES="0 500 1000"; DUMP_EVERY=50; FORCE=0; EXTRA=""
while getopts "n:c:b:d:x:f" opt; do
  case $opt in
    n) NSEED=$OPTARG ;; c) NCYC=$OPTARG ;; b) BRANCHES=$OPTARG ;;
    d) DUMP_EVERY=$OPTARG ;; x) EXTRA=$OPTARG ;; f) FORCE=1 ;;
    *) sed -n '2,22p' "$0"; exit 1 ;;
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

# Seeds waiting to run, over all run folders (submit_all.sh finds the same).
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
