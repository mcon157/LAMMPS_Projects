#!/bin/bash
# =====================================================================
# Build current stable LAMMPS with what the decks need, plus EXTRA-FIX
# for fix electron/stopping:
#   MANYBODY   tersoff/zbl (WC), eam/fs (W)      EXTRA-FIX  electron/stopping
#   core       table, hybrid/overlay, dt/reset   OPENMP     optional /omp styles
# Run on a login or interactive node:  bash build_lammps.sh
# Edit the module block first (names differ on every cluster).
# =====================================================================
set -euo pipefail

# ---- where things go -------------------------------------------------
SRC_ROOT=${SRC_ROOT:-$HOME/software}          # source checkout lives here
BRANCH=${BRANCH:-stable}                      # 'stable' = current stable release + bugfixes
NPROC=${NPROC:-8}                             # parallel compile jobs

# ---- modules: EDIT FOR YOUR CLUSTER (this set worked on one cluster) ----
# Needs CMake >= 3.20, a C++17 compiler and an MPI; use the same pair in your job scripts.
# Check:  module avail gcc ; module avail openmpi ; module avail cmake
module purge
module load gcc/10.2.0
module load openmpi/4.0.5-gcc
module load cmake/3.23.2-intel
module load hdf5/1.10.7-gcc
module load boost/1.78.0-gcc
module list

cmake --version | head -1
mpicxx --version | head -1

# ---- get the source --------------------------------------------------
mkdir -p "$SRC_ROOT"
cd "$SRC_ROOT"
if [ ! -d lammps ]; then
    git clone -b "$BRANCH" --depth 1 https://github.com/lammps/lammps.git lammps
else
    echo "lammps/ already exists - updating it"
    git -C lammps fetch --depth 1 origin "$BRANCH"
    git -C lammps checkout -B "$BRANCH" FETCH_HEAD
fi

# version string to quote in methods sections
grep LAMMPS_VERSION lammps/src/version.h

# ---- configure -------------------------------------------------------
cd lammps
BUILD_DIR=build-elstop
rm -rf "$BUILD_DIR"

cmake -S cmake -B "$BUILD_DIR" \
      -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_INSTALL_PREFIX="$SRC_ROOT/lammps-install" \
      -D BUILD_MPI=on \
      -D BUILD_OMP=on \
      -D PKG_MANYBODY=on \
      -D PKG_EXTRA-FIX=on \
      -D PKG_OPENMP=on \
      -D LAMMPS_MACHINE=elstop

# ---- compile + install ----------------------------------------------
cmake --build "$BUILD_DIR" --parallel "$NPROC"
cmake --install "$BUILD_DIR"

LMP="$SRC_ROOT/lammps-install/bin/lmp_elstop"
echo
echo "Executable: $LMP"

# ---- verify the styles you need are really in the binary -------------
echo "---- installed packages ----"
"$LMP" -h | sed -n '/Installed packages:/,/^$/p'
for s in electron/stopping tersoff/zbl dt/reset table eam/fs; do
    if "$LMP" -h | grep -qw "$s"; then
        echo "OK      $s"
    else
        echo "MISSING $s   <-- check the PKG_ flags above"
    fi
done

echo
echo "Add this to your job scripts (after the same module loads):"
echo "    LMP=$LMP"
echo "    srun \$LMP -in in_CIM_ES.step2 -var e 100"
