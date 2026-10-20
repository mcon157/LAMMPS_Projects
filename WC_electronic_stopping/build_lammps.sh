#!/bin/bash
# =====================================================================
# Build a current stable LAMMPS on an HPC cluster with the packages the
# W-C / He CIM decks need, plus EXTRA-FIX for fix electron/stopping.
#
#   MANYBODY  : pair tersoff/zbl (WC deck), pair eam/fs (pure-W deck)
#   EXTRA-FIX : fix electron/stopping, fix electron/stopping/fit
#   core      : pair table, hybrid/overlay, fix dt/reset, nve, nvt
#   OPENMP    : optional /omp accelerated styles (harmless to include)
#
# Run it on a LOGIN node (or better, an interactive compute node) from
# your home or project space:
#     bash build_lammps.sh
# Edit the MODULE block first - module names differ on every cluster.
# =====================================================================
set -euo pipefail

# ---- where things go -------------------------------------------------
SRC_ROOT=${SRC_ROOT:-$HOME/software}          # source checkout lives here
BRANCH=${BRANCH:-stable}                      # 'stable' = current yearly stable + bugfix updates
NPROC=${NPROC:-8}                             # parallel compile jobs (be polite on login nodes)

# ---- MODULES: EDIT FOR YOUR CLUSTER ---------------------------------
# Check what exists with:  module avail gcc ; module avail openmpi ; module avail cmake
# Needs: CMake >= 3.20, a C++17 compiler (a recent GCC/Intel/oneAPI), and an MPI.
# Use the SAME compiler+MPI pair at build time and in your job scripts.
module purge
module load gcc            # e.g. gcc/12.x
module load openmpi        # e.g. openmpi/4.1.x  (or intel-mpi / mpich)
module load cmake          # must report >= 3.20
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

# the version string you will quote in the paper / methods
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
