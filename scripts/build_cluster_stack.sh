#!/usr/bin/env bash

set -euo pipefail

# Cluster build helper for:
# - patched DFTB+ 25.1
# - LAMMPS 22Jul2025
# - DFTBP package from this repository

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

: "${WORKROOT:=$HOME/build_dftbp_stack}"
: "${DFTBPLUS_VERSION:=25.1}"
: "${LAMMPS_VERSION:=22Jul2025}"
: "${BUILD_JOBS:=8}"

DFTB_ARCHIVE="dftbplus-${DFTBPLUS_VERSION}.tar.xz"
LAMMPS_ARCHIVE="lammps-${LAMMPS_VERSION}.tar.gz"
DFTB_SRC="$WORKROOT/dftbplus-${DFTBPLUS_VERSION}"
LAMMPS_SRC="$WORKROOT/lammps-${LAMMPS_VERSION}"
DFTB_BUILD="$DFTB_SRC/_build_intelmpi"
DFTB_INSTALL="$DFTB_SRC/_install_intelmpi"
LAMMPS_BUILD="$LAMMPS_SRC/_build_dftbp_intelmpi"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing required command: $1" >&2
    exit 1
  }
}

print_step() {
  echo
  echo "==> $1"
}

print_step "Checking required commands"
require_cmd curl
require_cmd tar
require_cmd python3
require_cmd cmake
require_cmd make
require_cmd mpiicx
require_cmd mpiicpx
require_cmd mpiifx

print_step "Checking Intel MPI / MKL environment"
: "${MKLROOT:?MKLROOT is not set. Load your Intel MKL module first.}"
: "${I_MPI_ROOT:?I_MPI_ROOT is not set. Load your Intel MPI module first.}"

mkdir -p "$WORKROOT"
cd "$WORKROOT"

print_step "Downloading source archives if needed"
if [[ ! -f "$DFTB_ARCHIVE" ]]; then
  curl -L -O "https://github.com/dftbplus/dftbplus/releases/download/${DFTBPLUS_VERSION}/${DFTB_ARCHIVE}"
fi
if [[ ! -f "$LAMMPS_ARCHIVE" ]]; then
  curl -L -O "https://download.lammps.org/tars/${LAMMPS_ARCHIVE}"
fi

print_step "Extracting source archives if needed"
if [[ ! -d "$DFTB_SRC" ]]; then
  tar -xf "$DFTB_ARCHIVE"
fi
if [[ ! -d "$LAMMPS_SRC" ]]; then
  tar -xf "$LAMMPS_ARCHIVE"
fi

print_step "Patching DFTB+ source"
python3 "$REPO_ROOT/dftbplus_25patch/apply_patch.py" \
  --dftbplus-src "$DFTB_SRC"

print_step "Configuring DFTB+"
cmake -S "$DFTB_SRC" \
  -B "$DFTB_BUILD" \
  -DCMAKE_INSTALL_PREFIX="$DFTB_INSTALL" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=mpiicx \
  -DCMAKE_CXX_COMPILER=mpiicpx \
  -DCMAKE_Fortran_COMPILER=mpiifx \
  -DWITH_MPI=TRUE \
  -DWITH_API=TRUE \
  -DBUILD_SHARED_LIBS=TRUE \
  -DENABLE_DYNAMIC_LOADING=FALSE \
  -DWITH_PYTHON=FALSE

print_step "Building and installing DFTB+"
cmake --build "$DFTB_BUILD" --parallel "$BUILD_JOBS"
cmake --install "$DFTB_BUILD"

export DFTBPLUS_PREFIX="$DFTB_INSTALL"
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

print_step "Staging DFTBP into LAMMPS and building"
python3 "$REPO_ROOT/tools/integrate_lammps_dftbp.py" \
  --lammps-dir "$LAMMPS_SRC" \
  --dftbplus-prefix "$DFTBPLUS_PREFIX" \
  --build \
  --build-system cmake \
  --build-dir "$LAMMPS_BUILD" \
  --jobs "$BUILD_JOBS" \
  --cmake-arg=-D --cmake-arg=CMAKE_C_COMPILER=mpiicx \
  --cmake-arg=-D --cmake-arg=CMAKE_CXX_COMPILER=mpiicpx \
  --cmake-arg=-D --cmake-arg=BUILD_MPI=on

print_step "Validating final executable"
test -x "$LAMMPS_BUILD/lmp"
"$LAMMPS_BUILD/lmp" -h | grep -q 'DFTBP'

cat <<EOF

Build completed successfully.

DFTB+ install:
  $DFTB_INSTALL

LAMMPS executable:
  $LAMMPS_BUILD/lmp

Recommended runtime environment:
  export DFTBPLUS_PREFIX="$DFTB_INSTALL"
  export LD_LIBRARY_PATH="\$DFTBPLUS_PREFIX/lib:\$LD_LIBRARY_PATH"
  export OMP_NUM_THREADS=1

EOF
