#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

: "${WORKROOT:=$HOME/build_dftbp_stack}"
: "${DFTBPLUS_VERSION:=25.1}"
: "${LAMMPS_VERSION:=22Jul2025}"
: "${OMP_NUM_THREADS:=1}"

DFTBPLUS_PREFIX="${DFTBPLUS_PREFIX:-$WORKROOT/dftbplus-${DFTBPLUS_VERSION}/_install_intelmpi}"
LAMMPS_EXE="${LAMMPS_EXE:-$WORKROOT/lammps-${LAMMPS_VERSION}/_build_dftbp_intelmpi/lmp}"
EXAMPLE_DIR="${EXAMPLE_DIR:-$REPO_ROOT/examples/md}"
INPUT_FILE="${INPUT_FILE:-md.in}"

if [[ ! -x "$LAMMPS_EXE" ]]; then
  echo "LAMMPS executable not found: $LAMMPS_EXE" >&2
  exit 1
fi

if [[ ! -d "$EXAMPLE_DIR" ]]; then
  echo "Example directory not found: $EXAMPLE_DIR" >&2
  exit 1
fi

export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS

cd "$EXAMPLE_DIR"
exec "$LAMMPS_EXE" -in "$INPUT_FILE"
