# Porting Log

## Goal

Make this DFTBP package reusable against newer LAMMPS source trees with:

- legacy make-based package installation
- current CMake-based LAMMPS builds

## Current Test Target

- LAMMPS source tree: `lammps-22Jul2025`
- DFTB+ install used for later link tests:
  `/Users/8gj/Sources/dftbplus-25.1/_install`

## Findings So Far

### 2026-04-07

Code-level compatibility was tested first before dealing with external
library link details.

Observed compile errors against latest LAMMPS headers:

1. `virial_flag` is no longer available in the current `Fix` API.
2. `THERMO_ENERGY` is no longer a valid `setmask()` flag.

Applied source changes:

- replaced the old `virial_flag` usage with current fix flags:
  `energy_global_flag`, `virial_global_flag`, `thermo_energy`,
  `thermo_virial`
- removed `THERMO_ENERGY` from `setmask()`

Result:

- `src/fix_dftb.cpp` now compiles successfully as a standalone object
  against `lammps-22Jul2025/src`
- `src/dftbplus.h` header guard was cleaned up while doing the port

Manual make-based integration test:

- staged `src/` and `lib/dftbp/` into `lammps-22Jul2025`
- enabled the package with `make yes-dftbp`
- full `make mpi` build compiled the DFTBP package sources
- final link failed on unresolved DFTB+ / Fortran runtime symbols

Observed conclusion from the make flow:

- the main LAMMPS-side API compile errors appear fixed
- the remaining make-based blocker is external library linking

CMake integration work:

- added `cmake/Modules/Packages/DFTBP.cmake`
- added `tools/integrate_lammps_dftbp.py` to automate staging and build
  setup
- first CMake attempt failed because old make-generated files were still
  present in the source tree
- the staging script now runs `make no-all purge` before a CMake build
- a later CMake build configured successfully and compiled through the
  DFTBP sources
- final link still failed due to unresolved DFTB+ / Fortran runtime
  symbols on the current local test platform

Observed conclusion from the CMake flow:

- package registration inside current LAMMPS CMake now works
- compilation reaches the final executable link step
- the remaining blocker is DFTB+ link-line completeness, not package
  source registration

Documentation split:

- `README.md` is now intended for current porting and staging work
- `README.legacy.md` keeps the old `29Oct2020`-style instructions

DFTB+ companion patch helper:

- added `dftbplus_25patch/README.md`
- added `dftbplus_25patch/apply_patch.py`
- helper targets a fresh DFTB+ `25.1` source tree
- helper applies the MM API additions required by `fix dftbp`
- local checks covered script parsing, dry-run against `dftbplus-25.1`,
  and idempotent re-run behavior on a temporary patched tree

## Next Steps

1. Make DFTB+ link flags robust for static and shared installs
2. Record a known-good current build recipe once linking succeeds
3. Keep the staging script usable for freshly downloaded LAMMPS trees
4. Add a minimal smoke test after the link step is working
