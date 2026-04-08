# DFTBP

LAMMPS package for running DFTB+-based simulations through a custom
`fix dftbp` command.

This top-level README is for the current porting effort against newer
LAMMPS source trees. For the original older stack and manual install
notes, see [`README.legacy.md`](README.legacy.md).

## Document Split

- [`README.md`](README.md): current porting status, staging script, and
  known blockers on newer LAMMPS
- [`README.legacy.md`](README.legacy.md): legacy `29Oct2020`-era usage
  notes
- [`PORTING_LOG.md`](PORTING_LOG.md): chronological debug notes
- [`MAINTENANCE.md`](MAINTENANCE.md): branch and version policy

Version selection happens at `git clone` / `git checkout` time, not
during the LAMMPS install step.

## Current Scope

This repository adds a LAMMPS fix that wraps the DFTB+ API so that
molecular dynamics and relaxation can be performed with DFTB+ forces
inside LAMMPS.

The current goal on `main` is:

- make the package stage cleanly into newer LAMMPS source trees
- support both legacy make builds and current CMake builds
- keep porting work separate from the legacy-only usage notes

Reference:
Jung, Gang Seob, Stephan Irle, and Bobby G. Sumpter.
"Dynamic aspects of graphene deformation and fracture from approximate
density functional theory." *Carbon* 190 (2022): 183-193.
DOI: [10.1016/j.carbon.2022.01.002](https://doi.org/10.1016/j.carbon.2022.01.002)
Publisher page:
[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0008622322000021)

Related projects:
- Official DFTB+: <https://github.com/dftbplus/dftbplus>
- Recommended legacy fork:
  <https://github.com/gsjung0419/dftbplus>

## Current Test Status

Latest local smoke-test target:

- LAMMPS source tree: `22Jul2025`
- DFTB+ install used for tests: `25.1` with a local MM API patch

Status as of `2026-04-07`:

- `src/fix_dftb.cpp` was updated for the current LAMMPS `Fix` API
- standalone compilation against `lammps-22Jul2025/src` succeeds
- legacy make-based LAMMPS build compiles the DFTBP package sources
- CMake-based LAMMPS build now links successfully against a patched
  shared DFTB+ `25.1` install
- the current local executable build reports `DFTBP OPENMP` as installed
  packages
- the current shared-library build also produces a working Python-facing
  `liblammps_mpi.dylib`

The main known macOS gap is no longer DFTB+ linking for the shared CMake
path. The remaining rough edge is the default `install-python` target in
offline environments, since it tries to refresh Python wheel-build
dependencies from the network.

## DFTB+ 25.1 Patch Helper

If you start from a fresh DFTB+ `25.1` source tree, use the helper under
[`dftbplus_25patch`](dftbplus_25patch/README.md) before building and
installing DFTB+.

That helper adds the extra MM API entry points used by this repository:

- `dftbp_get_energyparts`
- `dftbp_get_virial`

## One-Shot Staging Script

Use [`tools/integrate_lammps_dftbp.py`](tools/integrate_lammps_dftbp.py)
to copy the package into a LAMMPS source tree and optionally build it.

What it does:

- copies `src/` into `LAMMPS/src/DFTBP`
- copies `lib/dftbp/` into `LAMMPS/lib/dftbp`
- installs the package CMake module into `LAMMPS/cmake/Modules/Packages`
- updates the LAMMPS make and CMake package lists
- auto-writes `lib/dftbp/Makefile.lammps*` from a DFTB+ install prefix
- can run either a make-based or CMake-based build

Stage only:

```bash
python3 tools/integrate_lammps_dftbp.py \
  --lammps-dir /path/to/lammps \
  --dftbplus-prefix /path/to/dftbplus-install
```

Stage and run legacy make build:

```bash
python3 tools/integrate_lammps_dftbp.py \
  --lammps-dir /path/to/lammps \
  --dftbplus-prefix /path/to/dftbplus-install \
  --build --build-system make --make-target mpi
```

Stage and run CMake build:

```bash
python3 tools/integrate_lammps_dftbp.py \
  --lammps-dir /path/to/lammps \
  --dftbplus-prefix /path/to/dftbplus-install \
  --build --build-system cmake \
  --build-dir /path/to/lammps/_build_dftbp
```

If your LAMMPS build needs extra compiler or MPI settings, pass them
through with repeated `--cmake-arg`.

## Shared Library and Python Usage

For Python use, there are still two pieces involved:

- the LAMMPS shared library itself
- the Python `lammps` package that loads that shared library

A Python wheel is only one packaging method for that second piece. It is
not a different API from the more direct Linux-style workflow where you
set `PYTHONPATH` and point the dynamic loader at a shared
`liblammps_*.so`.

Current local results on macOS:

- shared LAMMPS build succeeded with
  `liblammps_mpi.dylib`
- Python import and `lammps(name='mpi')` runtime construction succeeded
- `has_package('DFTBP')` returned `True`
- the default CMake `install-python` target did not complete in this
  offline environment because its helper tries to upgrade/install `pip`,
  `build`, `wheel`, and `setuptools`
- an offline local wheel build and install did succeed once the shared
  library was already present

So, in practice:

- on Linux, using the shared library directly remains a normal workflow
- on macOS, the same idea works; the wheel step is just a convenient way
  to place the Python files and bundled shared library together

## Requirements and Current Limitations

Before building or running, keep the following constraints in mind.

- `fix dftbp` currently works only in serial
- only `units metal` is supported
- the simulation must be 3D
- the recommended `fix dftbp` form is now `dftb_in.hsd off NULL`
- the optional implicit-GB slot accepts `off` or `on`
- the last argument must still currently be `NULL`
- the DFTB+ input file and Slater-Koster files must be available at run
  time
- set `OMP_NUM_THREADS` explicitly before running

For DFTB+, the original package notes listed these important options:

- `WITH_API=TRUE`
- `WITH_OMP=TRUE` or `FALSE`
- `WITH_DFTD3=TRUE` or `FALSE`
- `BUILD_SHARED_LIBS=TRUE`

## `fix dftbp` Syntax

The implemented command form is:

```lammps
fix ID group-ID dftbp dftb_in.hsd implicit-model coulomb-compute
```

Recommended current form:

```lammps
fix ID group-ID dftbp dftb_in.hsd off NULL
```

Optional implicit-GB form:

```lammps
fix ID group-ID dftbp dftb_in.hsd on NULL
```

Example from [`examples/md/md.in`](examples/md/md.in):

```lammps
fix             integ   all nve
fix             2       all dftbp dftb_in.hsd off NULL
fix_modify      2       energy yes
fix             therm   all langevin 10.0 10.0 0.1 699483
```

Notes:

- the fourth argument must be the DFTB+ HSD input file
- the fifth argument toggles implicit GB mode: `off` or `on`
- the sixth argument is reserved for a LAMMPS Coulomb compute ID and
  must currently be `NULL`
- with `on`, DFTBP uses the coordinate-only update path and reports a
  zero virial
- the legacy `NULL` and `GB` tokens are still accepted for backward
  compatibility and map to `off` and `on`
- the legacy one-extra-argument form is also still accepted when the
  final token is just `NULL`
- `fix_modify ... energy yes` includes the DFTB+ contribution in the
  reported energy

## Next Debug Target

The next step on newer LAMMPS is to make the Python-side installation
flow more robust across online and offline environments, while keeping
the shared-library path straightforward for Linux and macOS users.

## Credits

Environment setup and installation simplification were completed with
OpenAI Codex.

## Contact

Bug reports and comments can be sent to:

- `jungg@ornl.gov`
- `gs4phone@gmail.com`
