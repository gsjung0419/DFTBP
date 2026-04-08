# DFTBP Legacy README

This document keeps the original usage model for the older verified
stack. For current porting work on newer LAMMPS, see
[`README.md`](README.md).

## Supported Legacy Stack

Treat the following as the historical baseline:

- supported LAMMPS version: `29Oct2020`
- expected DFTB+ generation: legacy API compatible with the author's
  fork at <https://github.com/gsjung0419/dftbplus>
- reported upstream DFTB+ version in the original notes: `20.1`

The author notes that the code was tested with:

- Intel compiler 19 / 20
- MPICH 3.2.1 built with Intel compilers
- CMake 3.12 or 3.19
- LAMMPS `29Oct2020`

If you use newer LAMMPS or DFTB+ versions, assume that additional API
or linker porting work may be needed.

## Repository Layout

- `src/`: LAMMPS package source files
- `lib/dftbp/`: LAMMPS library makefiles and install helper
- `examples/md/`: graphene MD example input set

## Requirements and Limitations

- `fix dftbp` currently works only in serial
- only `units metal` is supported
- the simulation must be 3D
- the last argument of `fix dftbp` must currently be `NULL`
- the DFTB+ input file and Slater-Koster files must be available at run
  time
- set `OMP_NUM_THREADS` explicitly before running

The original package notes listed the following DFTB+ build options:

- `WITH_API=TRUE`
- `WITH_OMP=TRUE` or `FALSE`
- `WITH_DFTD3=TRUE` or `FALSE`
- `BUILD_SHARED_LIBS=TRUE`

## Building DFTB+

The recommended legacy path is to build the forked DFTB+ version linked
above, or to manually port the required API updates to a newer upstream
DFTB+ release.

Typical build steps:

1. Configure DFTB+ with API support enabled.
2. Build and install DFTB+.
3. Record the directory containing `libdftbplus`.

At minimum, make sure the installed DFTB+ build provides the API library
needed by this package.

## Installing the Package into LAMMPS

These instructions reflect the original manual installation flow for
LAMMPS `29Oct2020`.

1. Copy the library support files into the LAMMPS source tree:

   ```bash
   cp -R lib/dftbp /path/to/lammps/lib/
   ```

2. Copy this package source into the LAMMPS source tree:

   ```bash
   cp -R src /path/to/lammps/src/DFTBP
   ```

3. Edit [`lib/dftbp/Makefile.lammps.mpi`](lib/dftbp/Makefile.lammps.mpi)
   and point `dftbp_SYSLIB` to your DFTB+ installation, for example:

   ```make
   dftbp_SYSLIB = -L/path/to/dftbplus/lib -ldftbplus
   ```

4. Check that [`src/dftbplus.h`](src/dftbplus.h) matches the DFTB+ API
   header you built against.

5. Enable the package inside LAMMPS:

   ```bash
   cd /path/to/lammps/src
   make yes-dftbp
   make package-status
   ```

6. Build LAMMPS with your preferred target, for example:

   ```bash
   make mpi
   ```

The original README recommends building LAMMPS with the same compiler /
MPI toolchain used for DFTB+.

## `fix dftbp` Syntax

The implemented command form is:

```lammps
fix ID group-ID dftbp dftb_in.hsd NULL
```

Example from [`examples/md/md.in`](examples/md/md.in):

```lammps
fix             integ   all nve
fix             2       all dftbp dftb_in.hsd NULL
fix_modify      2       energy yes
fix             therm   all langevin 10.0 10.0 0.1 699483
```

Notes:

- the fourth argument must be the DFTB+ HSD input file
- the fifth argument must currently be `NULL`
- `fix_modify ... energy yes` is used so the DFTB+ contribution is
  included in the reported energy

## Example: Graphene MD

The example directory [`examples/md`](examples/md) contains:

- `md.in`: LAMMPS input script
- `graphene.data`: LAMMPS data file
- `graphene.gen`: DFTB+ geometry file
- `dftb_in.hsd`: DFTB+ input
- `C-C.skf`: Slater-Koster parameter file

Important:

- the atom order in `graphene.data` and `graphene.gen` must match
- review `dftb_in.hsd` and update it for your material system and SK
  files before production use

Suggested run command:

```bash
export OMP_NUM_THREADS=1
mpirun -np 1 lmp_mpi -in md.in
```

Use a single MPI rank unless you modify the code to support parallel
execution.

## Reported HPC Build Notes

The original notes included a user report from Iker Ortiz de Luzuriaga
at DIPC. The main points were:

- LAMMPS version used: `stable_29Oct2020`
- DFTB+ version used: `20.1`
- a forked DFTB+ build was also used
- on that system, one build provided `libmudpack.a` and the other
  provided a complete `libdftbplus.so`
- LAMMPS linker flags were adjusted manually to combine the required
  libraries
