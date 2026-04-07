# DFTBP

LAMMPS package for running DFTB+-based simulations through a custom
`fix dftbp` command.

This repository adds a LAMMPS fix that wraps the DFTB+ API so that
molecular dynamics and relaxation can be performed with DFTB+ forces
inside LAMMPS.

## Supported Legacy Stack

This repository should currently be treated as a legacy integration.

- Supported LAMMPS version: `29Oct2020`
- Expected DFTB+ generation: legacy API compatible with the author's
  fork at <https://github.com/gsjung0419/dftbplus>
- Reported upstream DFTB+ version in the original notes: `20.1`

If you are using newer LAMMPS or DFTB+ versions, assume that code or API
porting work will be required. Those combinations should be considered
unsupported unless they are tested and documented separately.

Version selection happens when cloning or checking out the repository,
not during the LAMMPS install step. If `main` is maintained as the
legacy-stable branch, a normal clone is the correct install path.

For repository maintenance guidance, see
[`MAINTENANCE.md`](MAINTENANCE.md).

Reference:
Gang Seob Jung, Stephan Irle, and Bobby Sumpter,
"Dynamic aspects of graphene deformation and fracture from approximate
density functional theory," *Carbon* (2022).

Related projects:
- Official DFTB+: <https://github.com/dftbplus/dftbplus>
- Recommended fork for this package:
  <https://github.com/gsjung0419/dftbplus>

## Status

This package was developed by GS Jung at ORNL based on
`fix_latte.cpp` / `fix_latte.h`.

The author notes that the code was tested with:
- Intel compiler 19 / 20
- MPICH 3.2.1 built with Intel compilers
- CMake 3.12 or 3.19
- LAMMPS 29 Oct 2020

Newer versions of LAMMPS or DFTB+ may require API updates. The author
specifically notes that the recommended DFTB+ fork contains API changes
needed by this package and that those changes were not merged upstream
at the time of writing. In practice, this repository should be read as
"works with the legacy stack above unless proven otherwise."

## Repository Layout

- `src/`: LAMMPS package source files (`fix_dftb.cpp`, `fix_dftb.h`,
  `dftbplus.h`)
- `lib/dftbp/`: LAMMPS library makefiles and install helper
- `examples/md/`: graphene MD example input set

## Requirements and Current Limitations

Before building or running, keep the following constraints in mind.

- `fix dftbp` currently works only in serial. The code aborts if
  LAMMPS is launched with more than one MPI rank.
- Only `units metal` is supported.
- The simulation must be 3D.
- The last argument of `fix dftbp` must currently be `NULL`.
  Passing a LAMMPS Coulomb potential is not supported.
- The DFTB+ input file and Slater-Koster files must be available at run
  time.
- Set `OMP_NUM_THREADS` explicitly before running. Otherwise DFTB+ may
  try to use all available cores.

For DFTB+, the README originally listed the following tested build
options:

- `WITH_API=TRUE` (required)
- `WITH_OMP=TRUE` or `FALSE`
- `WITH_DFTD3=TRUE` or `FALSE`
- `BUILD_SHARED_LIBS=TRUE`

`BUILD_SHARED_LIBS=TRUE` is especially important if you need a shared
library build (for example when using Python-based LAMMPS workflows).

## Building DFTB+

The recommended path is to build the forked DFTB+ version linked above,
or to manually port the required API updates to a newer upstream DFTB+
release.

Typical build steps:

1. Configure DFTB+ with API support enabled.
2. Build and install DFTB+.
3. Record the directory containing `libdftbplus`.

At minimum, make sure the installed DFTB+ build provides the API library
needed by this package.

## Installing the Package into LAMMPS

The original instructions were written for LAMMPS `29Oct2020`. The same
layout may work for nearby versions, but that has not been verified in
this repository.

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

4. Check that [`src/dftbplus.h`](src/dftbplus.h)
   matches the DFTB+ API header you built against.

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
- The fourth argument must be the DFTB+ HSD input file.
- The fifth argument must currently be `NULL`.
- `fix_modify ... energy yes` is used in the example so the DFTB+
  contribution is included in the reported energy.

## Example: Graphene MD

The example directory [`examples/md`](examples/md)
contains:

- `md.in`: LAMMPS input script
- `graphene.data`: LAMMPS data file
- `graphene.gen`: DFTB+ geometry file
- `dftb_in.hsd`: DFTB+ input
- `C-C.skf`: Slater-Koster parameter file

Important:
- The atom order in `graphene.data` and `graphene.gen` must match.
- Review `dftb_in.hsd` and update it for your material system and SK
  files before production use.

Suggested run command:

```bash
export OMP_NUM_THREADS=1
mpirun -np 1 lmp_mpi -in md.in
```

Use a single MPI rank unless you modify the code to support parallel
execution.

## Reported HPC Build Notes

The original README included a user report from Iker Ortiz de Luzuriaga
at DIPC. The main points were:

- LAMMPS version used: `stable_29Oct2020`
- DFTB+ version used: `20.1`
- A forked DFTB+ build was also used
- On that system, one build provided `libmudpack.a` and the other
  provided a complete `libdftbplus.so`
- LAMMPS linker flags were adjusted manually to combine the required
  libraries

If your environment already builds the forked DFTB+ cleanly and exposes
the needed API library, that is the simpler path.

## Contact

Bug reports and comments can be sent to:
- `jungg@ornl.gov`
- `gs4phone@gmail.com`
