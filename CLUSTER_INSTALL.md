# Cluster Install Recipe

This document records the final tested Linux cluster build flow for:

- `DFTB+ 25.1` with the API patch required by this repository
- `LAMMPS 22Jul2025`
- the `DFTBP` package from this repository
- Intel oneAPI compilers and Intel MPI

The final tested working combination is:

- `DFTB+ 25.1` built as `serial + OpenMP + API + s-dftd3`
- `LAMMPS 22Jul2025` built with `DFTBP + EXTRA-FIX`
- one MPI rank for LAMMPS and user-chosen OpenMP threads for DFTB+

## Assumptions

- Intel oneAPI is already installed on the cluster.
- Your environment provides:
  - `icx`
  - `icpx`
  - `ifx`
  - `mpiicx`
  - `mpiicpx`
  - `cmake`
  - `make`
  - `python3`
- `MKLROOT` and `I_MPI_ROOT` are available after you load the cluster modules.
- You are building on a login/build node, not inside a production queue job.

## Important Notes

- `fix dftbp` is still serial-only on the LAMMPS side.
- The tested working Linux path is not `MPI DFTB+ API`.
- If DFTB+ is built with `WITH_MPI=TRUE`, `fix dftbp` can fail at runtime with:
  `Cannot call MPI_INIT or MPI_INIT_THREAD more than once`
- The tested fix is to build DFTB+ with:
  - `WITH_MPI=FALSE`
  - `WITH_API=TRUE`
  - `WITH_SDFTD3=TRUE`
  - `BUILD_SHARED_LIBS=TRUE`
- `smd.in` from `TestDFTB` also needs the LAMMPS `EXTRA-FIX` package because it uses `fix smd`.
- A `DFTBP`-enabled LAMMPS executable requires the DFTB+ shared-library directory at runtime, not just at build time.
- This repository filters `NONE` and `-lNONE` out of `dftbplus.pc` link items inside `cmake/Modules/Packages/DFTBP.cmake`.

## Suggested Directory Layout

```bash
export WORKROOT=$HOME/build_dftbp_stack
mkdir -p "$WORKROOT"
cd "$WORKROOT"
```

## 1. Load the cluster toolchain

Replace these module commands with the ones used on your cluster.

```bash
module purge
module load intel-oneapi-compilers
module load intel-oneapi-mpi
module load intel-oneapi-mkl
module load cmake
module load python
```

Sanity checks:

```bash
which icx
which icpx
which ifx
which mpiicx
which mpiicpx
which cmake
echo "$MKLROOT"
echo "$I_MPI_ROOT"
```

## 2. Download the sources

```bash
cd "$WORKROOT"
curl -L -O https://github.com/dftbplus/dftbplus/releases/download/25.1/dftbplus-25.1.tar.xz
curl -L -O https://download.lammps.org/tars/lammps-22Jul2025.tar.gz
tar -xf dftbplus-25.1.tar.xz
tar -xf lammps-22Jul2025.tar.gz
git clone <YOUR-DFTBP-REPO-URL> DFTBP
```

If you already have this repository on the cluster, skip the clone and point the commands below at that path.

## 3. Patch DFTB+ 25.1

```bash
cd "$WORKROOT"
python3 DFTBP/dftbplus_25patch/apply_patch.py \
  --dftbplus-src "$WORKROOT/dftbplus-25.1"
```

## 4. Build and install DFTB+

Build DFTB+ as `serial + OpenMP + API + s-dftd3`.

```bash
cd "$WORKROOT"
cmake -S "$WORKROOT/dftbplus-25.1" \
  -B "$WORKROOT/dftbplus-25.1/_build_serial_omp_sdftd3" \
  -DCMAKE_INSTALL_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_CXX_COMPILER=icpx \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DWITH_MPI=FALSE \
  -DWITH_API=TRUE \
  -DBUILD_SHARED_LIBS=TRUE \
  -DENABLE_DYNAMIC_LOADING=FALSE \
  -DWITH_PYTHON=FALSE \
  -DWITH_SDFTD3=TRUE
cmake --build "$WORKROOT/dftbplus-25.1/_build_serial_omp_sdftd3" --parallel 8
cmake --install "$WORKROOT/dftbplus-25.1/_build_serial_omp_sdftd3"
```

Recommended runtime environment after install:

```bash
export DFTBPLUS_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3"
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:$LD_LIBRARY_PATH"
```

Shell startup reminder:

```bash
export DFTBPLUS_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3"
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:$LD_LIBRARY_PATH"
```

If `.bashrc` still points at an older MPI-enabled DFTB+ install, you can reintroduce the MPI double-init failure even when the new `lmp` binary is correct.

Quick check:

```bash
echo "$DFTBPLUS_PREFIX"
```

Optional quick check:

```bash
"$DFTBPLUS_PREFIX/bin/dftb+" --version
```

## 5. Stage the DFTBP package into LAMMPS and build

```bash
cd "$WORKROOT"
python3 "$WORKROOT/DFTBP/tools/integrate_lammps_dftbp.py" \
  --lammps-dir "$WORKROOT/lammps-22Jul2025" \
  --dftbplus-prefix "$DFTBPLUS_PREFIX" \
  --build \
  --build-system cmake \
  --build-dir "$WORKROOT/lammps-22Jul2025/_build_dftbp_sdftd3_extrafix" \
  --jobs 8 \
  --cmake-arg=-D --cmake-arg=CMAKE_C_COMPILER=mpiicx \
  --cmake-arg=-D --cmake-arg=CMAKE_CXX_COMPILER=mpiicpx \
  --cmake-arg=-D --cmake-arg=PKG_EXTRA-FIX=ON \
  --cmake-arg=-D --cmake-arg=BUILD_MPI=on
```

## 6. Runtime environment for LAMMPS

Before running `lmp`, keep the DFTB+ shared library visible:

```bash
export DFTBPLUS_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3"
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:$LD_LIBRARY_PATH"
export OMP_NUM_THREADS=<your_threads>
```

`OMP_NUM_THREADS=4` was only a test example. Use the thread count that matches your job allocation.

## 7. Minimal validation

Check that the executable exists:

```bash
ls -l "$WORKROOT/lammps-22Jul2025/_build_dftbp_sdftd3_extrafix/lmp"
```

Check that the DFTBP and EXTRA-FIX packages are compiled in:

```bash
"$WORKROOT/lammps-22Jul2025/_build_dftbp_sdftd3_extrafix/lmp" -h | rg 'DFTBP|EXTRA-FIX'
```

Run the tested example from this repository:

```bash
cd "$WORKROOT/DFTBP/TestDFTB"
"$WORKROOT/lammps-22Jul2025/_build_dftbp_sdftd3_extrafix/lmp" -in smd.in
```

## 8. Common Failure Modes

`cannot find -lNONE`

- Cause: `dftbplus.pc` can contain `NONE` in `Libs.private`.
- Status: handled in this repository by filtering `NONE` and `-lNONE` inside `cmake/Modules/Packages/DFTBP.cmake`.

`error while loading shared libraries: libdftbplus.so`

- Cause: `LD_LIBRARY_PATH` does not include the DFTB+ install `lib/`.
- Fix:

```bash
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:$LD_LIBRARY_PATH"
```

`error while loading shared libraries: libmpifx.so` or `libscalapackfx.so`

- Cause: those libraries came from an older MPI-enabled DFTB+ install.
- Fix: point `DFTBPLUS_PREFIX` and `LD_LIBRARY_PATH` to the serial/OpenMP install instead.

```bash
export DFTBPLUS_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3"
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:$LD_LIBRARY_PATH"
```

`Cannot call MPI_INIT or MPI_INIT_THREAD more than once`

- Cause: runtime environment still points at an MPI-enabled DFTB+ API build, while LAMMPS is already MPI-enabled.
- Fix: use the serial/OpenMP DFTB+ install and update shell startup files.

```bash
export DFTBPLUS_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3"
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:$LD_LIBRARY_PATH"
```

`DFTB+ compiled without support for s-dftd3 library`

- Cause: the input requests `s-dftd3`, but DFTB+ was built without `WITH_SDFTD3=TRUE`.
- Fix: rebuild DFTB+ with `-DWITH_SDFTD3=TRUE`.

`Unrecognized fix style 'smd'`

- Cause: the `EXTRA-FIX` package is not enabled in LAMMPS.
- Fix: rebuild LAMMPS with `-D PKG_EXTRA-FIX=ON`.

MPI wrapper mismatch

- Use `mpiicx` and `mpiicpx` consistently for the LAMMPS MPI build.
- Use `icx`, `icpx`, and `ifx` for the serial/OpenMP DFTB+ build.

Missing MKL or ScaLAPACK at configure time

- Make sure your oneAPI MKL module is loaded.
- Check:

```bash
echo "$MKLROOT"
```

## 9. Direct Commands Only

If you only want the command list without explanation:

```bash
export WORKROOT=$HOME/build_dftbp_stack
mkdir -p "$WORKROOT"
cd "$WORKROOT"

module purge
module load intel-oneapi-compilers
module load intel-oneapi-mpi
module load intel-oneapi-mkl
module load cmake
module load python

curl -L -O https://github.com/dftbplus/dftbplus/releases/download/25.1/dftbplus-25.1.tar.xz
curl -L -O https://download.lammps.org/tars/lammps-22Jul2025.tar.gz
tar -xf dftbplus-25.1.tar.xz
tar -xf lammps-22Jul2025.tar.gz
git clone <YOUR-DFTBP-REPO-URL> DFTBP

python3 DFTBP/dftbplus_25patch/apply_patch.py \
  --dftbplus-src "$WORKROOT/dftbplus-25.1"

cmake -S "$WORKROOT/dftbplus-25.1" \
  -B "$WORKROOT/dftbplus-25.1/_build_serial_omp_sdftd3" \
  -DCMAKE_INSTALL_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_CXX_COMPILER=icpx \
  -DCMAKE_Fortran_COMPILER=ifx \
  -DWITH_MPI=FALSE \
  -DWITH_API=TRUE \
  -DBUILD_SHARED_LIBS=TRUE \
  -DENABLE_DYNAMIC_LOADING=FALSE \
  -DWITH_PYTHON=FALSE \
  -DWITH_SDFTD3=TRUE
cmake --build "$WORKROOT/dftbplus-25.1/_build_serial_omp_sdftd3" --parallel 8
cmake --install "$WORKROOT/dftbplus-25.1/_build_serial_omp_sdftd3"

export DFTBPLUS_PREFIX="$WORKROOT/dftbplus-25.1/_install_serial_omp_sdftd3"
export LD_LIBRARY_PATH="$DFTBPLUS_PREFIX/lib:$LD_LIBRARY_PATH"

python3 "$WORKROOT/DFTBP/tools/integrate_lammps_dftbp.py" \
  --lammps-dir "$WORKROOT/lammps-22Jul2025" \
  --dftbplus-prefix "$DFTBPLUS_PREFIX" \
  --build \
  --build-system cmake \
  --build-dir "$WORKROOT/lammps-22Jul2025/_build_dftbp_sdftd3_extrafix" \
  --jobs 8 \
  --cmake-arg=-D --cmake-arg=CMAKE_C_COMPILER=mpiicx \
  --cmake-arg=-D --cmake-arg=CMAKE_CXX_COMPILER=mpiicpx \
  --cmake-arg=-D --cmake-arg=PKG_EXTRA-FIX=ON \
  --cmake-arg=-D --cmake-arg=BUILD_MPI=on

export OMP_NUM_THREADS=<your_threads>

cd "$WORKROOT/DFTBP/TestDFTB"
"$WORKROOT/lammps-22Jul2025/_build_dftbp_sdftd3_extrafix/lmp" -in smd.in
```
