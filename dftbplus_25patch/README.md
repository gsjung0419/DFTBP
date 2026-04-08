# DFTB+ 25.1 Patch Helper

This folder packages the DFTB+ source changes needed by the current
`DFTBP` integration on newer LAMMPS builds.

It is intended for a fresh DFTB+ `25.1` source tree. The helper patches
the MM API so that the following functions exist for `fix dftbp`:

- `dftbp_get_energyparts`
- `dftbp_get_virial`

## What It Changes

The script modifies these files inside a DFTB+ source tree:

- `src/dftbp/api/mm/capi.F90`
- `src/dftbp/api/mm/dftbplus.h`
- `src/dftbp/api/mm/mmapi.F90`
- `src/dftbp/dftbplus/mainapi.F90`

It does not build DFTB+ for you. It only patches the source tree.

## Usage

Apply the patch to a DFTB+ `25.1` source directory:

```bash
python3 /path/to/DFTBP/dftbplus_25patch/apply_patch.py \
  --dftbplus-src /path/to/dftbplus-25.1
```

Dry run:

```bash
python3 /path/to/DFTBP/dftbplus_25patch/apply_patch.py \
  --dftbplus-src /path/to/dftbplus-25.1 \
  --dry-run
```

By default, backups are written under:

```text
/path/to/dftbplus-25.1/.dftbplus_25patch_backup/
```

## After Patching

1. Build and install DFTB+ with API support enabled.
2. Point the LAMMPS-side DFTBP integration to that DFTB+ install.
3. Stage/build the LAMMPS package with
   [`tools/integrate_lammps_dftbp.py`](../tools/integrate_lammps_dftbp.py).

## Notes

- The helper is intentionally strict and expects the DFTB+ `25.1`
  source layout.
- If the expected code blocks are not found, the script stops instead of
  guessing.
- Re-running the script is safe. Already-patched files are detected and
  left unchanged.
