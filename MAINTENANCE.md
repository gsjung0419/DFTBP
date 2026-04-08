# Maintenance Guide

This repository now has two separate documentation tracks:

- `main`: current porting work for newer LAMMPS
- `legacy/lammps-29Oct2020`: frozen legacy baseline

## Recommended Version Policy

Treat the following stack as the supported baseline:

- LAMMPS: `29Oct2020`
- DFTB+: author's fork at
  <https://github.com/gsjung0419/dftbplus>
- Upstream DFTB+ generation referenced in earlier notes: `20.1`

If another version combination is made to work, record it explicitly
instead of silently broadening compatibility claims.

## Recommended Branch Strategy

Use the following split:

- `main`: latest porting work and reusable install automation
- `legacy/lammps-29Oct2020`: frozen maintenance branch for the old stack
- tags: known-good snapshots on either branch

Suggested tag style:

- `legacy-lammps29Oct2020-dftbplus20.1-v1`
- `port-lammps22Jul2025-dftbplus25.1-v1`

This keeps the older verified stack readable while allowing `main` to
move toward newer LAMMPS releases.

## How Users Should Install It

Version selection should happen at clone / checkout time, not during the
LAMMPS install step.

- If a user wants the current porting branch, a normal clone is enough.
- If a user wants the old verified stack, they should clone the legacy
  branch or check out a legacy tag first.
- The LAMMPS install step itself does not choose between legacy and
  modern code paths in this repository.

Examples:

```bash
git clone git@github.com:gsjung0419/DFTBP.git
```

If you want the legacy branch directly:

```bash
git clone -b legacy/lammps-29Oct2020 git@github.com:gsjung0419/DFTBP.git
```

Or clone first and then check out a known-good tag:

```bash
git clone git@github.com:gsjung0419/DFTBP.git
cd DFTBP
git checkout legacy-lammps29Oct2020-dftbplus20.1-v1
```

## Recommended Documentation Practice

Keep version information in four places:

- [`README.md`](README.md): current porting status and install flow
- [`README.legacy.md`](README.legacy.md): old-stack usage notes
- [`PORTING_LOG.md`](PORTING_LOG.md): chronological debug notes
- this file: maintenance policy and branch/tag rules

When you verify a working environment, record:

- exact LAMMPS version
- exact DFTB+ source or commit
- compiler version
- MPI version
- whether the build was serial-only

## Practical Advice

If you verify a new current stack, record:

1. exact LAMMPS version
2. exact DFTB+ release or commit
3. whether the DFTB+ library was static or shared
4. compiler / MPI / Fortran runtime details
5. whether `fix dftbp` remained serial-only
