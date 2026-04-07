# Maintenance Guide

This repository is best maintained as a legacy package tied to an older
LAMMPS / DFTB+ combination.

## Recommended Version Policy

Treat the following stack as the supported baseline:

- LAMMPS: `29Oct2020`
- DFTB+: author's fork at
  <https://github.com/gsjung0419/dftbplus>
- Upstream DFTB+ generation referenced in earlier notes: `20.1`

If another version combination is made to work, record it explicitly
instead of silently broadening compatibility claims.

## Recommended Branch Strategy

If this repository will remain legacy-only for now, the simplest option
is:

- Keep `main` as the stable legacy branch
- Make only bug fixes or documentation updates on `main`
- Tag known-good snapshots

Suggested tag style:

- `legacy-lammps29Oct2020-dftbplus20.1-v1`
- `legacy-lammps29Oct2020-dftbplus20.1-v2`

If you plan to attempt a port to newer LAMMPS or newer DFTB+, use this
split instead:

- `main`: future porting work
- `legacy/lammps-29Oct2020`: frozen maintenance branch for the old stack

That avoids mixing compatibility fixes for the old code with API porting
changes for newer versions.

## How Users Should Install It

Version selection should happen at clone / checkout time, not during the
LAMMPS install step.

- If `main` is the legacy-stable branch, a normal clone is enough.
- If a separate legacy branch exists later, users should clone that
  branch or check out the legacy tag first.
- The LAMMPS install step itself does not choose between legacy and
  modern code paths in this repository.

Examples:

```bash
git clone git@github.com:gsjung0419/DFTBP.git
```

If you later split branches:

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

Keep version information in two places:

- [`README.md`](README.md): one short "Supported Legacy Stack" section
- this file: maintenance policy and branch/tag rules

When you verify a working environment, record:

- exact LAMMPS version
- exact DFTB+ source or commit
- compiler version
- MPI version
- whether the build was serial-only

## Practical Advice

If your real goal is only to preserve the currently working code, do not
start with a new porting branch yet. Instead:

1. Keep `main` as legacy-stable.
2. Add a tag to the last known-good state.
3. Only create a separate modernization branch when you are ready to
   test newer LAMMPS or DFTB+ APIs.

That is usually the lowest-maintenance setup.
