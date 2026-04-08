#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


class PatchError(RuntimeError):
    pass


@dataclass
class FilePlan:
    path: Path
    relpath: str
    original: str
    updated: str
    changed: bool


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def replace_once(text: str, old: str, new: str, marker: str, label: str, path: Path) -> tuple[str, bool]:
    if marker in text:
        return text, False
    if old not in text:
        raise PatchError(f"{path}: expected block not found for {label}")
    return text.replace(old, new, 1), True


def patch_dftbplus_header(text: str, path: Path) -> tuple[str, bool]:
    changed = False

    old = """#ifdef __cplusplus
extern \"C\" {
#endif
"""
    new = """#ifdef __cplusplus
extern \"C\" {
  typedef bool _Bool;
#else
  #include <stdbool.h>
#endif
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "#include <stdbool.h>",
        "stdbool compatibility",
        path,
    )
    changed |= step_changed

    old = """void dftbp_get_energy(DftbPlus *instance, double *mermin_energy);


/**
 * Queries the gradients of the current geometry.
"""
    new = """void dftbp_get_energy(DftbPlus *instance, double *mermin_energy);


/**
 * Queries the individual energy contributions of the current geometry.
 *
 * \\param[inout] instance Handler of the DFTB+ instance.
 *
 * \\param[out] mermin_energy Mermin free energy of the current geometry.
 * \\param[out] repulsive_energy Repulsive energy contribution.
 * \\param[out] electronic_energy Electronic energy contribution.
 */
void dftbp_get_energyparts(DftbPlus *instance, double *mermin_energy,
    double *repulsive_energy, double *electronic_energy);


/**
 * Queries the gradients of the current geometry.
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "void dftbp_get_energyparts(",
        "energy parts declaration",
        path,
    )
    changed |= step_changed

    old = """void dftbp_get_stress_tensor(DftbPlus *instance, double *stresstensor);


/**
 * Queries the net populations on the atoms.
"""
    new = """void dftbp_get_stress_tensor(DftbPlus *instance, double *stresstensor);


/**
 * Queries the virial tensor of the current periodic box.
 *
 * \\param[inout] instance Handler of the DFTB+ instance.
 *
 * \\param[out] virial Virial tensor for the periodic box. Shape: [3, 3].
 */
void dftbp_get_virial(DftbPlus *instance, double *virial);


/**
 * Queries the net populations on the atoms.
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "void dftbp_get_virial(",
        "virial declaration",
        path,
    )
    changed |= step_changed

    return text, changed


def patch_capi(text: str, path: Path) -> tuple[str, bool]:
    old = """  end subroutine handleOutputFile_


end module dftbp_capi
"""
    new = """  end subroutine handleOutputFile_


  !> Obtain the DFTB+ energy parts (total, repulsive, electronic)
  subroutine c_DftbPlus_getEnergyParts(handler, merminEnergy, repulsiveEnergy, electronicEnergy) &
      & bind(C, name='dftbp_get_energyparts')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Resulting energies
    real(c_double), intent(out) :: merminEnergy
    real(c_double), intent(out) :: repulsiveEnergy
    real(c_double), intent(out) :: electronicEnergy

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call instance%getEnergyParts(merminEnergy, repulsiveEnergy, electronicEnergy)

  end subroutine c_DftbPlus_getEnergyParts


  !> Obtain the virial tensor (3x3)
  subroutine c_DftbPlus_getVirial(handler, virial) bind(C, name='dftbp_get_virial')

    !> Handler for the calculation
    type(c_DftbPlus), intent(inout) :: handler

    !> Virial tensor
    real(c_double), intent(out) :: virial(3, 3)

    type(TDftbPlusC), pointer :: instance

    call c_f_pointer(handler%instance, instance)
    call instance%getVirial(virial(:, 1:3))

  end subroutine c_DftbPlus_getVirial


end module dftbp_capi
"""
    text, changed = replace_once(
        text,
        old,
        new,
        "bind(C, name='dftbp_get_energyparts')",
        "C API bindings",
        path,
    )
    return text, changed


def patch_mmapi(text: str, path: Path) -> tuple[str, bool]:
    changed = False

    old = """  use dftbp_dftbplus_mainapi, only : checkSpeciesNames, doOneTdStep, finalizeTimeProp,&
      & getAtomicMasses, getCM5Charges, getCutOff, getElStatPotential, getEnergy,&
      & getExtChargeGradients, getGradients, getGrossCharges, getLocalKS, getRefCharges,&
      & getStressTensor, getTdForces, initializeTimeProp, nrOfAtoms, nrOfKPoints, nrOfLocalKS,&
      & nrOfSpin, setExternalCharges, setExternalEfield, setExternalPotential, setGeometry,&
      & setNeighbourList, setQDepExtPotProxy, setRefCharges, setTdCoordsAndVelos,&
      & setTdElectricField, updateDataDependentOnSpeciesOrdering
"""
    new = """  use dftbp_dftbplus_mainapi, only : checkSpeciesNames, doOneTdStep, finalizeTimeProp,&
      & getAtomicMasses, getCM5Charges, getCutOff, getElStatPotential, getEnergy,&
      & getEnergyParts, getExtChargeGradients, getGradients, getGrossCharges, getLocalKS,&
      & getRefCharges, getStressTensor, getTdForces, getVirial, initializeTimeProp,&
      & nrOfAtoms, nrOfKPoints, nrOfLocalKS, nrOfSpin, setExternalCharges,&
      & setExternalEfield, setExternalPotential, setGeometry, setNeighbourList,&
      & setQDepExtPotProxy, setRefCharges, setTdCoordsAndVelos, setTdElectricField,&
      & updateDataDependentOnSpeciesOrdering
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "getEnergyParts, getExtChargeGradients",
        "main API imports",
        path,
    )
    changed |= step_changed

    old = """    !> Obtain the DFTB+ energy
    procedure :: getEnergy => TDftbPlus_getEnergy
    !> Obtain the DFTB+ gradients
    procedure :: getGradients => TDftbPlus_getGradients
    !> Obtain the DFTB+ stress tensor
    procedure :: getStressTensor => TDftbPlus_getStressTensor
"""
    new = """    !> Obtain the DFTB+ energy
    procedure :: getEnergy => TDftbPlus_getEnergy
    !> Obtain the DFTB+ energy parts (total/repulsive/electronic)
    procedure :: getEnergyParts => TDftbPlus_getEnergyParts
    !> Obtain the DFTB+ gradients
    procedure :: getGradients => TDftbPlus_getGradients
    !> Obtain the DFTB+ stress tensor
    procedure :: getStressTensor => TDftbPlus_getStressTensor
    !> Obtain the DFTB+ virial tensor
    procedure :: getVirial => TDftbPlus_getVirial
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "procedure :: getEnergyParts => TDftbPlus_getEnergyParts",
        "type-bound procedures",
        path,
    )
    changed |= step_changed

    old = """  end subroutine TDftbPlus_getEnergy


  !> Returns the gradient on the atoms in the system.
"""
    new = """  end subroutine TDftbPlus_getEnergy


  !> Return the energy parts of the current system.
  subroutine TDftbPlus_getEnergyParts(this, merminEnergy, repulsiveEnergy, electronicEnergy)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Mermin free energy.
    real(dp), intent(out) :: merminEnergy

    !> Repulsive energy contribution.
    real(dp), intent(out) :: repulsiveEnergy

    !> Electronic energy contribution.
    real(dp), intent(out) :: electronicEnergy

    call this%checkInit()

    call getEnergyParts(this%env, this%main, merminEnergy, repulsiveEnergy, electronicEnergy)

  end subroutine TDftbPlus_getEnergyParts


  !> Returns the gradient on the atoms in the system.
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "subroutine TDftbPlus_getEnergyParts(",
        "energy parts wrapper",
        path,
    )
    changed |= step_changed

    old = """  end subroutine TDftbPlus_getStressTensor


  !> Returns the gradients on the external charges.
"""
    new = """  end subroutine TDftbPlus_getStressTensor


  !> Returns the virial tensor of the periodic system.
  subroutine TDftbPlus_getVirial(this, virial)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Virial tensor.
    real(dp), intent(out) :: virial(:,:)

    call this%checkInit()

    call getVirial(this%env, this%main, virial)

  end subroutine TDftbPlus_getVirial


  !> Returns the gradients on the external charges.
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "subroutine TDftbPlus_getVirial(",
        "virial wrapper",
        path,
    )
    changed |= step_changed

    return text, changed


def patch_mainapi(text: str, path: Path) -> tuple[str, bool]:
    changed = False

    old = """  public :: checkSpeciesNames, doOneTdStep, finalizeTimeProp, getAtomicMasses, getCM5Charges
  public :: getCutOff, getElStatPotential, getEnergy, getExtChargeGradients, getGradients
  public :: getGrossCharges, getLocalKS, getRefCharges, getStressTensor, getTdForces
  public :: initializeTimeProp, nrOfAtoms, nrOfKPoints, nrOfLocalKS, nrOfSpin, setElectronNumber
  public :: setExternalCharges, setExternalEfield, setExternalPotential, setGeometry
  public :: setNeighbourList, setQDepExtPotProxy, setRefCharges, setTdCoordsAndVelos
  public :: setTdElectricField, updateDataDependentOnSpeciesOrdering
"""
    new = """  public :: checkSpeciesNames, doOneTdStep, finalizeTimeProp, getAtomicMasses, getCM5Charges
  public :: getCutOff, getElStatPotential, getEnergy, getEnergyParts, getExtChargeGradients
  public :: getGradients, getGrossCharges, getLocalKS, getRefCharges, getStressTensor
  public :: getTdForces, getVirial, initializeTimeProp, nrOfAtoms, nrOfKPoints
  public :: nrOfLocalKS, nrOfSpin, setElectronNumber, setExternalCharges
  public :: setExternalEfield, setExternalPotential, setGeometry, setNeighbourList
  public :: setQDepExtPotProxy, setRefCharges, setTdCoordsAndVelos
  public :: setTdElectricField, updateDataDependentOnSpeciesOrdering
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "getEnergy, getEnergyParts, getExtChargeGradients",
        "public exports",
        path,
    )
    changed |= step_changed

    old = """  end subroutine getEnergy


  !> get forces on atoms
"""
    new = """  end subroutine getEnergy


  !> get total, repulsive, and electronic energies
  subroutine getEnergyParts(env, main, merminEnergy, repulsiveEnergy, electronicEnergy)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> Resulting energies
    real(dp), intent(out) :: merminEnergy
    real(dp), intent(out) :: repulsiveEnergy
    real(dp), intent(out) :: electronicEnergy

    integer :: iDet

    call recalcGeometry(env, main)
    iDet = size(main%dftbEnergy)
    merminEnergy = main%dftbEnergy(iDet)%EMermin
    repulsiveEnergy = main%dftbEnergy(iDet)%Erep
    electronicEnergy = main%dftbEnergy(iDet)%Eelec

  end subroutine getEnergyParts


  !> get forces on atoms
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "subroutine getEnergyParts(env, main, merminEnergy, repulsiveEnergy, electronicEnergy)",
        "energy parts implementation",
        path,
    )
    changed |= step_changed

    old = """  end subroutine getStressTensor


  !> get the gross (Mulliken projected) charges for atoms wrt neutral atoms
"""
    new = """  end subroutine getStressTensor


  !> get virial tensor for unit cell
  subroutine getVirial(env, main, virial)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> resulting virial tensor
    real(dp), intent(out) :: virial(:,:)

    if (.not. main%tStress) then
      call error("Virial not available, you must initialise your calculator with&
          & this property enabled.")
    end if

    call recalcGeometry(env, main)
    virial(:,:) = main%totalStress(:,:) * main%cellVol

  end subroutine getVirial


  !> get the gross (Mulliken projected) charges for atoms wrt neutral atoms
"""
    text, step_changed = replace_once(
        text,
        old,
        new,
        "subroutine getVirial(env, main, virial)",
        "virial implementation",
        path,
    )
    changed |= step_changed

    return text, changed


PATCHERS = {
    "src/dftbp/api/mm/capi.F90": patch_capi,
    "src/dftbp/api/mm/dftbplus.h": patch_dftbplus_header,
    "src/dftbp/api/mm/mmapi.F90": patch_mmapi,
    "src/dftbp/dftbplus/mainapi.F90": patch_mainapi,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply the DFTB+ 25.1 MM API patch required by the DFTBP LAMMPS integration."
    )
    parser.add_argument(
        "--dftbplus-src",
        required=True,
        help="Path to the DFTB+ source tree to patch",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would change without writing files",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not save original files under .dftbplus_25patch_backup",
    )
    return parser.parse_args()


def make_backup(path: Path, source_root: Path, backup_root: Path) -> None:
    backup_path = backup_root / path.relative_to(source_root)
    backup_path.parent.mkdir(parents=True, exist_ok=True)
    if not backup_path.exists():
        shutil.copy2(path, backup_path)


def build_plan(source_root: Path) -> list[FilePlan]:
    plan: list[FilePlan] = []

    for relpath, patcher in PATCHERS.items():
        path = source_root / relpath
        if not path.exists():
            raise PatchError(f"Required file not found: {path}")

        original = read_text(path)
        updated, changed = patcher(original, path)
        plan.append(
            FilePlan(
                path=path,
                relpath=relpath,
                original=original,
                updated=updated,
                changed=changed,
            )
        )

    return plan


def main() -> int:
    args = parse_args()
    source_root = Path(args.dftbplus_src).resolve()

    if not source_root.exists():
        raise SystemExit(f"DFTB+ source tree does not exist: {source_root}")

    plan = build_plan(source_root)
    changed_files = [item for item in plan if item.changed]

    if not changed_files:
        print("DFTB+ 25.1 patch already applied; no files changed.")
        return 0

    print(f"Target source tree: {source_root}")
    for item in plan:
        status = "PATCH" if item.changed else "SKIP"
        print(f"[{status}] {item.relpath}")

    if args.dry_run:
        print("Dry run only; no files were modified.")
        return 0

    backup_root = source_root / ".dftbplus_25patch_backup"
    if not args.no_backup:
        for item in changed_files:
            make_backup(item.path, source_root, backup_root)
        print(f"Backups written under: {backup_root}")

    for item in changed_files:
        write_text(item.path, item.updated)

    print("Patch applied successfully.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except PatchError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
