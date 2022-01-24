# DFTBP
LAMMPS codes for DFTB+

Reference: Jung, Gang Seob, Stephan Irle, and Bobby Sumpter. "Dynamic aspects of graphene deformation and fracture from approximate density functional theory." Carbon (2022)

Please check DFTB+ code https://github.com/dftbplus/dftbplus

This package provides a fix dftbp command which is a wrapper on the DFTB+ DFTB code, so that molecular dynamics can be run with LAMMPS
using density-functional tight-binding calculations by DFTB+.

GS JUNG@ORNL made this package based on fix_latte.cpp and fix_latte.h (developed by LANL)

Please report any bug/commetns to jungg@ornl.gov or gs4phone@gmail.com

The necessary functions for DFTB+ API have been implemented (update is allowed by dftb+).

Until that, you can check the details of changes in https://github.com/gsjung0419/dftbplus

0. The code is only tested with following environments.

 -. Intel compiler (19 and 20)
 
 -. MPICH 3.2.1 version compiled by intel compiler
 
 -. CMAKE 3.12 or 3.19
 
 -. Matching the versions of programs are recommended (The author did not check other possible combinations). 


1. DFTB+ options (tested) (See your confg.cmake in dftb+)

 -. WITH_OMP: TRUE or FALSE
 
 -. WITH_DFTD3: TRUE or FALSE
 
 -. WITH_API: TRUE (Mandatory for static library, e.g., library.a)
 
 -. BUILD_SHARED_LIBS: TURE (Mandatory for shared library, e.g., library.so, mandatory for python-based lammps users)
 
 -. You may have compile issues if you try with other environments.


2.1. Installation of DFTB+ (recommended for the most fundamental version)

 -. Compile the dftb+ (WITH_OMP, WITH_API: TRUE)
 
 -. Check libdftbplus.a is generated.
 
 -. Locate libdftbplus.a in preffered folder (ex., ~/applic/lib/dftbplus)


2.2. Installation of LAMMPS (recommended for the most fundamental version)

 -. Use LAMMPS version 29OCT20, the newest version did not work for "fix" based modification.
 
 -. Copy lib/dftbp to LAMMPS_sourcecode/lib
 
 -. Edit Makefile.lammps.mpi in lib/dftbp, "location of library" you set in DFTB+ installation (ex., ~/applic/lib/dftbplus).
 
 -. Move src as lammps/src/DFTBP 
 
 !! Make sure your dftbplus.h is the same as provided from dftbplus (there might be updates for API functions).
 
 -. Edit Makefile in lammps/src folder;include "dftbp" in PACKAGE.
 
 -. Check "make package-status" whether you have DFTBP package
 
 -. Make available with "make yes-dftbp"
 
 -. Compile lammps with "make mpi" in src folder (compile with intel compiler + mpich3.2.1)

3 Usage and Example

 !! Before you are running the simulation. Check OMP_NUM_THREADS= 1 or specify a number
 
 !! Without specified the OMP number, DFTB+ tries to use all available cores, making it slow.
 

 3.1 Graphene md (example C-C.skf is MIO)
 
  -. data file for lammps and dftb+ should be prepared separately.
  
  -. The order of input atoms should be matched.
  
  -. Check dftb_in.hsd for options for dftb+. Most options are turned off except force calculations
  


