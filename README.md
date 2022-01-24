# DFTBP
LAMMPS codes for DFTB+
Please check DFTB+ code https://github.com/dftbplus/dftbplus

This package provides a fix dftbp command which is a wrapper on the
DFTB+ DFTB code, so that molecular dynamics can be run with LAMMPS
using density-functional tight-binding quantum forces calculated by
DFTB+.

GS JUNG@ORNL made this package based on fix_latte.cpp and fix_latte.h (developed by LANL)
Please report any bug/commetns to jungg@ornl.gov or gs4phone@gmail.com
The necessary functions for DFTB+ API have been implemented (update is allowed by dftb+)
Until that, you can check the details of changes in https://github.com/gsjung0419/dftbplus

0. The code is only tested with following environments.
 -. Intel compiler (19 and 20)
 -. MPICH 3.2.1 version compiled by intel compiler
 -. CMAKE 3.12 or 3.19
 -. Matching the versions of programs are recommended (The author did not check other possible combinations). 

1. DFTB+ options (tested) (See your confg.cmake in dftb+)
 -. WITH_OMP: TRUE or FALSE
 -. WITH_DFTD3: TRUE or FALSE
 -. WITH_API: TRUE (Mandatory for static library, e.g., lirbary.a)
 -. BUILD_SHARED_LIBS: TURE (Mandatory for shared library, e.g., lirbary.so, mandatory for python-based lammps users)
 -. You may have compile issues if you try with other environments.

2-1. Installation of DFTB+ (recommended for the most fundamental version)
 -. Compile the dftb+ (WITH_OMP, WITH_API: TRUE)
 -. Check libdftbplus.a is generated. 
 -. Locate libdftbplus.a in preffered folder (ex., ~/applic/lib/dftbplus)
 
2.2. Installation of LAMMPS (recommended for the most fundamental version)
 -. Use LAMMPS version 29OCT20, the newest version did not work for "fix" based modification. 
 -. Copy lib/dftbp to LAMMPS_sourcecode/lib
 -. Edit Makefile.lammps.mpi in lib/dftbp, "location of library" you set in DFTB+ installation (ex., ~/applic/lib/dftbplus).
 -. Move DFTBP into your lammps/src folder
 !! Make sure your dftbplus.h is the same as provided from dftbplus (there might be updates for API functions).
 -. Edit Makefile in lammps/src folder;include "dftbp" in PACKAGE.
 -. Check "make package-status" whether you have DFTBP package and
 -. Make available with "make yes-dftbp" 
 -. Compile lammps with "make mpi" in src folder (compile with intel compiler + mpich3.2.1)

3. Usage and Example
 -. 


