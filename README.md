# DFTBP
LAMMPS codes for DFTB+

Reference: Gang Seob Jung, Stephan Irle, and Bobby Sumpter. "Dynamic aspects of graphene deformation and fracture from approximate density functional theory." Carbon (2022)

Please check DFTB+ code https://github.com/dftbplus/dftbplus

This package provides a fix dftbp command which is a wrapper on the DFTB+ DFTB code, so that molecular dynamics can be run with LAMMPS
using density-functional tight-binding calculations by DFTB+.

GS JUNG@ORNL made this package based on fix_latte.cpp and fix_latte.h (developed by LANL)

Please report any bug/commetns to jungg@ornl.gov or gs4phone@gmail.com

The necessary functions for DFTB+ API have been implemented (update needs to be allowed by dftb+, please note that this forked version was not merged).

You should use the forked version in https://github.com/gsjung0419/dftbplus, or you should update the API in the recente version of DFTB+ (was not tested). 

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
  
-----------------------------------------------------------------------------------------------------
Notes from Iker Ortiz de Luzuriaga @ DIPC

Installation in HPC of the DIPC (Donostia International Physics Center).

Software involved:


LAMMPS version Stable 29 Oct 2020.

dftbplus, version 20.1.

dftbplus, forked version by the author of DFTBP.

Steps:

Download DFTBP:

 git clone https://github.com/gsjung0419/DFTBP

Download dftbplus:

 wget https://github.com/dftbplus/dftbplus/archive/refs/tags/20.1.tar.gz
 
Download forked dftbplus:

 git clone https://github.com/gsjung0419/dftbplus
 
Download LAMMPS:

 wget https://github.com/lammps/lammps/archive/refs/tags/stable_29Oct2020.tar.gz
 
tar zxvf stable_29Oct2020.tar.gz

Compile original dftbplus with foss/2023a toolchain:

 ml load foss/2023a
 
tar zxvf 20.1.tar.gz

cd 20.1

./utils/get_opt_externals ALL

mkdir _build && cd _build

cmake -DCMAKE_INSTALL_PREFIX=/scratch/$USER/opt/dftb+ -DCMAKE_TOOLCHAIN_FILE=../sys/intel.cmake -DWITH_API=TRUE -DWITH_OMP=TRUE -DBUILD_SHARED_LIBS=TRUE ..

make

make install

This will compile an incomplete version of libdftbplus.so and libmudpack.a, which is required later.

Compile forked dftbplus with intel/2023a:

 ml load intel/2023a
 
cd dftbplus

./utils/get_opt_externals ALL

mkdir _build && cd _build

cmake -DCMAKE_INSTALL_PREFIX=/scratch/$USER/opt/dftb+ -DCMAKE_TOOLCHAIN_FILE=../sys/intel.cmake -DWITH_API=TRUE -DWITH_OMP=TRUE -DBUILD_SHARED_LIBS=TRUE ..

make

make install

This will compile an complete version of libdftbplus.so but no libmudpack.a

Check that you have both files on /scratch/$USER/opt/dftb+

Copy lib/dftbp to LAMMPS_sourcecode/lib:

 cp -R DFTBP/lib/dftbp lammps-stable_29Oct2020/lib/

Copy src to lammps/src/DFTBP

 cp -R DFTBP/src lammps-stable_29Oct2020/src/DFTBP

Enter the LAMMPS directory:

 cd lammps-stable_29Oct2020

Edit Makefile.lammps.mpi in lib/dftbp. It must contain the following (possibly not required) line:

 dftbp_SYSLIB = -L/scratch/cperezmig/opt/dftb+/lib -ldftbplus

Edit Makefile in lammps/src folder;include "dftbp" in PACKAGE.

Install the package:

 make yes-dftbp

Check that the package has been included. Dftbp must appear as YES:

 make package-status

Modify src/MAKE/Makefile.mpi to change compilation and linker flags. Note the gcc and intel flags. We are compiling with both toolchains activated

 CCFLAGS =   -g -O3 -fopenmp

LINKFLAGS =	-g -O3 -fopenmp -lirc -L/scratch/cperezmig/opt/dftb+/lib -ldftbplus -lgfortran -lmpi -lmudpack -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -lifcore -limf -lsvml -fopenmp

Compile the mpi objective:

 cd src

make mpi

Now you should have two files: lmp_mpi and liblammps_mpi.a. Install them where you want.

...comments from GS: "if possible, just installation of forked version (DFTB+) is recommended"....
