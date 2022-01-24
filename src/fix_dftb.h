/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(dftbp,FixDFTBP)

#else

#ifndef LMP_FIX_DFTBP_H
#define LMP_FIX_DFTBP_H

#include "fix.h"
#include "dftbplus.h"

namespace LAMMPS_NS {

class FixDFTBP : public Fix {
public:
  FixDFTBP(class LAMMPS *, int, char **);
  virtual ~FixDFTBP();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void min_setup(int);
  void setup_pre_reverse(int, int);
  void initial_integrate(int);
  void pre_reverse(int, int);
  void post_force(int);
  double compute_vector(int);
  
  void min_post_force(int);
  void final_integrate();
  void reset_dt();
  double compute_scalar();
  double memory_usage();

  void write_restart(FILE *);
  void restart(char *);

protected:
  char *id_pe;
  char *infilename;
  int coulomb,pbcflag,pe_peratom,virial_global,virial_peratom,neighflag;
  int eflag_caller;
  
  int nmax,newsystem;
  double *qpotential;
  //double **flatte;
  double **fdftbp;
  //double dftbp_energy;  
  double dftbp_energy,dftbp_repulsiveE,dftbp_electronicE;

  
  DftbPlusInput *inputs;
  DftbPlus *dftbplus;

  double mermin_energy, mermin_energy_total,r_energy,e_energy;
  
  double const sunitconv=1.0/0.367493245336341E-01;
  double const funitconv=1.0/0.194469064593167E-01;
  double const eunitconv=1.0/0.367493245336341E-01;
  double const lunitconv=0.188972598857892E+01;
  double *gradients, *gradients_total, *gross_charges, *gross_charges_total;
  double *latvecs;
  double *fcoords;
  double *fstress;

  class NeighList *list;
  class Compute *c_pe;

  double energy_all[9];
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Must use units metal with fix dftbp command

UNDOCUMENTED

E: Fix dftbp currently runs only in serial

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix latte does not yet support a LAMMPS calculation of a Coulomb potential

UNDOCUMENTED

E: Could not find fix latte compute ID

UNDOCUMENTED

E: Fix dftbp compute ID does not compute pe/atom

UNDOCUMENTED

E: Fix dftbp requires 3d problem

UNDOCUMENTED

E: Fix dftbp cannot compute Coulomb potential

UNDOCUMENTED

E: Fix dftbp requires 3d simulation

UNDOCUMENTED

W: Fix dftbp should come after all other integration fixes

UNDOCUMENTED

E: Internal DFTB+ problem

UNDOCUMENTED

*/
