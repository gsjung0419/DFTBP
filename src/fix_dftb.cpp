/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: GS JUNG (ORNL)
------------------------------------------------------------------------- */

#include "fix_dftb.h"
#include <cstdio>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "modify.h"
#include "compute.h"
#include "memory.h"
#include "error.h"

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

//It is from latte
/*extern "C" {
  void dftbp(int *, int *, double *, int *, int *,
             double *, double *, double *, double *,
             double *, double *, double *, int *,
             double *, double *, double *, double *, int * , bool *);
  //int dftbp_abiversion();
  }*/

// the ABIVERSION number here must be kept consistent
// with its counterpart in the DFTBP library and the
// prototype above. We want to catch mismatches with
// a meaningful error messages, as they can cause
// difficult to debug crashes or memory corruption.

#define DFTBP_ABIVERSION 20180622
#define INVOKED_PERATOM 8

/* ---------------------------------------------------------------------- */

FixDFTBP::FixDFTBP(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Must use units metal with fix dftbp command");
  
  if (comm->nprocs != 1)
    error->all(FLERR,"Fix dftbp currently runs only in serial");

  if (narg != 5) error->all(FLERR,"Illegal fix dftbp command");

  //scalar_flag = 1;
  //global_freq = 1;
  //extscalar = 1;
  virial_flag = 1;
  thermo_virial = 1;
  vector_flag = 1;
  size_vector = 9;
  // store ID of compute pe/atom used to generate Coulomb potential for DFTBP
  // NULL means DFTBP will compute Coulombic potential

  coulomb = 0;
  dftbplus=NULL;
  inputs=NULL;
  id_pe = NULL;
  infilename=NULL;
  fcoords=NULL;
  latvecs=NULL;
  gradients = NULL;
  fstress =NULL;
  qpotential = NULL;
  fdftbp = NULL;
  
  if (strcmp(arg[3],"NULL") == 0) {
    error->all(FLERR,"Fix dftbp requires dftb_in.hsd file and SK files");
  }else{
    int n = strlen(arg[3]) + 1;
    infilename = new char[n];
    strcpy(infilename,arg[3]);
  }
  
  if (strcmp(arg[4],"NULL") != 0) {
    coulomb = 1;
    error->all(FLERR,"Fix dftbp does not yet support a LAMMPS calculation "
               "of a Coulomb potential");

    int n = strlen(arg[4]) + 1;
    id_pe = new char[n];
    strcpy(id_pe,arg[4]);

    int ipe = modify->find_compute(id_pe);
    if (ipe < 0) error->all(FLERR,"Could not find fix dftbp compute ID");
    if (modify->compute[ipe]->peatomflag == 0)
    error->all(FLERR,"Fix dftbp compute ID does not compute pe/atom");
  }

  //If this constructor moved to init() function, minimization does not work..why?
  //ans:
  inputs=new DftbPlusInput;
  dftbplus=new DftbPlus;
  char outfilename[]="output";
  //Initializes a DFTB+ calculator
  dftbp_init(dftbplus,outfilename);
  //Fills up a DFTB+ input tree from a HSD input file
  dftbp_get_input_from_file(dftbplus,infilename,inputs);
  dftbp_process_input(dftbplus,inputs);

  int nAtoms = dftbp_get_nr_atoms(dftbplus);
  latvecs=new double[3*3];
  fstress=new double[3*3];
  gradients = new double[3*nAtoms];
  fcoords = new double[3*nAtoms];
  
  //std::cout<<"DFTB Debug: FixDFTB constructor called: "<<std::endl;  
  
  nmax = 0;
  dftbp_energy = 0.0;

  
}

/* ---------------------------------------------------------------------- */

FixDFTBP::~FixDFTBP()
{
  delete [] id_pe;
  delete [] infilename;  
  dftbp_final(dftbplus);
  delete inputs;
  delete [] latvecs;
  delete [] fstress;
  delete [] fcoords;
  delete [] gradients;
  //free(fcoords);
  //free(gradients);
  
  memory->destroy(qpotential);
  memory->destroy(fdftbp);
  //std::cout<<"DFTB Debug: FixDFTB destructor called: "<<std::endl;  
}

/* ---------------------------------------------------------------------- */

int FixDFTBP::setmask()
{
  int mask = 0;
  //mask |= INITIAL_INTEGRATE;
  //mask |= FINAL_INTEGRATE;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;  
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDFTBP::init()
{
  // error checks

  if (domain->dimension == 2)
    error->all(FLERR,"Fix dftbp requires 3d problem");

  if (coulomb) {
    if (atom->q_flag == 0 || force->pair == NULL || force->kspace == NULL)
      error->all(FLERR,"Fix dftbp cannot compute Coulomb potential");

    int ipe = modify->find_compute(id_pe);
    if (ipe < 0) error->all(FLERR,"Could not find fix dftbp compute ID");
    c_pe = modify->compute[ipe];
  }

  // must be fully periodic or fully non-periodic

  if (domain->nonperiodic == 0) pbcflag = 1;
  else if (!domain->xperiodic && !domain->yperiodic && !domain->zperiodic)
    pbcflag = 0;
  else error->all(FLERR,"Fix dftbp requires 3d simulation");

  // create fdftbp charges if needed 
  
  // for now, assume nlocal will never change
  
  if (coulomb && qpotential == NULL) {
    memory->create(qpotential,atom->nlocal,"dftbp:qpotential");
    memory->create(fdftbp,atom->nlocal,3,"dftbp:fdftbp");
  }

  //gradients = (double *) calloc(3 * nAtoms, sizeof(double));
  //fcoords= (double *) calloc(3 * nAtoms, sizeof(double));

  //std::cout<<"DFTB Debug: FixDFTB init() called: "<<std::endl;  

}

/* ---------------------------------------------------------------------- */

void FixDFTBP::init_list(int /*id*/, NeighList * /*ptr*/)
{
  // list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixDFTBP::setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
  //std::cout<<"DFTB Debug: FixDFTB setup() called: "<<std::endl;  
  
}

/* ---------------------------------------------------------------------- */

void FixDFTBP::min_setup(int vflag)
{
  newsystem = 1;
  post_force(vflag);
  newsystem = 0;
  //std::cout<<"DFTB Debug: FixDFTB min_setup() called: "<<std::endl;  
}

/* ---------------------------------------------------------------------- */

void FixDFTBP::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
  //std::cout<<"DFTB Debug: FixDFTB setup_pre_reverse() called: "<<std::endl;  
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixDFTBP::initial_integrate(int /*vflag*/) {}

/* ----------------------------------------------------------------------
   store eflag, so can use it in post_force to tally per-atom energies
------------------------------------------------------------------------- */

void FixDFTBP::pre_reverse(int eflag, int /*vflag*/)
{
  eflag_caller = eflag;
  //std::cout<<"DFTB Debug: FixDFTB pre_reverse() called: "<<std::endl;    
}

/* ---------------------------------------------------------------------- */

void FixDFTBP::post_force(int vflag){
  //std::cout<<"DFTB Debug: FixDFTB post_force() called: "<<std::endl;
  //std::cout<<"eflag and vflag: "<<eflag_caller<<" "<<vflag<<std::endl;  


  bigint ntimestep = update->ntimestep;
  
  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // hardwire these unsupported flags for now

  int coulombflag = 0;
  pe_peratom = 0;
  virial_global = 1;              // set via vflag_global at some point
  virial_peratom = 0;
  neighflag = 0;
  //set flags used by DFTBP
  // NOTE: DFTBP does not compute per-atom energies or virials
  int flags[6];

  flags[0] = pbcflag;         // 1 for fully periodic, 0 for fully non-periodic
  flags[1] = coulombflag;     // 1 for LAMMPS computes Coulombics, 0 for DFTBP
  flags[2] = eflag_atom;      // 1 to return per-atom energies, 0 for no
  flags[3] = vflag_global && thermo_virial;    // 1 to return global/per-atom
  flags[4] = vflag_atom && thermo_virial;      //   virial, 0 for no
  flags[5] = neighflag;       // 1 to pass neighbor list to DFTBP, 0 for no

  // setup DFTBP arguments
  int nAtoms = atom->nlocal;
  double *coords = &atom->x[0][0];
  int *type = atom->type;
  int ntypes = atom->ntypes;
  double *mass = &atom->mass[1];
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *forces;
  bool dftbperror = 0;
  //if (coulomb) forces = &fdftbp[0][0];
  forces = &atom->f[0][0];
  int maxiter = -1;
  
  latvecs[0]=(boxhi[0]-boxlo[0])*lunitconv;
  //latvecs[1]=domain->xy*1.889725989;
  //latvecs[2]=domain->xz*1.889725989;
  //latvecs[3]=0.0;
  latvecs[1]=0.0;
  latvecs[2]=0.0;
  latvecs[3]=domain->xy*lunitconv;

  latvecs[4]=(boxhi[1]-boxlo[1])*lunitconv;

  //latvecs[5]=domain->yz*1.889725989;  
  //latvecs[6]=0.0;
  //latvecs[7]=0.0;
  latvecs[5]=0.0;
  latvecs[6]=domain->xz*lunitconv;
  latvecs[7]=domain->yz*lunitconv;  

  latvecs[8]=(boxhi[2]-boxlo[2])*lunitconv;
  
  for (int i=0;i<nAtoms;i++){
    for(int j=0;j<3;j++){
      fcoords[i*3+j]=coords[i*3+j]*lunitconv;      
    }
  }

  //std::cout<<"DFTB Debug: check point after coordinate set: "<<nAtoms <<std::endl;    
  //dftbp_set_coords(dftbplus,fcoords);

  dftbp_set_coords_and_lattice_vecs(dftbplus,fcoords,latvecs);  

  //std::cout<<"DFTB Debug: check point after coordinate set: "<<std::endl;  

  //energy update only eflag>1
  if(eflag){
    //dftbp_get_energy(dftbplus, &mermin_energy);
    dftbp_get_energyparts(dftbplus, &mermin_energy,&r_energy,&e_energy);    
    dftbp_energy=mermin_energy*eunitconv;
    dftbp_repulsiveE=r_energy*eunitconv;
    dftbp_electronicE=e_energy*eunitconv;
   }

  //std::cout<<"DFTB Debug: check point after energy called: "<<std::endl;    
  //Force should be always updated! 

  
  dftbp_get_gradients(dftbplus,gradients);

  //std::cout<<"DFTB Debug: check point after graidient called: "<<std::endl;      
  for (int i=0;i<nAtoms;i++){
    for(int j=0;j<3;j++){
      forces[i*3+j]-=gradients[i*3+j]*funitconv;
    }
  }
  
  //std::cout<<"DFTB Debug: check point after force set: "<<std::endl;    
  //energy update only vflag>1  
  //if(vflag>0){
    dftbp_get_virial(dftbplus,fstress);
    virial[0]=fstress[0]*sunitconv;
    virial[1]=fstress[4]*sunitconv;
    virial[2]=fstress[8]*sunitconv;
    virial[3]=0.5*(fstress[1]+fstress[3])*sunitconv;
    virial[4]=0.5*(fstress[2]+fstress[6])*sunitconv;
    virial[5]=0.5*(fstress[5]+fstress[7])*sunitconv;      
    //}
    //std::cout<<"DFTB Debug: check point after virial set: "<<std::endl;    
  return;
}

/*void FixDFTBP::post_force(int vflag)
{

  int eflag = eflag_caller;
  ev_init(eflag,vflag);

  // compute Coulombic potential = pe[i]/q[i]
  // invoke compute pe/atom
  // wrap with clear/add and trigger pe/atom calculation every step

  if (coulomb) {
    modify->clearstep_compute();

    if (!(c_pe->invoked_flag & INVOKED_PERATOM)) {
      c_pe->compute_peratom();
      c_pe->invoked_flag |= INVOKED_PERATOM;
    }

    modify->addstep_compute(update->ntimestep+1);

    double *pe = c_pe->vector_atom;
    double *q = atom->q;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (q[i]) qpotential[i] = pe[i]/q[i];
      else qpotential[i] = 0.0;
  }

  // hardwire these unsupported flags for now

  int coulombflag = 0;
  // pe_peratom = 0;
  // virial_global = 1;              // set via vflag_global at some point
  // virial_peratom = 0;
  neighflag = 0;

  // set flags used by DFTBP
  // NOTE: DFTBP does not compute per-atom energies or virials

  int flags[6];

  flags[0] = pbcflag;         // 1 for fully periodic, 0 for fully non-periodic
  flags[1] = coulombflag;     // 1 for LAMMPS computes Coulombics, 0 for DFTBP
  flags[2] = eflag_atom;      // 1 to return per-atom energies, 0 for no
  flags[3] = vflag_global && thermo_virial;    // 1 to return global/per-atom
  flags[4] = vflag_atom && thermo_virial;      //   virial, 0 for no
  flags[5] = neighflag;       // 1 to pass neighbor list to DFTBP, 0 for no

  // setup DFTBP arguments

  int natoms = atom->nlocal;
  double *coords = &atom->x[0][0];
  int *type = atom->type;
  int ntypes = atom->ntypes;
  double *mass = &atom->mass[1];
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *forces;
  bool dftbperror = 0;
  if (coulomb) forces = &fdftbp[0][0];
  else forces = &atom->f[0][0];
  int maxiter = -1;

  //dftbp(flags,&natoms,coords,type,&ntypes,mass,boxlo,boxhi,&domain->xy,
  //      &domain->xz,&domain->yz,forces,&maxiter,&dftbp_energy,
  //      &atom->v[0][0],&update->dt,virial,&newsystem,&dftbperror);

  if (dftbperror) error->all(FLERR,"Internal DFTBP problem");

  // sum DFTBP forces to LAMMPS forces
  // e.g. LAMMPS may compute Coulombics at some point

  if (coulomb) {
    double **f = atom->f;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) {
      f[i][0] += fdftbp[i][0];
      f[i][1] += fdftbp[i][1];
      f[i][2] += fdftbp[i][2];
    }
  }
  }*/

/* ---------------------------------------------------------------------- */

void FixDFTBP::min_post_force(int vflag)
{
  //To perform energy minimization, eflag should be set > 0 
  eflag_caller=1;
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   integrate electronic degrees of freedom
------------------------------------------------------------------------- */

void FixDFTBP::final_integrate() {}

/* ---------------------------------------------------------------------- */

void FixDFTBP::reset_dt()
{
  //dtv = update->dt;
  //dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   DFTB energy from DFTBP
------------------------------------------------------------------------- */

double FixDFTBP::compute_scalar()
{
  return dftbp_energy;
}


void FixDFTBP::write_restart(FILE *fp){
#define RESTART_ITEM 6
  double buf[RESTART_ITEM+1];
  if(comm->me == 0){
    buf[0]=dftbp_energy;
    for (int i=0;i<6;i++){
      buf[i+1]=virial[i];
    } 
  }
}

void FixDFTBP::restart(char *buf){
  double *restored = (double *)buf;
  dftbp_energy=restored[0];
  for (int i=0;i<6;i++){
    virial[i]=restored[i+1];
  }
  delete [] restored;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double FixDFTBP::memory_usage()
{
  //std::cout<<"DFTB Debug: FixDFTB memory_usage() called: "<<std::endl;  
  
  double bytes = 0.0;
  if (coulomb) bytes += nmax * sizeof(double);
  if (coulomb) bytes += nmax*3 * sizeof(double);
  return bytes;
}


double FixDFTBP::compute_vector(int n)
{
  energy_all[0]=dftbp_energy;
  energy_all[1]=dftbp_repulsiveE;
  energy_all[2]=dftbp_electronicE;
  energy_all[3]=virial[0];
  energy_all[4]=virial[1];
  energy_all[5]=virial[2];
  energy_all[6]=virial[3];
  energy_all[7]=virial[4];
  energy_all[8]=virial[5];  
  
  return energy_all[n];
}
