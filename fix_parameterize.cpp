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

#include <stdio.h>
#include <string.h>
#include "fix_parameterize.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

#include "random_mars.h" //random number generator, Marsaglia style
#include "pair_tersoff.h"
#include "pair_lj_cut_coul_inout.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixParameterize::FixParameterize(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5)
    error->all(FLERR,"Illegal fix parameterize command: requires input file (target forces) and random seed");
  //get inputs from user
  char* input_filename = arg[3];
  random_seed = force->inumeric(FLERR,arg[4]);
  //read file containing target forces (expecting to be in a flat list, one number on each line)
  FILE *input_file = fopen(input_filename, "r");
  if(input_file==NULL) error->all(FLERR,"No such target-forces file exists");
  char buffer[80];
  n_target_forces = 0;
  while(fscanf(input_file, "%s\n", &buffer) != EOF) //first pass: count how many forces there are
    n_target_forces++;
  rewind(input_file);
  target_forces = new double[n_target_forces]; //allocate enough memory for all the forces
  int i = 0;
  while(fscanf(input_file, "%s\n", &buffer) != EOF) //second pass: read the forces
    target_forces[i++] = force->numeric(FLERR,buffer);
  fclose(input_file);

  read_tersoff_bounds_file();

  dynamic_group_allow = 1; //probably not needed for this fix
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixParameterize::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixParameterize::init()
{
  random = new RanMars(lmp, random_seed + comm->me);

  tersoff = (PairTersoff*) force->pair_match("tersoff",1,0);
  if(tersoff==NULL) error->all(FLERR,"Can't find Tersoff pair style in this run");
  
  printf("tersoff nparams = %d\n", tersoff->nparams);
  for(int i=0; i<tersoff->nparams; i++) {
    printf("tersoff elements = %d %d %d, A = %f\n", tersoff->params[i].ielement, tersoff->params[i].jelement, tersoff->params[i].kelement, tersoff->params[i].biga);
  }
  //copy starting Tersoff parameters
  best_tersoff_params = new PairTersoff::Param[tersoff->nparams];
  memcpy(best_tersoff_params, tersoff->params, tersoff->nparams*sizeof(PairTersoff::Param));
  //print starting Tersoff parameters
  for(int i=0; i<tersoff->nparams; i++) {
    printf("Tersoff param A = %f, %f\n", tersoff->params[i].biga, best_tersoff_params[i].biga);
  }
  //copy starting charges and LJ parameters by atom type
  best_charges = new double[atom->ntypes];
  best_lj_sigma = new double[atom->ntypes];
  best_lj_epsilon = new double[atom->ntypes];
  
  lj = (PairLJCutCoulInOut*) force->pair_match("lj/cut/coul/inout",1,0);
  if(lj==NULL) error->all(FLERR,"Can't find lj/cut/coul/inout pair style in this run");
  
  for(int type=0; type<atom->ntypes; type++) {
    //charges belong to atoms, so we must find an atom with type i. LAMMPS types start at 1, not 0. 
    int atom_of_type = 0;
    while(atom_of_type < atom->natoms) {
      if( atom->type[atom_of_type] == type+1 )
        break;
      atom_of_type++;
    }
    best_charges[type] = atom->q[atom_of_type];
    //get LJ params direct from pair style
    best_lj_sigma[type] = lj->sigma[type+1][type+1];
    best_lj_epsilon[type] = lj->epsilon[type+1][type+1];
    
    printf("Type %d: q=%-5.2f, r=%-5.2f, e=%-5.2f\n", type, best_charges[type], best_lj_sigma[type], best_lj_epsilon[type]);
  }
  printf("Atom types: %d\n", atom->ntypes);
  
  //get ready to record results
  best_error = -1;
  counter_since_last_file_write = 0;
  ready_to_write_file = 0;
}

double FixParameterize::calculate_error()
{
  double **f = atom->f;
  //compare target forces with force vector **f
  double squared_error = 0.0;
  for(int i=0; i<n_target_forces; i+=3) {
    double dfx = f[i/3][0] - target_forces[i];
    double dfy = f[i/3][1] - target_forces[i+1];
    double dfz = f[i/3][2] - target_forces[i+2];
    squared_error += dfx*dfx + dfy*dfy + dfz*dfz;
  }
  return squared_error;
}

void FixParameterize::initial_integrate(int vflag) //modify parameters before the step
{
  tersoff->cutmax = 0.0; //changes if R or D are changed
  
  for(int i=0; i<tersoff->nparams; i++) {
    if(tersoff->params[i].ielement==0 && tersoff->params[i].jelement==1 && tersoff->params[i].kelement==1) {
      
      tersoff->params[i].biga = best_tersoff_params[i].biga * (1.0 + random->gaussian()*0.01);
      tersoff->params[i].bigb = best_tersoff_params[i].bigb * (1.0 + random->gaussian()*0.01);
      
      //recalculate resulting parameters: cut, cutsq, c1, c2, c3, c4, and cutmax
      tersoff->params[i].cut = tersoff->params[i].bigr + tersoff->params[i].bigd;
      tersoff->params[i].cutsq = tersoff->params[i].cut*tersoff->params[i].cut;
      tersoff->params[i].c1 = pow(2.0*tersoff->params[i].powern*1.0e-16,-1.0/tersoff->params[i].powern);
      tersoff->params[i].c2 = pow(2.0*tersoff->params[i].powern*1.0e-8,-1.0/tersoff->params[i].powern);
      tersoff->params[i].c3 = 1.0/tersoff->params[i].c2;
      tersoff->params[i].c4 = 1.0/tersoff->params[i].c1;
      if(tersoff->cutmax < tersoff->params[i].cut)
        tersoff->cutmax = tersoff->params[i].cut;
    }
  }
  //modify charges
  for(int a=0; a < atom->natoms; a++) {
    for(int type=0; type < atom->ntypes; type++) {
      if(atom->type[a] == type+1)
        atom->q[a] = best_charges[type] * (1.0 + random->gaussian()*0.01);
        //todo: make sure there are no resulting parameters based on charge that need to be recalculated after changing charge
        //todo: check whether atoms are indexed from 1 or 0
    }
  }
  //modify LJ parameters
  for(int type=0; type < atom->ntypes; type++) {
    lj->sigma[type+1][type+1] = best_lj_sigma[type] * (1.0 + random->gaussian()*0.01);
    lj->epsilon[type+1][type+1] = best_lj_epsilon[type] * (1.0 + random->gaussian()*0.01);
  }
  //recalculate resulting parameters
  for(int i=1; i <= atom->ntypes; i++)
    for(int j=1; j <= atom->ntypes; j++)
      lj->init_one(i, j);
}

void FixParameterize::final_integrate() //check the results after the step
{
  if(best_error<0) best_error = calculate_error();
  
  double new_error = calculate_error();
  
  if(new_error < best_error) {
    memcpy(best_tersoff_params, tersoff->params, tersoff->nparams*sizeof(PairTersoff::Param));
    best_error = new_error;
    ready_to_write_file = 1;
  }
  counter_since_last_file_write++;
  if(ready_to_write_file && counter_since_last_file_write>1000)
      write_tersoff_file();
  
  printf("%f %f\n", best_error, new_error);
}

/*
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4;
    int ielement,jelement,kelement;
    int powermint;
    double Z_i,Z_j;              // added for TersoffZBL
    double ZBLcut,ZBLexpscale;
    double c5,ca1,ca4;           // added for TersoffMOD
    double powern_del;
  };
*/

void FixParameterize::write_tersoff_file() {
  //open file
  FILE *output = fopen("best.tersoff", "w");
  //output error
  fprintf(output, "# Error: %g\n", best_error);
  //output charges on one line, listed by type
  fprintf(output, "# Charges: ");
  for(int type=0; type < atom->ntypes; type++) fprintf(output, "%9.6g", best_charges[type]);
  fprintf(output, "\n");
  //output LJ sigma on one line, listed by type
  fprintf(output, "# LJ-sigma: ");
  for(int type=0; type < atom->ntypes; type++) fprintf(output, "%9.6g", best_lj_sigma[type]);
  fprintf(output, "\n");
  //output LJ epsilon on one line, listed by type
  fprintf(output, "# LJ-epsilon: ");
  for(int type=0; type < atom->ntypes; type++) fprintf(output, "%9.6g", best_lj_epsilon[type]);
  fprintf(output, "\n");
  //output Tersoff parameters by pair
  for(int i=0; i<tersoff->nparams; i++) {
    PairTersoff::Param t = best_tersoff_params[i];
    fprintf(output, "%-2s %-2s %-2s ", tersoff->elements[t.ielement], tersoff->elements[t.jelement], tersoff->elements[t.kelement]);
    fprintf(output, "%9.6g %9.6g %9.6g %9.6g %9.6g %9.6g\n", t.powerm, t.gamma, t.lam3, t.c, t.d, t.h);
    fprintf(output, "         %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g\n\n", t.powern, t.beta, t.lam2, t.bigb, t.bigr, t.bigd, t.lam1, t.biga);
  }
  fclose(output);
  counter_since_last_file_write = 0;
  ready_to_write_file = 0;
}

void FixParameterize::read_tersoff_bounds_file() {
  char **args = new char*[5];
  args[0] = "*";
  args[1] = "*";
  args[2] = "bounds.tersoff";
  args[3] = "Pb";
  args[4] = "I";
  
  PairTersoff *tersoff_upper_bound = new PairTersoff(lmp);
  tersoff_upper_bound->coeff(5, args);
  tersoff_upper_bound->read_file((char *)"bounds.tersoff");
  
  PairTersoff *tersoff_lower_bound = new PairTersoff(lmp);
  tersoff_lower_bound->coeff(5, args);
  tersoff_lower_bound->read_file((char *)"bounds.tersoff");
  
  //exit(0);
}



