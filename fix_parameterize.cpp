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

#include <vector> //for std:vector
#include <algorithm> //for std:copy

#include <fstream> //for file reading
#include <string> //for file reading

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

void FixParameterize::read_params_from_comments(const char *filename, std::vector<double> &charges, std::vector<double> &lj_sigma, std::vector<double> &lj_epsilon) {
  charges.assign(atom->ntypes, 0.0);
  lj_sigma.assign(atom->ntypes, 0.0);
  lj_epsilon.assign(atom->ntypes, 0.0);
  
  std::ifstream infile(filename);
  std::string line;
  while(std::getline(infile, line)) {
    if(line.rfind("# Charges:", 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
          charges[type] = force->numeric(FLERR,strtok( (type==0?((char*)line.c_str()+10):NULL), " ")); //strtok needs first arg to be NULL on all calls after the first
      }
    }
    if(line.rfind("# LJ-sigma:", 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
          lj_sigma[type] = force->numeric(FLERR,strtok( (type==0?((char*)line.c_str()+11):NULL), " ")); //strtok needs first arg to be NULL on all calls after the first
      }
    }
    if(line.rfind("# LJ-epsilon:", 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
          lj_epsilon[type] = force->numeric(FLERR,strtok( (type==0?((char*)line.c_str()+13):NULL), " ")); //strtok needs first arg to be NULL on all calls after the first
      }
    }
  }
}

FixParameterize::FixParameterize(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8)
    error->all(FLERR,"Illegal fix parameterize command: requires input file (target forces), upper bounds file, lower bounds file, output file, and a random seed 1-2^31");
  //get inputs from user
  char* input_forces_filename = arg[3];
  char* upper_bounds_filename = arg[4];
  char* lower_bounds_filename = arg[5];
  output_filename = arg[6];
  random_seed = force->inumeric(FLERR,arg[7]);
  
  //read file containing target forces (expecting to be in a flat list, one number on each line)
  FILE *input_file = fopen(input_forces_filename, "r");
  if(input_file==NULL) error->all(FLERR,"No such target-forces file exists");
  char buffer[100];
  while(fscanf(input_file, "%s\n", (char*)&buffer) != EOF) //read the forces
    target_forces.push_back( force->numeric(FLERR,buffer) );
  fclose(input_file);

  //read Tersoff bounds
  tersoff_upper_bound = new PairTersoff(lmp);
  tersoff_lower_bound = new PairTersoff(lmp);
  const char *tersoff_bounds_args[] = {"*", "*", upper_bounds_filename, "Pb", "Cl"}; //todo: must use real elements
  tersoff_upper_bound->coeff(5, (char**)tersoff_bounds_args);
  tersoff_bounds_args[2] = lower_bounds_filename;
  tersoff_lower_bound->coeff(5, (char**)tersoff_bounds_args);
  
  //read other upper bounds
  read_params_from_comments(upper_bounds_filename, charges_upper, lj_sigma_upper, lj_epsilon_upper);
  read_params_from_comments(lower_bounds_filename, charges_lower, lj_sigma_lower, lj_epsilon_lower);
  read_params_from_comments((char*) "input.tersoff", charges_current, lj_sigma_current, lj_epsilon_current); //todo: change "input.tersoff" to a user-input variable
  
  //variables for LAMMPS fix
  dynamic_group_allow = 1; //probably not needed for this fix
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixParameterize::setmask()
{
  int mask = 0;
  mask |= FINAL_INTEGRATE; //say that we need the final_integrate() function
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
  //copy starting charges and LJ parameters by atom type  
  lj = (PairLJCutCoulInOut*) force->pair_match("lj/cut/coul/inout",1,0);
  if(lj==NULL) error->all(FLERR,"Can't find lj/cut/coul/inout pair style in this run");
  
  pack_params();
  params_best = params_current;
  unpack_params(params_current);
  
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
  for(unsigned int i=0; i<target_forces.size(); i+=3) {
    double dfx = f[i/3][0] - target_forces[i];
    double dfy = f[i/3][1] - target_forces[i+1];
    double dfz = f[i/3][2] - target_forces[i+2];
    double target_magnitude = target_forces[i]*target_forces[i] + target_forces[i+1]*target_forces[i+1] + target_forces[i+2]*target_forces[i+2];
    double small_magnitude = 1e-4;
    squared_error += (dfx*dfx + dfy*dfy + dfz*dfz) / (target_magnitude+small_magnitude);
    //printf("%e\n%e\n%e\n", f[i/3][0], f[i/3][1], f[i/3][2]);
  }
  //exit(0);
  return squared_error;
}

void FixParameterize::final_integrate() //check the results after the step
{
  if(best_error<0) {
    best_error = calculate_error();
    write_tersoff_file();
  }
  
  double new_error = calculate_error();
  
  //printf("%.1e    %.1e\n", best_error, new_error);
  //printf("%d  %d\n", update->nsteps, update->ntimestep);
  
  if(new_error < best_error) {
    //copy params_current into params_best
    std::copy ( params_current.begin(), params_current.end(), params_best.begin() );
    best_error = new_error;
    ready_to_write_file = 1;
    printf("New best: %f\n", sqrt(best_error/atom->natoms) );
  }
  counter_since_last_file_write++;
  if(ready_to_write_file && counter_since_last_file_write>1000)
      write_tersoff_file();
      
  //modify params_current to make a new guess for the next step
  
  int RANDOM = 0;
  int DDS = 1;
  int method = DDS;
  
  if(method==DDS) { //Dynamically Dimensioned Search
    for(unsigned int i=0; i<params_current.size(); i++) {
      double fraction_done = 1.0*update->ntimestep/update->nsteps;
      if( random->uniform() < fraction_done-1.0/params_current.size() ) continue; //skip this param (Dynamically Dimensioned Search)
      double dx = (params_upper[i] - params_lower[i]) * 0.2 * random->gaussian();
      params_current[i] = params_best[i] + dx;
      if(params_current[i] > params_upper[i]) params_current[i] = 2*params_upper[i] - params_current[i]; //reflect across bound
      if(params_current[i] < params_lower[i]) params_current[i] = 2*params_lower[i] - params_current[i]; //reflect across bound
    }
  }
  else if(method==RANDOM) { //Pick a random point within the bounds
    for(unsigned int i=0; i<params_current.size(); i++) {
      params_current[i] = params_lower[i] + (params_upper[i]-params_lower[i])*random->uniform();
    }
  }
  
  //put parameters into LAMMPS objects so as to calculate forces next step
  unpack_params(params_current);
}

/*
  struct Param {
    double lam1,lam2,lam3;
    double c,d,h;
    double gamma,powerm;
    double powern,beta;
    double biga,bigb,bigd,bigr;
    double cut,cutsq;
    double c1,c2,c3,c4; //derived params
    int ielement,jelement,kelement;
    int powermint; //derived params
    double Z_i,Z_j;              // added for TersoffZBL
    double ZBLcut,ZBLexpscale;
    double c5,ca1,ca4;           // added for TersoffMOD
    double powern_del;
  };
*/

void FixParameterize::write_tersoff_file() {
  //get parameters from flat array back into LAMMPS objects
  unpack_params(params_best);
  //open file
  FILE *output = fopen(output_filename, "w"); //todo: make this filename a variable
  //label
  fprintf(output, "# Generated by fix_parameterize.cpp. LAMMPS units = real.\n");
  //output error
  fprintf(output, "# Error: %g\n", sqrt(best_error/atom->natoms) );
  //output charges on one line, listed by type
  fprintf(output, "# Charges:    ");
  for(int type=0; type < atom->ntypes; type++) fprintf(output, " %9.6g ", charges_current[type]);
  fprintf(output, "\n");
  //output LJ sigma on one line, listed by type
  fprintf(output, "# LJ-sigma:   ");
  for(int type=0; type < atom->ntypes; type++) fprintf(output, " %9.6g ", lj_sigma_current[type]);
  fprintf(output, "\n");
  //output LJ epsilon on one line, listed by type
  fprintf(output, "# LJ-epsilon: ");
  for(int type=0; type < atom->ntypes; type++) fprintf(output, " %9.6g ", lj_epsilon_current[type]);
  fprintf(output, "\n");
  //output legend
  fprintf(output, "# i, j, k,       m,        gamma,    lambda3,  c,        d,        costheta0,\n#                n,        beta,     lambda2,  B,        R,        D,        lambda1,  A\n");
  //output Tersoff parameters by pair
  for(int i=0; i<tersoff->nparams; i++) {
    PairTersoff::Param t = tersoff->params[i];
    fprintf(output, "%-2s %-2s %-2s ", tersoff->elements[t.ielement], tersoff->elements[t.jelement], tersoff->elements[t.kelement]);
    fprintf(output, "%9.6g %9.6g %9.6g %9.6g %9.6g %9.6g\n", t.powerm, t.gamma, t.lam3, t.c, t.d, t.h);
    fprintf(output, "         %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g\n\n", t.powern, t.beta, t.lam2, t.bigb, t.bigr, t.bigd, t.lam1, t.biga);
  }
  fclose(output);
  counter_since_last_file_write = 0;
  ready_to_write_file = 0;
}

//Macro to reduce the repetitive code below. This is needed because C++ lacks "reflection" - there is no way to loop through a list of an object's members
#define pack_param(X) if( hi.X > lo.X ) { params_upper.push_back(hi.X); params_lower.push_back(lo.X); params_current.push_back(tt.X); }

#define pack_typewise_param(HI,LO,X) if(HI[type] > LO[type]) { params_upper.push_back(HI[type]); params_lower.push_back(LO[type]); params_current.push_back(X[type]); }

void FixParameterize::pack_params() {
  //pack Tersoff parameters into a flat array
  for(int i=0; i<tersoff->nparams; i++) {
    PairTersoff::Param hi = tersoff_upper_bound->params[i];
    PairTersoff::Param lo = tersoff_lower_bound->params[i];
    PairTersoff::Param tt = tersoff->params[i];
    pack_param(lam1);
    pack_param(lam2);
    pack_param(lam3);
    pack_param(c);
    pack_param(d);
    pack_param(h);
    pack_param(gamma);
    pack_param(powerm);
    pack_param(powern);
    pack_param(beta);
    pack_param(biga);
    pack_param(bigb);
    pack_param(bigd);
    pack_param(bigr);
    pack_param(cut);
  }
  /*
  //get other parameters
  for(int type=0; type < atom->ntypes; type++) {
    //charges belong to atoms, so we must find an atom with type i. LAMMPS types start at 1, not 0. 
    int atom_of_type = 1;
    while(atom_of_type < atom->natoms) {
      if( atom->type[atom_of_type] == type+1 )
        break;
      atom_of_type++;
    }
    charges_current.push_back( atom->q[atom_of_type] );
    //get LJ params direct from pair style
    lj_sigma_current.push_back( lj->sigma[type+1][type+1] );
    lj_epsilon_current.push_back( lj->epsilon[type+1][type+1] );
  }
  */
  //pack other parameters into a flat array
  for(int type=0; type < atom->ntypes; type++) {
    pack_typewise_param(charges_upper, charges_lower, charges_current);
    pack_typewise_param(lj_sigma_upper, lj_sigma_lower, lj_sigma_current);
    pack_typewise_param(lj_epsilon_upper, lj_epsilon_lower, lj_epsilon_current);
  }
}

#define unpack_param(X,I) if( tersoff_upper_bound->params[i].X > tersoff_lower_bound->params[i].X ) { tersoff->params[i].X = pp[I]; I++; }

#define unpack_typewise_param(HI,LO,X,I) if( HI[type]>LO[type] ) { X[type] = pp[I]; I++; }

void FixParameterize::unpack_params(std::vector<double> pp) {
  int which = 0;
  //unpack Tersoff parameters from the flat array to the LAMMPS object
  for(int i=0; i<tersoff->nparams; i++) {
    unpack_param(lam1, which);
    unpack_param(lam2, which);
    unpack_param(lam3, which);
    unpack_param(c, which);
    unpack_param(d, which);
    unpack_param(h, which);
    unpack_param(gamma, which);
    unpack_param(powerm, which);
    unpack_param(powern, which);
    unpack_param(beta, which);
    unpack_param(biga, which);
    unpack_param(bigb, which);
    unpack_param(bigd, which);
    unpack_param(bigr, which);
    unpack_param(cut, which);
    
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
  //unpack other parameters from the flat array to the LAMMPS object
  for(int type=0; type < atom->ntypes; type++) {
    unpack_typewise_param(charges_upper, charges_lower, charges_current, which);
    unpack_typewise_param(lj_sigma_upper, lj_sigma_lower, lj_sigma_current, which);
    unpack_typewise_param(lj_epsilon_upper, lj_epsilon_lower, lj_epsilon_current, which);
  }
  //modify charges
  for(int a=0; a < atom->natoms; a++) {
    for(int type=0; type < atom->ntypes; type++) {
      if(atom->type[a+1] == type+1)
        atom->q[a+1] = charges_current[type];
        //todo: make sure there are no resulting parameters based on charge that need to be recalculated after changing charge
    }
  }
  //modify LJ parameters
  for(int type=0; type < atom->ntypes; type++) {
    lj->sigma[type+1][type+1] = lj_sigma_current[type];
    lj->epsilon[type+1][type+1] = lj_epsilon_current[type];
  }
  //recalculate resulting parameters
  for(int i=1; i <= atom->ntypes; i++)
    for(int j=1; j <= atom->ntypes; j++)
      lj->init_one(i, j);
}
