/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   This file is a (heavy) modification of min_fire.cpp
------------------------------------------------------------------------- */

#include <math.h>
#include "min_params.h"
#include "universe.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "output.h"
#include "timer.h"
#include "error.h"

#include "random_mars.h" //random number generator, Marsaglia style
#include "pair_tersoff.h"
#include "pair_lj_cut_coul_inout.h"
#include "compute.h" //for compute pe/atom

#include <vector> //for std:vector
#include <algorithm> //for std:copy and std::fill
#include <fstream> //for file reading
#include <string> //for file reading
#include <string.h> //for strtok

#include <assert.h>

#include <nlopt.h> //NLopt Non-Linear Optimization library. Compilation flags: gcc -I/fs/home/jms875/install/include -L/fs/home/jms875/install/lib -lnlopt -lm

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

/* ---------------------------------------------------------------------- */

MinParams::MinParams(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */

void MinParams::modify_params(int narg, char **arg)
{
  if (narg != 3)
    error->all(FLERR,"Illegal min_modify command: min_style 'params' takes three arguments (run_name, optimization method (0-2), and random_seed)");

  run_name = std::string(arg[0]); //prefix for all output files
  optimization_method = force->inumeric(FLERR,arg[1]); //int 0-2, eg RANDOM=1
  random_seed = force->inumeric(FLERR,arg[2]);
  
  printf("Run name = %s\n", run_name.c_str());
  printf("Optimization method = %d\n", optimization_method);
  printf("Random seed = %d\n", random_seed);
}

void MinParams::init()
{
  Min::init();
  
  if(run_name.size()==0)
    error->all(FLERR,"With min_style params, must call min_modify with two arguments (run_name and random_seed) before starting minimize");
  std::string input_forces_filename = run_name+"_forces.txt";
  std::string input_energies_filename = run_name+"_energies.txt";
  std::string upper_bounds_filename = run_name+"_upper.tersoff";
  std::string lower_bounds_filename = run_name+"_lower.tersoff";
  std::string input_tersoff_filename = run_name+"_input.tersoff";
  
  //read file containing target energies (expecting to be in a flat list, one number on each line)
  FILE *input_file = fopen(input_energies_filename.c_str(), "r");
  if(input_file==NULL) error->all(FLERR,"No such target-energies file exists");
  char buffer[100];
  while(fscanf(input_file, "%s\n", (char*)&buffer) != EOF) //read the forces
    target_energies.push_back( force->numeric(FLERR,buffer) );
  fclose(input_file);
  current_energies = target_energies;
  
  //read file containing target forces (expecting to be in a flat list, one number on each line)
  input_file = fopen(input_forces_filename.c_str(), "r");
  if(input_file==NULL) error->all(FLERR,"No such target-forces file exists");
  while(fscanf(input_file, "%s\n", (char*)&buffer) != EOF) //read the forces
    target_forces.push_back( force->numeric(FLERR,buffer) );
  fclose(input_file);

  //read Tersoff bounds
  tersoff_upper_bound = new PairTersoff(lmp);
  tersoff_lower_bound = new PairTersoff(lmp);
  const char *tersoff_bounds_args[] = {"*", "*", upper_bounds_filename.c_str(), "Pb", "Cl"}; //todo: must use real elements
  tersoff_upper_bound->coeff(5, (char**)tersoff_bounds_args);
  tersoff_bounds_args[2] = lower_bounds_filename.c_str();
  tersoff_lower_bound->coeff(5, (char**)tersoff_bounds_args);
  
  //read other upper bounds
  read_params_from_comments(upper_bounds_filename, charges_upper, lj_sigma_upper, lj_epsilon_upper);
  read_params_from_comments(lower_bounds_filename, charges_lower, lj_sigma_lower, lj_epsilon_lower);
  read_params_from_comments(input_tersoff_filename, charges_current, lj_sigma_current, lj_epsilon_current); //todo: change "input.tersoff" to a user-input variable
  
  random = new RanMars(lmp, random_seed + comm->me);

  tersoff = (PairTersoff*) force->pair_match("tersoff",1,0);
  if(tersoff==NULL) error->all(FLERR,"Can't find Tersoff pair style in this run");
  
  int id = modify->find_compute("atom_pe"); //must use this name for the compute!
  if (id < 0) error->all(FLERR,"Minimization could not find pe/atom compute (must be named atom_pe)");
  compute_pe = modify->compute[id];
  
  printf("tersoff nparams = %d\n", tersoff->nparams);
  for(int i=0; i<tersoff->nparams; i++) {
    printf("tersoff elements = %d %d %d, A = %f\n", tersoff->params[i].ielement, tersoff->params[i].jelement, tersoff->params[i].kelement, tersoff->params[i].biga);
  }
  for(int type=0; type < atom->ntypes; type++) { //for testing only
	  printf("Pairwise: type %d, e=%f, s=%f, q=%f\n", type, lj_epsilon_current[type], lj_sigma_current[type], charges_current[type]);
  }
  
  //get LAMMPS lj pair object
  lj = (PairLJCutCoulInOut*) force->pair_match("lj/cut/coul/inout",1,0);
  if(lj==NULL) error->all(FLERR,"Can't find lj/cut/coul/inout pair style in this run");
  
  pack_params(); //read parameters from LAMMPS objects into flat arrays
  params_best = params_current;
  unpack_params(params_current); //put parameters into LAMMPS objects so as to calculate forces next step
  
  printf("Parameters to optimize: %lu\n", params_current.size());
  
  //get ready to record results
  best_error = -1;
  counter_since_last_file_write = 0;
  ready_to_write_file = 0;
}

/* ---------------------------------------------------------------------- */

void MinParams::setup_style()
{
  double **v = atom->v;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    v[i][0] = v[i][1] = v[i][2] = 0.0;
}

/* ----------------------------------------------------------------------
   set current vector lengths and pointers
   called after atoms have migrated
------------------------------------------------------------------------- */

void MinParams::reset_vectors()
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

/* ---------------------------------------------------------------------- */

int MinParams::iterate(int maxiter)
{
  bigint ntimestep;

  for (int iter = 0; iter<maxiter; iter++) {

    ntimestep = ++update->ntimestep;
    niter++;

    ecurrent = energy_force(0);
    neval++;
    if(neval>update->max_eval) {
      return MAXEVAL;
    }


    update->eflag_atom = update->ntimestep; //needed to avoid error message. TODO: make this less of a hack. 
    compute_pe->compute_peratom(); //also sort of a hack



    if(best_error<0) {
      best_error = calculate_error();
      write_tersoff_file();
    }
    double new_error = calculate_error();

    if(0) { //for testing
      for(unsigned int i=0; i<params_current.size(); i++) {
        printf("%.3g->%.3g ", params_best[i], params_current[i]);
      }
      puts("");
      printf("Error: F=%g  E=%g  Sum=%g\n", force_error, energy_error, best_error);
    }
    
    if(new_error < best_error) {
      //copy params_current into params_best
      std::copy ( params_current.begin(), params_current.end(), params_best.begin() );
      best_error = new_error;
      ready_to_write_file = 1;
      printf("New best: %f\n", best_error );
    }
    counter_since_last_file_write++;
    if(ready_to_write_file && counter_since_last_file_write>1000)
        write_tersoff_file();
        
    //modify params_current to make a new guess for the next step
    
    int RANDOM = 0;
    int DDS = 1;
    int HILL_CLIMB = 2;

    double step_size = 0.2; //default in DDS = 0.2
    
    if(optimization_method==RANDOM) { //Pick a random point within the bounds
      for(unsigned int i=0; i<params_current.size(); i++) {
        params_current[i] = params_lower[i] + (params_upper[i]-params_lower[i])*random->uniform();
      }
    }
    if(optimization_method==DDS) { //Dynamically Dimensioned Search
      for(unsigned int i=0; i<params_current.size(); i++) {
        double fraction_done = 1.0*update->ntimestep/update->nsteps;
        if( random->uniform() < fraction_done-1.0/params_current.size() ) continue; //skip this param (Dynamically Dimensioned Search)
        double dx = (params_upper[i] - params_lower[i]) * step_size * random->gaussian();
        params_current[i] = params_best[i] + dx;
        if(params_current[i] > params_upper[i]) params_current[i] = 2*params_upper[i] - params_current[i]; //reflect across bound
        if(params_current[i] < params_lower[i]) params_current[i] = 2*params_lower[i] - params_current[i]; //reflect across bound
      }
    }
    else if(optimization_method==HILL_CLIMB) { //Change just one param at a time
      unsigned int i = int(random->uniform()*params_current.size()); 
      double dx = (params_upper[i] - params_lower[i]) * step_size * random->gaussian();
      params_current[i] = params_best[i] + dx;
      if(params_current[i] > params_upper[i]) params_current[i] = 2*params_upper[i] - params_current[i]; //reflect across bound
      if(params_current[i] < params_lower[i]) params_current[i] = 2*params_lower[i] - params_current[i]; //reflect across bound
    }
    
    //put parameters into LAMMPS objects so as to calculate forces next step
    unpack_params(params_current);

    //convergence criteria
    if(energy_error < update->etol) {
       return ETOL;
    }
    if(force_error<update->ftol) {
      return FTOL;
    }

    // output for thermo, dump, restart files
    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }

  return MAXITER;
}


void MinParams::read_params_from_comments(std::string filename, std::vector<double> &charges, std::vector<double> &lj_sigma, std::vector<double> &lj_epsilon) {
  charges.assign(atom->ntypes, 0.0);
  lj_sigma.assign(atom->ntypes, 0.0);
  lj_epsilon.assign(atom->ntypes, 0.0);
  
  std::ifstream infile(filename.c_str());
  std::string line;
  std::string charge_label = "# Charges:";
  std::string lj_sigma_label = "# LJ-sigma:";
  std::string lj_epsilon_label = "# LJ-epsilon:";
  while(std::getline(infile, line)) {
    if(line.rfind(charge_label, 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
        char *token = strtok( (type==0?((char*)line.c_str()+charge_label.size()):NULL), " ");
        charges[type] = force->numeric(FLERR,token); //strtok needs first arg to be NULL on all calls after the first
      }
    }
    if(line.rfind(lj_sigma_label, 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
        char *token = strtok( (type==0?((char*)line.c_str()+lj_sigma_label.size()):NULL), " ");
        lj_sigma[type] = force->numeric(FLERR,token); //strtok needs first arg to be NULL on all calls after the first
      }
    }
    if(line.rfind(lj_epsilon_label, 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
        char *token = strtok( (type==0?((char*)line.c_str()+lj_epsilon_label.size()):NULL), " ");
        lj_epsilon[type] = force->numeric(FLERR,token); //strtok needs first arg to be NULL on all calls after the first
      }
    }
  }
}


double MinParams::calculate_error()
{
  if(compute_pe->vector_atom==NULL) {
    error->all(FLERR,"compute_pe->vector_atom is NULL. Energies have not been computed!");
  }

  double **f = atom->f;
  double **x = atom->x;
  //compare target forces with force vector **f
  force_error = 0.0;
  for(unsigned int i=0; i<target_forces.size(); i+=3) {
    double dfx = f[i/3][0] - target_forces[i];
    double dfy = f[i/3][1] - target_forces[i+1];
    double dfz = f[i/3][2] - target_forces[i+2];
    double target_magnitude = target_forces[i]*target_forces[i] + target_forces[i+1]*target_forces[i+1] + target_forces[i+2]*target_forces[i+2];
    double small_magnitude = 1e-4;
    force_error += (dfx*dfx + dfy*dfy + dfz*dfz) / (target_magnitude+small_magnitude);
  }
  //get current energies
  std::fill(current_energies.begin(), current_energies.end(), 0.0);
  for(unsigned int i=0; i<atom->natoms; i++) {
    unsigned int which_system = int((x[i][0]+100.0)/1000.0);
    current_energies[ which_system ] += compute_pe->vector_atom[i];
  }
  //set baseline for each type of system, denoted by E=0.0 in energy input file
  double baseline_energy = 0.0;
  for(unsigned int i=0; i<target_energies.size(); i++) {
    if(target_energies[i]==0.0)
      baseline_energy = current_energies[i];
    current_energies[i] -= baseline_energy;
  }
  //compare target energies
  double kT = 0.6; //kcal/mol
  energy_error = 0.0;
  for(unsigned int i=0; i<target_energies.size(); i++) {
    double error = (target_energies[i]-current_energies[i])/(target_energies[i]+kT);
    energy_error += error*error;
  }
  //for testing
  if(0) {
    for(unsigned int i=0; i<target_energies.size(); i++) {
      printf("E: %g %g\n", target_energies[i], current_energies[i]);
    }
    puts("Done printing energies");
    exit(0);
  }
  
  force_error = sqrt(force_error/atom->natoms); //Normalize error
  energy_error = sqrt(energy_error/target_energies.size()); //Normalize error
  return 0.0*force_error + energy_error;
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


void MinParams::write_tersoff_file() {
  //get parameters from flat array back into LAMMPS objects
  unpack_params(params_best);
  //open file
  std::string output_filename = run_name+"_best.tersoff";
  FILE *output = fopen(output_filename.c_str(), "w"); //todo: make this filename a variable
  //label
  fprintf(output, "# Generated by fix_parameterize.cpp. LAMMPS units = real.\n");
  //output error
  fprintf(output, "# Error metric=%g, Energy error=%.2f%%, Force error=%.2f%%\n", best_error, energy_error*100, force_error*100);
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

void MinParams::pack_params() {
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
  
  //pack other parameters into a flat array
  for(int type=0; type < atom->ntypes; type++) {
    pack_typewise_param(charges_upper, charges_lower, charges_current);
    pack_typewise_param(lj_sigma_upper, lj_sigma_lower, lj_sigma_current);
    pack_typewise_param(lj_epsilon_upper, lj_epsilon_lower, lj_epsilon_current);
  }
  
}

#define unpack_param(X,I) if( tersoff_upper_bound->params[i].X > tersoff_lower_bound->params[i].X ) { tersoff->params[i].X = pp[I]; I++; }

#define unpack_typewise_param(HI,LO,X,I) if( HI[type]>LO[type] ) { X[type] = pp[I]; I++; }

void MinParams::unpack_params(std::vector<double> pp) {
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
    //printf("my cutsq %d %d %d = %f\n", tersoff->params[i].ielement, tersoff->params[i].jelement, tersoff->params[i].kelement, tersoff->params[i].cutsq);
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
      if(atom->type[a] == type+1) //atom->type[] is zero-indexed
        atom->q[a] = charges_current[type]; //atom->q[] is zero-indexed
        //todo: make sure there are no resulting parameters based on charge that need to be recalculated after changing charge
    }
	//printf("Atom %d: t=%d, q=%f\n", a, atom->type[a], atom->q[a]);
  }
  
  //modify LJ parameters
  for(int type=0; type < atom->ntypes; type++) {
    lj->sigma[type+1][type+1] = lj_sigma_current[type]; //lj->sigma and lj->epsilon are 1-indexed
    lj->epsilon[type+1][type+1] = lj_epsilon_current[type];
  }
  //recalculate resulting parameters
  for(int i=1; i <= atom->ntypes; i++) {
    for(int j=i; j <= atom->ntypes; j++) {
      lj->setflag[i][j] = 0; //unset
      lj->init_one(i, j);
    }
  }
  /*for(int i=1; i <= atom->ntypes; i++) { //for testing only
    for(int j=i; j <= atom->ntypes; j++) {
      printf("Pair %d %d: e=%f s=%f\n", i, j, lj->epsilon[i][j], lj->sigma[i][j]);
    }
  }*/
}

double MinParams::NLopt_target_function(unsigned params_count, const double *params, double *gradient, void *optional_data) {
  return 0.0;
}

