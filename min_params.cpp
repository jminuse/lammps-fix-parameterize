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
#include <string> //for C++ string syntax
#include <string.h> //for strtok() function

#include <assert.h>

#include <nlopt.h> //NLopt Non-Linear Optimization library. Compilation flags: gcc -I/fs/home/jms875/install/include -L/fs/home/jms875/install/lib -lnlopt -lm

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

/* ---------------------------------------------------------------------- */

MinParams::MinParams(LAMMPS *lmp) : Min(lmp) {}

#define OPT_PRINT_ALL NLOPT_LD_LBFGS //used as code for "PRINT"
#define OPT_RAND_SBPLX NLOPT_LD_SLSQP //used as code for RAND_SBPLX

/* ---------------------------------------------------------------------- */

void MinParams::modify_params(int narg, char **arg)
{
  if (narg != 3)
    error->all(FLERR,"Illegal min_modify command: min_style 'params' takes three arguments (run_name, optimization method (default=SBPLX), and random_seed)");

  run_name = std::string(arg[0]); //prefix for all output files
  std::string optimization_method = std::string(arg[1]); //name of algorithm
  random_seed = force->inumeric(FLERR,arg[2]);
  
  if(optimization_method == "DIRECT") { // DIviding RECTangles algorithm, deterministic
    algorithm = NLOPT_GN_DIRECT;
  } else if(optimization_method == "DIRECT-L") { // DIviding RECTangles algorithm, Locally-biased, non-deterministic
    algorithm = NLOPT_GN_DIRECT_L_RAND;
  } else if(optimization_method == "CRS") { // Controlled Random Search
    algorithm = NLOPT_GN_CRS2_LM;
  } else if(optimization_method == "ISRES") { // Improved Stochastic Ranking Evolution Strategy. Mutation + local simplex optimization. 
    algorithm = NLOPT_GN_ISRES;
  } else if(optimization_method == "COBYLA") { // "Error -4: Halted because roundoff errors limited progress."
    algorithm = NLOPT_LN_COBYLA;
  } else if(optimization_method == "ESCH") { // evolutionary algorithm. Certainly slow, possibly steady. 
    algorithm = NLOPT_GN_ESCH;
  } else if(optimization_method == "BOBYQA") { // "Error -4: Halted because roundoff errors limited progress."
    algorithm = NLOPT_LN_BOBYQA;
  } else if(optimization_method == "NELDERMEAD") { // Nelder-Mead simplex. "Error -4: Halted because roundoff errors limited progress."
    algorithm = NLOPT_LN_NELDERMEAD;  
  } else if(optimization_method == "PRINT") { // Print energy and force results at current parameters
    algorithm = OPT_PRINT_ALL;
  } else if(optimization_method == "RAND_SBPLX") { // Print energy and force results at current parameters
    algorithm = OPT_RAND_SBPLX;
  } else {
    optimization_method = "SBPLX";  // SBPLX (local, based on Nelder-Mead simplex)
    algorithm = NLOPT_LN_SBPLX;
  }
  
  printf("Run name = %s\n", run_name.c_str());
  printf("Optimization method = %s\n", optimization_method.c_str());
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
  
  printf("Parameters to optimize: %lu\n", params_current.size());
  
  //get ready to record results
  best_error = -1;
  counter_since_last_file_write = 0;
  ready_to_write_file = 0;
}

/* ---------------------------------------------------------------------- */

int MinParams::iterate(int maxiter)
{
  run_NLopt();
  return MAXEVAL;
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
        char *token = strtok( (type==0?((char*)line.c_str()+charge_label.size()):NULL), " "); //strtok needs first arg to be NULL on all calls after the first
        charges[type] = force->numeric(FLERR,token);
      }
    }
    if(line.rfind(lj_sigma_label, 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
        char *token = strtok( (type==0?((char*)line.c_str()+lj_sigma_label.size()):NULL), " "); //strtok needs first arg to be NULL on all calls after the first
        lj_sigma[type] = force->numeric(FLERR,token);
      }
    }
    if(line.rfind(lj_epsilon_label, 0) == 0) {
      for(int type=0; type < atom->ntypes; type++) {
        char *token = strtok( (type==0?((char*)line.c_str()+lj_epsilon_label.size()):NULL), " "); //strtok needs first arg to be NULL on all calls after the first
        lj_epsilon[type] = force->numeric(FLERR,token);
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
  
  if(algorithm == OPT_PRINT_ALL) { //for testing: set as true to print force and energy comparison, then exit
    for(unsigned int i=0; i<target_forces.size(); i+=3) {
      double dfx = f[i/3][0] - target_forces[i];
      double dfy = f[i/3][1] - target_forces[i+1];
      double dfz = f[i/3][2] - target_forces[i+2];
      printf("F: %g %g %g vs. %g %g %g\n", f[i/3][0], f[i/3][1], f[i/3][2], target_forces[i], target_forces[i+1], target_forces[i+2]);
    }
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

#define iff(X) if( !( tersoff_upper_bound->params[index_ijk].X > tersoff_lower_bound->params[index_ijk].X ) )

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
  }
  //enforce Tersoff constraints: use mixing rules from https://www.quantumwise.com/documents/manuals/latest/ReferenceManual/index.html/ref.tersoffmixitpotential.html
  for(int i = 0; i < tersoff->nelements; i++) { //assumes Tersoff elements come first in type list!
    for(int j = 0; j < tersoff->nelements; j++) {
      for(int k = 0; k < tersoff->nelements; k++) {
        int index_iii = tersoff->elem2param[i][i][i];
        int index_jjj = tersoff->elem2param[j][j][j];
        int index_ijj = tersoff->elem2param[i][j][j];
        int index_ijk = tersoff->elem2param[i][j][k];
        //pairwise mixing: 
        iff(lam1) tersoff->params[index_ijk].lam1 = 0.5*(tersoff->params[index_iii].lam1 + tersoff->params[index_jjj].lam1);
        iff(lam2) tersoff->params[index_ijk].lam2 = 0.5*(tersoff->params[index_iii].lam2 + tersoff->params[index_jjj].lam2);
        iff(biga) tersoff->params[index_ijk].biga = sqrt(tersoff->params[index_iii].biga * tersoff->params[index_jjj].biga);
        iff(bigb) tersoff->params[index_ijk].bigb = sqrt(tersoff->params[index_iii].bigb * tersoff->params[index_jjj].bigb);
        iff(bigd) tersoff->params[index_ijk].bigd = sqrt(tersoff->params[index_iii].bigd * tersoff->params[index_jjj].bigd);
        iff(bigr) tersoff->params[index_ijk].bigr = sqrt(tersoff->params[index_iii].bigr * tersoff->params[index_jjj].bigr);
		lj->cut_inner[i+1][j+1] = tersoff->params[index_ijk].bigr; //copy Tersoff cutoff into pair inout
        //directly copied parameters
        iff(beta) tersoff->params[index_ijk].beta = tersoff->params[index_iii].beta;
        iff(powern) tersoff->params[index_ijk].powern = tersoff->params[index_iii].powern;
        iff(c) tersoff->params[index_ijk].c = tersoff->params[index_iii].c;
        iff(d) tersoff->params[index_ijk].d = tersoff->params[index_iii].d;
        iff(h) tersoff->params[index_ijk].h = tersoff->params[index_iii].h;
        //3-wise:
        iff(gamma) tersoff->params[index_ijk].gamma = tersoff->params[index_ijj].gamma;
        iff(lam3) tersoff->params[index_ijk].lam3 = sqrt(tersoff->params[index_iii].lam3*tersoff->params[index_jjj].lam3); //tersoff->params[index_ijj].lam3 ?
        iff(powerm) tersoff->params[index_ijk].powerm = tersoff->params[index_ijj].powerm;
      }
    }
  }
  //recalculate resulting parameters in Tersoff param objects: cut, cutsq, c1, c2, c3, c4, and cutmax
  for(int i=0; i<tersoff->nparams; i++) {
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
  //enforce charge constraint: Only charge on type 0 matters, q_1 = -0.5*q_0
  charges_current[1] = -0.5*charges_current[0];
  //put new charge on each atom
  for(int a=0; a < atom->natoms; a++) {
    atom->q[a] = charges_current[atom->type[a]-1]; //todo: make sure there are no resulting parameters based on charge that need to be recalculated after changing charge
  }
  //modify LJ parameters
  for(int type=0; type < atom->ntypes; type++) {
    lj->sigma[type+1][type+1] = lj_sigma_current[type]; //lj->sigma and lj->epsilon are 1-indexed
    lj->epsilon[type+1][type+1] = lj_epsilon_current[type];
	//modify rcut inner here?
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

double MinParams::NLopt_target_function(unsigned params_count, const double *params) {
  bigint ntimestep;
  ntimestep = ++update->ntimestep;
  niter++;
  
  unpack_params(params_current); //put parameters into LAMMPS objects

  ecurrent = energy_force(0);
  neval++;


  update->eflag_atom = update->ntimestep; //needed to avoid error message. TODO: make this less of a hack. 
  compute_pe->compute_peratom(); //also sort of a hack


  if(best_error<0) {
    best_error = calculate_error();
    for(unsigned int i=0; i<params_best.size(); i++) 
      params_best[i] = params[i]; //copy into params_best
    write_tersoff_file();
  }
  double new_error = calculate_error();
  
  if(new_error < best_error) {
    best_error = new_error;
    ready_to_write_file = 1;
    //printf("New best: %f\n", best_error );
  }
  counter_since_last_file_write++;
  if(ready_to_write_file && counter_since_last_file_write>10000) {
    for(unsigned int i=0; i<params_best.size(); i++) 
      params_best[i] = params[i]; //copy into params_best
    write_tersoff_file();
  }

  // output for thermo, dump, restart files
  if(output->next == ntimestep) {
    timer->stamp();
    output->write(ntimestep);
    timer->stamp(Timer::OUTPUT);
  }
  
  return new_error;
}

double NLopt_target_function_wrapper(unsigned params_count, const double *params, double *gradient, void *optional_data) {
  MinParams *me = (MinParams *) optional_data;
  return me->NLopt_target_function(params_count, params);
}

void MinParams::run_NLopt() {
  puts("Running NLopt optimization...");
  
  nlopt_opt opt;
  if(algorithm==OPT_RAND_SBPLX) {
    for(unsigned int i=0; i<params_current.size(); i++)
      params_current[i] = params_lower[i] + random->uniform()*(params_upper[i]-params_lower[i]);
    opt = nlopt_create(NLOPT_LN_SBPLX, params_current.size());
  }
  else
    opt = nlopt_create(algorithm, params_current.size());
  nlopt_srand(random_seed);
  
  //enforce requirement of NLopt, that guess must be within the bounds
  for(unsigned int i=0; i<params_current.size(); i++) {
    assert(params_upper[i]>=params_current[i]);
    assert(params_current[i]>=params_lower[i]);
  }
  double tolerance = 1e-15;
  nlopt_set_lower_bounds(opt, &params_lower[0]);
  nlopt_set_upper_bounds(opt, &params_upper[0]);
  nlopt_set_maxeval(opt, update->max_eval);
  nlopt_set_ftol_abs(opt, tolerance);
  nlopt_set_min_objective(opt, NLopt_target_function_wrapper, this);
  
  double error; // the minimum objective value upon return
  
  int return_code = nlopt_optimize(opt, &params_current[0], &error);
  if(return_code < 0) {
    printf("NLopt failed with code %d\n", return_code);
  }
  else {
    printf("NLopt finished with code %d\n", return_code);
  }
}

