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
#define mark_as_unused_for_compiler(x) (void)(x)

// same as in other min classes

enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,ZEROQUAD};

/* ---------------------------------------------------------------------- */

MinParams::MinParams(LAMMPS *lmp) : Min(lmp) {}

#define OPT_PRINT_ALL NLOPT_LD_LBFGS //used as code for "PRINT"
#define OPT_RAND_SBPLX NLOPT_LD_SLSQP //used as code for RAND_SBPLX
#define OPT_GAUSS_SBPLX NLOPT_LN_PRAXIS //used as code for GAUSS_SBPLX

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
  } else if(optimization_method == "RAND_SBPLX") { // Generate uniformly distributed random solutions and locally optimize with SBPLX
    algorithm = OPT_RAND_SBPLX;
  } else if(optimization_method == "GAUSS_SBPLX") { // Generate Gaussian distributed random solutions near the current best and locally optimize with SBPLX
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

  //make random number generator
  random = new RanMars(lmp, random_seed + comm->me);
  
  //get Tersoff calculator (must be declared in LAMMPS script)
  tersoff = (PairTersoff*) force->pair_match("tersoff",1,0);
  if(tersoff==NULL) error->all(FLERR,"Can't find Tersoff pair style in this run");

  //read Tersoff bounds
  tersoff_upper_bound = new PairTersoff(lmp);
  tersoff_lower_bound = new PairTersoff(lmp);
  const char *tersoff_bounds_args[atom->ntypes+3];
  tersoff_bounds_args[0] = "*";
  tersoff_bounds_args[1] = "*";
  tersoff_bounds_args[2] = upper_bounds_filename.c_str();
  for(int type=1; type <= atom->ntypes; type++) {
    int tersoff_type = tersoff->map[type];
    if(tersoff_type>=0)
      tersoff_bounds_args[type+2] = tersoff->elements[tersoff_type];
    else
      tersoff_bounds_args[type+2] = "NULL";
  }
  tersoff_upper_bound->coeff(atom->ntypes+3, (char**)tersoff_bounds_args);
  tersoff_bounds_args[2] = lower_bounds_filename.c_str();
  tersoff_lower_bound->coeff(atom->ntypes+3, (char**)tersoff_bounds_args);
  
  //read other upper bounds
  read_params_from_comments(upper_bounds_filename, charges_upper, lj_sigma_upper, lj_epsilon_upper);
  read_params_from_comments(lower_bounds_filename, charges_lower, lj_sigma_lower, lj_epsilon_lower);
  read_params_from_comments(input_tersoff_filename, charges_current, lj_sigma_current, lj_epsilon_current);
  
  int id = modify->find_compute("atom_pe"); //must use this name for the compute!
  if (id < 0) error->all(FLERR,"Minimization could not find pe/atom compute (must be named atom_pe)");
  compute_pe = modify->compute[id];

  printf("tersoff nparams = %d\n", tersoff->nparams);
  for(int i=0; i<tersoff->nparams; i++) {
    printf("tersoff elements = %s,%s,%s, A = %f\n", tersoff->elements[ tersoff->params[i].ielement ], tersoff->elements[ tersoff->params[i].jelement ], tersoff->elements[ tersoff->params[i].kelement ], tersoff->params[i].biga);
  }
  for(int type=0; type < tersoff->nelements; type++) { //for testing only
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
  maxiter = std::max(update->max_eval, maxiter);
  run_NLopt(maxiter);
  return maxiter;
}


void MinParams::read_params_from_comments(std::string filename, std::vector<double> &charges, std::vector<double> &lj_sigma, std::vector<double> &lj_epsilon) {
  charges.assign(tersoff->nelements, 0.0);
  lj_sigma.assign(tersoff->nelements, 0.0);
  lj_epsilon.assign(tersoff->nelements, 0.0);
  
  std::ifstream infile(filename.c_str());
  std::string line;
  std::string charge_label = "# Charges:";
  std::string lj_sigma_label = "# LJ-sigma:";
  std::string lj_epsilon_label = "# LJ-epsilon:";
  while(std::getline(infile, line)) {
    if(line.rfind(charge_label, 0) == 0) {
      for(int type=0; type < tersoff->nelements; type++) {
        char *token = strtok( (type==0?((char*)line.c_str()+charge_label.size()):NULL), " "); //strtok needs first arg to be NULL on all calls after the first
        charges[type] = force->numeric(FLERR,token);
      }
    }
    if(line.rfind(lj_sigma_label, 0) == 0) {
      for(int type=0; type < tersoff->nelements; type++) {
        char *token = strtok( (type==0?((char*)line.c_str()+lj_sigma_label.size()):NULL), " "); //strtok needs first arg to be NULL on all calls after the first
        lj_sigma[type] = force->numeric(FLERR,token);
      }
    }
    if(line.rfind(lj_epsilon_label, 0) == 0) {
      for(int type=0; type < tersoff->nelements; type++) {
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
    double small_magnitude = 1e-4; //kcal/mol/Angstrom
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
    error += 1e6*isnan(current_energies[i]); //massive penalty for NaN errors
    energy_error += error*error;
  }
  
  if(algorithm == OPT_PRINT_ALL) { //for testing: set as true to print force and energy comparison, then exit
    for(unsigned int i=0; i<target_forces.size(); i+=3) {
      printf("F: %g %g %g vs. %g %g %g\n", f[i/3][0], f[i/3][1], f[i/3][2], target_forces[i], target_forces[i+1], target_forces[i+2]);
    }
    for(unsigned int i=0; i<target_energies.size(); i++) {
      printf("E: %g %g %d\n", target_energies[i], current_energies[i], isnan(current_energies[i]));
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
  FILE *output = fopen(output_filename.c_str(), "w");
  //label
  fprintf(output, "# Generated by min_style params, run = '%s'. LAMMPS units = real.\n", run_name.c_str());
  //output error
  fprintf(output, "# Error metric=%g\n", best_error);
  //output charges on one line, listed by type
  fprintf(output, "# Charges:    ");
  for(int i=1; i <= atom->ntypes; i++) if(tersoff->map[i]>=0) fprintf(output, " %9.6g ", charges_current[tersoff->map[i]]);
  fprintf(output, "\n");
  //output LJ sigma on one line, listed by type
  fprintf(output, "# LJ-sigma:   ");
  for(int i=1; i <= atom->ntypes; i++) if(tersoff->map[i]>=0) fprintf(output, " %9.6g ", lj_sigma_current[tersoff->map[i]]);
  fprintf(output, "\n");
  //output LJ epsilon on one line, listed by type
  fprintf(output, "# LJ-epsilon: ");
  for(int i=1; i <= atom->ntypes; i++) if(tersoff->map[i]>=0) fprintf(output, " %9.6g ", lj_epsilon_current[tersoff->map[i]]);
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
  for(int type=0; type < tersoff->nelements; type++) {
    pack_typewise_param(charges_upper, charges_lower, charges_current);
    pack_typewise_param(lj_sigma_upper, lj_sigma_lower, lj_sigma_current);
    pack_typewise_param(lj_epsilon_upper, lj_epsilon_lower, lj_epsilon_current);
  }
  
}

#define unpack_param(X) if( tersoff_upper_bound->params[i].X > tersoff_lower_bound->params[i].X ) { tersoff->params[i].X = params[param_i]; param_i++; }

#define iff(X) if( !( tersoff_upper_bound->params[index_ijk].X > tersoff_lower_bound->params[index_ijk].X ) )

void MinParams::unpack_params(std::vector<double> params) {
  int param_i = 0;
  //unpack Tersoff parameters from the flat array to the LAMMPS object
  for(int i=0; i<tersoff->nparams; i++) {
    unpack_param(lam1);
    unpack_param(lam2);
    unpack_param(lam3);
    unpack_param(c);
    unpack_param(d);
    unpack_param(h);
    unpack_param(gamma);
    unpack_param(powerm);
    unpack_param(powern);
    unpack_param(beta);
    unpack_param(biga);
    unpack_param(bigb);
    unpack_param(bigd);
    unpack_param(bigr);
    unpack_param(cut);
  }
  //enforce Tersoff constraints
  for(int i_type = 1; i_type <= atom->ntypes; i_type++) {
    for(int j_type = 1; j_type <= atom->ntypes; j_type++) {
      for(int k_type = 1; k_type <= atom->ntypes; k_type++) {
        int i_tersoff = tersoff->map[i_type];
        int j_tersoff = tersoff->map[j_type];
        int k_tersoff = tersoff->map[k_type];
        if(i_tersoff<0 || j_tersoff<0 || k_tersoff<0) continue;
        int index_iii = tersoff->elem2param[i_tersoff][i_tersoff][i_tersoff];
        int index_jjj = tersoff->elem2param[j_tersoff][j_tersoff][j_tersoff];
        int index_ijj = tersoff->elem2param[i_tersoff][j_tersoff][j_tersoff];
        int index_ikk = tersoff->elem2param[i_tersoff][k_tersoff][k_tersoff];
        int index_ijk = tersoff->elem2param[i_tersoff][j_tersoff][k_tersoff];
        //char *element_i = tersoff->elements[i_tersoff];
        //char *element_j = tersoff->elements[j_tersoff];
        //char *element_k = tersoff->elements[k_tersoff];

        if(tersoff->params[index_iii].biga > 0.0) { // mix iii and jjj parameters, using standard mixing rules per https://www.quantumwise.com/documents/manuals/latest/ReferenceManual/index.html/ref.tersoffmixitpotential.html
          //pairwise mixing:
          iff(lam1) tersoff->params[index_ijk].lam1 = 0.5*(tersoff->params[index_iii].lam1 + tersoff->params[index_jjj].lam1);
          iff(lam2) tersoff->params[index_ijk].lam2 = 0.5*(tersoff->params[index_iii].lam2 + tersoff->params[index_jjj].lam2);
          iff(biga) tersoff->params[index_ijk].biga = sqrt(tersoff->params[index_iii].biga * tersoff->params[index_jjj].biga);
          iff(bigb) tersoff->params[index_ijk].bigb = sqrt(tersoff->params[index_iii].bigb * tersoff->params[index_jjj].bigb);
          iff(bigd) tersoff->params[index_ijk].bigd = sqrt(tersoff->params[index_iii].bigd * tersoff->params[index_jjj].bigd);
          iff(bigr) tersoff->params[index_ijk].bigr = sqrt(tersoff->params[index_iii].bigr * tersoff->params[index_jjj].bigr);
          //directly copied parameters:
          iff(beta) tersoff->params[index_ijk].beta = tersoff->params[index_iii].beta;
          iff(powern) tersoff->params[index_ijk].powern = tersoff->params[index_iii].powern;
          iff(c) tersoff->params[index_ijk].c = tersoff->params[index_iii].c;
          iff(d) tersoff->params[index_ijk].d = tersoff->params[index_iii].d;
          iff(h) tersoff->params[index_ijk].h = tersoff->params[index_iii].h;
          //3-wise parameters, treated as ijj
          iff(gamma) tersoff->params[index_ijk].gamma = tersoff->params[index_ijj].gamma;
          iff(lam3) tersoff->params[index_ijk].lam3 = tersoff->params[index_ijj].lam3;
          iff(powerm) tersoff->params[index_ijk].powerm = tersoff->params[index_ijj].powerm;
        }
        else if(tersoff->params[index_ijj].biga > 0.0 && tersoff->params[index_ikk].biga > 0.0) { // use ijj parameters mixed with ikk parameters
          //printf("Bonds: %s--%s--%s = mix(%s--%s--%s, %s--%s--%s)\n", element_j, element_i, element_k,    element_j, element_i, element_j,    element_k, element_i, element_k);
          //
          iff(lam1) tersoff->params[index_ijk].lam1 = sqrt(tersoff->params[index_ijj].lam1 * tersoff->params[index_ikk].lam1);
          iff(lam2) tersoff->params[index_ijk].lam2 = sqrt(tersoff->params[index_ijj].lam2 * tersoff->params[index_ikk].lam2);
          iff(biga) tersoff->params[index_ijk].biga = sqrt(tersoff->params[index_ijj].biga * tersoff->params[index_ikk].biga);
          iff(bigb) tersoff->params[index_ijk].bigb = sqrt(tersoff->params[index_ijj].bigb * tersoff->params[index_ikk].bigb);
          iff(bigd) tersoff->params[index_ijk].bigd = 0.5*(tersoff->params[index_ijj].bigd + tersoff->params[index_ikk].bigd);
          iff(bigr) tersoff->params[index_ijk].bigr = sqrt(tersoff->params[index_ijj].bigr * tersoff->params[index_ikk].bigr);
          //
          iff(beta) tersoff->params[index_ijk].beta = sqrt(tersoff->params[index_ijj].beta * tersoff->params[index_ikk].beta);
          iff(powern) tersoff->params[index_ijk].powern = 0.5*(tersoff->params[index_ijj].powern + tersoff->params[index_ikk].powern);
          iff(c) tersoff->params[index_ijk].c = sqrt(tersoff->params[index_ijj].c * tersoff->params[index_ikk].c);
		  iff(d) tersoff->params[index_ijk].d = 0.5*(tersoff->params[index_ijj].d + tersoff->params[index_ikk].d);
          iff(h) tersoff->params[index_ijk].h = 0.5*(tersoff->params[index_ijj].h + tersoff->params[index_ikk].h);
          //3-wise parameters, mixed ijj and ikk
          iff(gamma) tersoff->params[index_ijk].gamma = sqrt(tersoff->params[index_ijj].gamma * tersoff->params[index_ikk].gamma);
          iff(lam3) tersoff->params[index_ijk].lam3 = sqrt(tersoff->params[index_ijj].lam3 * tersoff->params[index_ikk].lam3);
          iff(powerm) tersoff->params[index_ijk].powerm = sqrt(tersoff->params[index_ijj].powerm * tersoff->params[index_ikk].powerm);
        }
        else {
          //printf("No bonds: %s--%s--%s\n", element_j, element_i, element_k);
        }
      }
    }
  }
  //puts("Exiting");
  //exit(0);
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
  //copy Tersoff cutoff into pair inout, if Tersoff cutoff exists. Use ijj. 
  for(int i_type = 1; i_type <= atom->ntypes; i_type++) {
    for(int j_type = 1; j_type <= atom->ntypes; j_type++) {
      int i_tersoff = tersoff->map[i_type];
      int j_tersoff = tersoff->map[j_type];
      if(i_tersoff<0 || j_tersoff<0) continue;
      int index_ijj = tersoff->elem2param[i_tersoff][j_tersoff][j_tersoff];
      if(tersoff->params[index_ijj].cut > 0.0)
		  lj->cut_inner[i_type][j_type] = tersoff->params[index_ijj].cut;
    }
  }
  //unpack other parameters from the flat array to the LAMMPS object
  for(int atom_type=1; atom_type <= tersoff->nelements; atom_type++) {
    int tersoff_type = tersoff->map[atom_type];
    if(tersoff_type >= 0) {
      if( charges_upper[tersoff_type]    > charges_lower[tersoff_type]    ) { charges_current[tersoff_type]    = params[param_i]; param_i++; }
      if( lj_sigma_upper[tersoff_type]   > lj_sigma_lower[tersoff_type]   ) { lj_sigma_current[tersoff_type]   = params[param_i]; param_i++; }
      if( lj_epsilon_upper[tersoff_type] > lj_epsilon_lower[tersoff_type] ) { lj_epsilon_current[tersoff_type] = params[param_i]; param_i++; }
    }
  }
  //enforce charge constraint: Only charge on type 0 matters, Q_1 = -0.5*Q_0. Assumes system is PbCl2.
  charges_current[1] = -0.5*charges_current[0];
  //put new charge on each atom
  for(int a=0; a < atom->natoms; a++) {
    atom->q[a] = charges_current[atom->type[a]-1]; //todo: make sure there are no resulting parameters based on charge that need to be recalculated after changing charge
  }
  //modify LJ parameters
  for(int atom_type=1; atom_type <= atom->ntypes; atom_type++) {
    int tersoff_type = tersoff->map[atom_type];
    if(tersoff_type >= 0) {
      lj->sigma[atom_type][atom_type] = lj_sigma_current[tersoff_type]; //lj->sigma and lj->epsilon are 1-indexed
      lj->epsilon[atom_type][atom_type] = lj_epsilon_current[tersoff_type];
	  //lj->cut_inner modified above in Tersoff loop (bigr)
	}
  }
  //recalculate resulting parameters
  for(int i=1; i <= atom->ntypes; i++) {
    for(int j=i; j <= atom->ntypes; j++) {
      lj->setflag[i][j] = 0; //unset, or init_one (below) is ignored
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
  mark_as_unused_for_compiler(params_count); //mark as unused for compiler
  bigint ntimestep = ++update->ntimestep;
  niter++;
  neval++;
  
  unpack_params(params_current); //put parameters into LAMMPS objects

  modify->addstep_compute_all(update->ntimestep);
  ecurrent = energy_force(0);
  compute_pe->compute_peratom(); //sort of a hack
  

  //printf("force->pair->eatom[0] = %g\n", force->pair->eatom[0] );

  if(best_error<0) {
    best_error = calculate_error();
    for(unsigned int i=0; i<params_best.size(); i++) 
      params_best[i] = params[i]; //copy into params_best
    write_tersoff_file();
  }
  double new_error = calculate_error();
  //printf("%f %f %f\n", best_error, new_error, ecurrent);
  if(new_error < best_error) {
    best_error = new_error;
    ready_to_write_file = 1;
    for(unsigned int i=0; i<params_best.size(); i++) 
      params_best[i] = params[i]; //copy into params_best
    //printf("New best: %f\n", best_error );
  }
  counter_since_last_file_write++;
  if(ready_to_write_file && counter_since_last_file_write>10000) {
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
  mark_as_unused_for_compiler(gradient);
  MinParams *me = (MinParams *) optional_data;
  return me->NLopt_target_function(params_count, params);
}

void MinParams::run_NLopt(int maxiter) {
  puts("Running NLopt optimization...");
  
  nlopt_opt opt;
  if(algorithm==OPT_RAND_SBPLX || algorithm==OPT_GAUSS_SBPLX)
    opt = nlopt_create(NLOPT_LN_SBPLX, params_current.size());
  else
    opt = nlopt_create(algorithm, params_current.size());
  nlopt_srand(random_seed);
  
  //enforce requirement of NLopt, that guess must be within the bounds
  for(unsigned int i=0; i<params_current.size(); i++) {
    if(params_upper[i]<params_current[i] || params_current[i]<params_lower[i]) {
      printf("Condition %g>%g>%g violated for param %d\n", params_upper[i], params_current[i], params_lower[i], i);
      error->all(FLERR,"Should be upper > current > lower");
    }
  }
  double tolerance = 1e-15;
  nlopt_set_lower_bounds(opt, &params_lower[0]);
  nlopt_set_upper_bounds(opt, &params_upper[0]);
  nlopt_set_maxeval(opt, maxiter);
  nlopt_set_ftol_abs(opt, tolerance);
  nlopt_set_min_objective(opt, NLopt_target_function_wrapper, this);
  
  double error; // the minimum objective value upon return
  int return_code; 
  
  if(algorithm==OPT_RAND_SBPLX || algorithm==OPT_GAUSS_SBPLX) {
    while(1) {
      nlopt_set_ftol_abs(opt, 1e-6); //set relatively loose tolerance, to allow multiple local optimizations
      for(unsigned int i=0; i<params_current.size(); i++) {
        if(algorithm==OPT_RAND_SBPLX)
          params_current[i] = params_lower[i] + random->uniform()*(params_upper[i]-params_lower[i]); //uniform random
        else {
          params_current[i] = params_best[i] + random->gaussian()*0.2*(params_upper[i]-params_lower[i]); //Gaussian random, standard deviation = 0.2 (as in DDS algorithm)
          if(params_upper[i]<params_current[i]) params_current[i] = params_upper[i];
          if(params_current[i]<params_lower[i]) params_current[i] = params_lower[i];
        }
      }
      return_code = nlopt_optimize(opt, &params_current[0], &error);
      printf("Return code=%d, error=%g\n", return_code, error);
    }
  }
  else
    return_code = nlopt_optimize(opt, &params_current[0], &error);
  
  if(return_code < 0) {
    printf("NLopt failed with code %d\n", return_code);
  }
  else {
    printf("NLopt finished with code %d\n", return_code);
  }
}

/* ---------------------------------------------------------------------- */

void MinParams::setup_style() //needed by Min.cpp
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

void MinParams::reset_vectors() //needed by Min.cpp
{
  // atomic dof

  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}
