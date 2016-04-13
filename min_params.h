/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   This file is a (heavy) modification of min_fire.h
------------------------------------------------------------------------- */

#ifdef MINIMIZE_CLASS

MinimizeStyle(params,MinParams)

#else

#ifndef LMP_MIN_PARAMS_H
#define LMP_MIN_PARAMS_H

#include "min.h"

#include "random_mars.h"
#include "pair_tersoff.h"
#include "pair_lj_cut_coul_inout.h"

#include <vector>
#include <string>

#include <nlopt.h> //NLopt Non-Linear Optimization library. Compilation flags: gcc -I/fs/home/jms875/install/include -L/fs/home/jms875/install/lib -lnlopt -lm

namespace LAMMPS_NS {

class MinParams : public Min {
 public:
  MinParams(class LAMMPS *);
  ~MinParams() {}
  void init();
  void setup_style();
  void reset_vectors();
  int iterate(int);
  void modify_params(int, char **); // to set run_name, etc
  double NLopt_target_function(unsigned params_count, const double *params); //function which gets evaluated by the NLopt library

 protected:
  int random_seed;
  RanMars *random;
  nlopt_algorithm algorithm;
  double best_error;
  double force_error, energy_error;
  std::vector<double> target_forces;
  std::vector<double> target_energies;
  std::vector<double> current_energies;
  //double* atomwise_energies;
  int counter_since_last_file_write;
  int ready_to_write_file;
  std::string run_name;
  
  Compute *compute_pe;
  
  //tersoff values, per pair of atom types
  PairTersoff *tersoff;
  
  //Charges and LJ, per atom type
  PairLJCutCoulInOut *lj;
  std::vector<double> charges_current;
  std::vector<double> lj_sigma_current;
  std::vector<double> lj_epsilon_current;
  
  //bounds
  PairTersoff *tersoff_upper_bound;
  std::vector<double> charges_upper;
  std::vector<double> lj_sigma_upper;
  std::vector<double> lj_epsilon_upper;
  PairTersoff *tersoff_lower_bound;
  std::vector<double> charges_lower;
  std::vector<double> lj_sigma_lower;
  std::vector<double> lj_epsilon_lower;
  
  //Flat array of all parameters made by pack_params
  std::vector<double> params_best;
  std::vector<double> params_upper;
  std::vector<double> params_lower;
  std::vector<double> params_current;
  
  //Functions
  double calculate_error();
  void write_tersoff_file();
  void pack_params();
  void unpack_params(std::vector<double> pp);
  void read_params_from_comments(std::string filename, std::vector<double> &charges, std::vector<double> &lj_sigma, std::vector<double> &lj_epsilon);
  void run_NLopt(int maxiter);
};
}

#endif
#endif
