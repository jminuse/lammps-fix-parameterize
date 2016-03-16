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

FixStyle(parameterize,FixParameterize)

#else

#ifndef LMP_FIX_PARAMETERIZE_H
#define LMP_FIX_PARAMETERIZE_H

#include "fix.h"

#include "random_mars.h"
#include "pair_tersoff.h"
#include "pair_lj_cut_coul_inout.h"

#include <vector>

namespace LAMMPS_NS {

class FixParameterize : public Fix {
 public:
  FixParameterize(class LAMMPS *, int, char **);
  virtual ~FixParameterize() {}
  int setmask();
  virtual void init();
  virtual void final_integrate();

 protected:
  int random_seed;
  RanMars *random;
  double best_error;
  std::vector<double> target_forces;
  std::vector<double> target_energies;
  std::vector<double> current_energies;
  double* atomwise_energies;
  int counter_since_last_file_write;
  int ready_to_write_file;
  char *output_filename;
  
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
  void read_params_from_comments(const char *filename, std::vector<double> &charges, std::vector<double> &lj_sigma, std::vector<double> &lj_epsilon);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
