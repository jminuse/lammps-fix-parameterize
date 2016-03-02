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

namespace LAMMPS_NS {

class FixParameterize : public Fix {
 public:
  FixParameterize(class LAMMPS *, int, char **);
  virtual ~FixParameterize() {}
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:
  int random_seed;
  int n_target_forces;
  double *target_forces;
  RanMars *random;
  double best_error;
  int counter_since_last_file_write;
  int ready_to_write_file;
  
  //tersoff values, per pair of atom types
  PairTersoff *tersoff;
  PairTersoff::Param *best_tersoff_params;
  
  //Charges and LJ, per atom type
  PairLJCutCoulInOut *lj;
  double *best_charges;
  double *best_lj_sigma;
  double *best_lj_epsilon;
  
  PairTersoff *tersoff_upper_bound;
  PairTersoff *tersoff_lower_bound;
  
  double calculate_error();
  void write_tersoff_file();
  void pack_params();
  void unpack_params();
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
