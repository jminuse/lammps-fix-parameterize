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

#ifdef PAIR_CLASS

PairStyle(lj/cut/coul/inout,PairLJCutCoulInOut)

#else

#ifndef LMP_PAIR_LJ_CUT_COUL_IN_OUT_H
#define LMP_PAIR_LJ_CUT_COUL_IN_OUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulInOut : public Pair {
 friend class MinParams; //allows MinParams to modify protected parameters
 
 public:
  PairLJCutCoulInOut(class LAMMPS *);
  ~PairLJCutCoulInOut();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_inner_global;
  double **cut_inner,**cut_inner_sq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;

  double cut_outer,cut_outer_sq;
  double alpha;
  double f_shift,e_shift;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/cut/coul/inout requires atom attribute q

The atom style defined does not have these attributes.

*/
