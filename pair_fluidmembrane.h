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

/* ----------------------------------------------------------------------
   Contributing author: Hongyan Yuan (URI), June, 2015
   Modified by Changhao Li for modeling of soft gel
   Last modified: 07/13/2020
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(fluidmembrane,Pairfluidmembrane)

#else

#ifndef LMP_PAIR_fluidmembrane_H
#define LMP_PAIR_fluidmembrane_H

#include "pair.h"

namespace LAMMPS_NS {

class Pairfluidmembrane : public Pair {
 public:
  Pairfluidmembrane(LAMMPS *lmp);
  virtual ~Pairfluidmembrane();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);


 protected:

  double cut_global;

  double **epsilon,**sigma,**cut,**zeta,**mu,**beta;  // model parameter values for atom-type pairs

  class AtomVecEllipsoid *avec;

  void allocate();
  double fluidmembrane_analytic(const int i, const int j, double a1[3][3],
                           double a2[3][3], double *r12,
                           const double rsq, double *fforce, double *ttor,
                           double *rtor);


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

E: Pair fluidmembrane requires atom style ellipsoid

Self-explanatory.

*/
