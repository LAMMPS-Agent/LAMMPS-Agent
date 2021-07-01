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
   Modified by Changhao Li (czl478@psu.edu, changhaoli1997@gmail.com) for
   modeling of bacteria film. Last modified date: 06/30/2020
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nve/body/agent,FixNVEBodyAgent)

#else

#ifndef LMP_FIX_NVE_BODY_AGENT_H
#define LMP_FIX_NVE_BODY_AGENT_H

#include <random>
#include <vector>
#include "fix_nve.h"

#include "util_fibrin.h"


namespace LAMMPS_NS {

class FixNVEBodyAgent : public FixNVE {
 public:
  FixNVEBodyAgent(class LAMMPS *, int, char **);
  void init();
  int setmask();
  void initial_integrate(int);
  void pre_exchange();
  void final_integrate();


 private:
  double dtq;                    // timestep
  double growth_rate;            // expectation of growth rate, unit is 1/[T]
  double standard_dev;           // standard derivation of growth rate
  double L_critical;             // critical length for proliferation
  double nu_0;                   // damping constant of ambient environments

  double del_height;             // if a cell is above z=del_height plane, then delete it
  double frozen_radius;          // radius of the frozen area (cell can move but not grow)
  double m_radial;               // keep radial alignments
  double mutant_fraction;         // fraction of mutant (zero friction) cell
  double noise_level;            // pre-defined noise level applying on both force and moment vector

  int nsteps;                    // recording current timestep

  double element[4][2] = {{2, 2}, {-2, 2}, {-2, -2}, {2, -2}};

  std::default_random_engine generator;
  std::normal_distribution<double> distribution;
  std::vector<double> list_growth_rate;

  class AtomVecBody *avec;
  class AtomVec *avec_current;
  void apply_damping_force(int, double*, double**, double**);
  void body2space(double*, double*, double*);
  void grow_single_body(int ibody, double growth_ratio);
  void translate_single_body(int ibody, double* vec);
  double length(double* coords);
  double radius(double* data, int i);
  double volume(double r, double L);
  bool out_of_z_plane(Cell* icell, double height);
  bool need_reneighbor();
  double* rot_inertia(double r, double L, double* mom);
  void grow_all_body(double given_growth_ratio = 0);
  void proliferate_all_body();
  void add_noise(double* f, double* mom, double noise_level);
  void add_tension(double* f, double* mom, double A);
  void set_force(int ibody, double fx, double fy, double fz, double tx, double ty, double tz);
  bool freeze(int ibody, double R);
  void radial_moment(int ibody, double M);

  // benchmarking
  void deform_single_element_tensile(double element[][2], double* r1, double* r2, double* c1, double* c2, double dl);
  void apply_gel_cell_adhesion(double* f, double *torque, double* r1, double* r2, double* c1, double* c2);
  void apply_hypothetical_confinements(int ibody, double R, double H);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Fix nve/body requires atom style body

Self-explanatory.

E: Fix nve/body requires bodies

This fix can only be used for particles that are bodies.

*/
