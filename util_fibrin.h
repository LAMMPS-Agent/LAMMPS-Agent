/* ----------------------------------------------------------------------
   Author: Farzan Beroz, Princeton University
   Github: https://github.com/farzanb/Agent-based-biofilm-model
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Modified by Changhao Li (czl478@psu.edu) for modeling of bacteria film
   Last modified: 07/13/2020
------------------------------------------------------------------------- */
#pragma once
#include <iostream>
#include <complex>
#include <fstream>
#include <string>
#include <math.h>
#include <random>
#include <iomanip>
#include <cmath>
#include <map>
#include <sys/stat.h>
#include <vector>
#include <numeric>

using namespace std;

// The model parameters

// // The following values determine the units:
static const double R = 1; // spherocylinder radius
// static const double E_0 = 1 * R * R; // elastic modulus of the cell
static const double A = 1;	// growth rate

// // Cell-cell and cell-surface interaction parameters:
// //const double A_0 = E_0/1000.0; 	// cell-cell adhesion energy
// static const double A_0 = 0;
// static double E_1 = 100;  	// the elastic modulus of cell-surface interaction (cylindrical term)
// static const double A_1 = E_1/1000; 	// cell-surface adhesion energy (cylindrical term)
// static const double E_2 = E_1;			// the elastic modulus of cell-surface interaction (spherical term)
// static const double A_2 = A_1;			// cell-surface adhesion energy (spherical term)

// // Viscosity parameters:
// static double NU_1 = 0.1 * (E_0) / 1; // surface drag coefficient
// static double NU_0 = NU_1 * E_0 / 100; // ambient (positional) drag

// Additional cell parameters:
static const double Adev  = 0.2; // growth rate noise
static const double L0 = 2*R; // initial length of the first cell
static double Lf; // mean final length

// Simulation parameters
static const double DOF_GROW = 7; // degrees of freedom (growing cell)
static const int MAX_DOF = 5000; // Total allowed degrees of freedom. N * DOF_GROW must be less than MAX_DOF, else segfault!
static const double noise_level = 1e-8; // A small, symmetry breaking noise is chosen from a distribution of this width
// and added to the generalized forces at each step during the evolution.

static const double h0 = 1e-4; // absolute accuracy of diff. eq. solver

static const double T = 12; // Integrate for this amount of total time.
static const double tstep = 0.1; // for output purposes only. Outputs data after every interval of time "tstep"

// Note: pi is given by M_PI

static double xy = 0; // Confine to xy plane?
static const double xz = 0; // Confine to xz plane?

static const int D = 5; // The physical degrees of freedom of a single cell.

static const double hd0 = 0.01;
static const double k0 = 4.5;


struct overlapVec {
	// For two contacting cells, this struct specifies the locations of the points along
	// the cell centerlines that are separated by the smallest distance.
	// The distance is given by "overlap."
    double rix1, riy1, riz1;
	double rix2, riy2, riz2;
	double overlap;
};
struct myVec {
	// Vector
	double x, y, z;
};
struct my6Vec {
	// Cell center-of-mass position (x, y, z) and orientation (nx, ny, nz)
	double x, y, z, nx, ny, nz;
};
struct derivFunc{
	// This struct is used when numerically integrating the equation of motion
	// to store the generalized forces acting on the cell coordinates
	double xd, yd, zd, nxd, nyd, nzd;
};
struct pop_info {
	// When tracking the verticalization of a cell,
	// this struct is used to store the orientation and
	// cell-to-cell contact forces provided by neighboring cells
	double vx, vy, vz;
	double fx, fy, fz;
	double nz;
};

inline double cot(double i) { return(1 / tan(i)); } // Cotangent function
inline double csc(double i) { return(1 / sin(i)); } // Cosecant function

inline double cross(double x1, double y1, double z1, double x2, double y2, double z2, int dof) {
	// Cross product of r1 = (x1, y1, z1) and r2 (x2, y2, z2).
	// Returns the degree of freedom "dof" of the vector r1 x r2
	if (dof == 0) {
		return -(y2*z1) + y1*z2;
	} else if (dof == 1) {
		return x2*z1 - x1*z2;
	} else if (dof == 2) {
		return -(x2*y1) + x1*y2;
	}
	return 0;
}

inline double vnorm(double x, double y, double z) {
	// Return the magnitude of a vector (x, y, z)
	return sqrt(x*x + y*y + z*z);
}

class Cell {
	// One cell. Stores the coordinates and creates daughter cells when time to divide
	private:
		double x, y, z, nx, ny, nz, l, Lfi, ai; //initial node position on grid
		double tborn;
		string lineage;
		int fpop;
		double current_nz;
		int current_sa;

		double dx, dy, dz, dnx, dny, dnz;
		double m1x, m1y, m1z, m2x, m2y, m2z;

		double mf1, mf2, mf3, max_l;
		double of1, of2, of3;
		int confine;
		//Cell* d1_cell;
		//Cell* d2_cell;

		vector<pop_info> pop_vec;

	public:
		Cell (double, double, double, double, double, double, double);
		void set_pos(double mx, double my, double mz) {
			x = mx;
			y = my;
			z = mz;
		}
		void set_n(double mnx, double mny, double mnz) {
			nx = mnx;
			ny = mny;
			nz = mnz;
		}
		void set_l(double ml) {
			l = ml;
		}
		void set_confine() { confine = 1; }
		void set_ai(double mai) { ai = mai; }
		void set_x(double mx) { x = mx; }
		void set_y(double mx) { y = mx; }
		void set_z(double mx) { z = mx; }
		void set_nx(double mx) { nx = mx; }
		void set_ny(double mx) { ny = mx; }
		void set_nz(double mx) { nz = mx; }
		void set_m() {
			double ax = rand()/(static_cast<double>(RAND_MAX));
			double ay = rand()/(static_cast<double>(RAND_MAX));
			double az = rand()/(static_cast<double>(RAND_MAX));
			double cx = cross(nx, ny, nz, ax, ay, az, 0);
			double cy = cross(nx, ny, nz, ax, ay, az, 1);
			double cz = cross(nx, ny, nz, ax, ay, az, 2);
			double cn = vnorm(cx, cy, cz);
			while (cn < 1e-6) {
				ax = rand()/(static_cast<double>(RAND_MAX));
				ay = rand()/(static_cast<double>(RAND_MAX));
				az = rand()/(static_cast<double>(RAND_MAX));
				cx = cross(nx, ny, nz, ax, ay, az, 0);
				cy = cross(nx, ny, nz, ax, ay, az, 1);
				cz = cross(nx, ny, nz, ax, ay, az, 2);
				cn = vnorm(cx, cy, cz);
			}
			m1x = cx / cn;
			m1y = cy / cn;
			m1z = cz / cn;

			m2x = cross(nx, ny, nz, m1x, m1y, m1z, 0);
			m2y = cross(nx, ny, nz, m1x, m1y, m1z, 1);
			m2z = cross(nx, ny, nz, m1x, m1y, m1z, 2);

			/*
			m1x = 0;
			m1y = 1;
			m1z = 0;

			m2x = 0;
			m2y = 0;
			m2z = 1;
			*/
		}
		void set_m1(double pp) {
			nx = nx + m1x * pp;
			ny = ny + m1y * pp;
			nz = nz + m1z * pp;
		}
		void set_m2(double pp) {
			nx = nx + m2x * pp;
			ny = ny + m2y * pp;
			nz = nz + m2z * pp;
		}
		void set_noise() {
		    dx = noise_level * (rand()/(static_cast<double>(RAND_MAX)) - 0.5);
		    dy = noise_level * (rand()/(static_cast<double>(RAND_MAX)) - 0.5);
		    dz = noise_level * (rand()/(static_cast<double>(RAND_MAX)) - 0.5);
		    dnx = noise_level * (rand()/(static_cast<double>(RAND_MAX)) - 0.5);
		    dny = noise_level * (rand()/(static_cast<double>(RAND_MAX)) - 0.5);
		    dnz = noise_level * (rand()/(static_cast<double>(RAND_MAX)) - 0.5);
		}
		void set_tborn(double ti) {
			tborn = ti;
		}
		void set_lin(string xb) {
			lineage = lineage + xb;
		}
		void set_fpop() { fpop = 1; }
		void set_current_nz(double mz) {
			current_nz = mz;
		}
		void set_current_sa(int sa) { current_sa = sa; }
		void set_mf1(double mf) {
			mf1 = mf;
		}
		void set_mf2(double mf) {
			mf2 = mf;
		}
		void set_mf3(double mf) {
			mf3 = mf;
		}
		void set_of1(double mf) {
			of1 = mf;
		}
		void set_of2(double mf) {
			of2 = mf;
		}
		void set_of3(double mf) {
			of3 = mf;
		}
		void set_max_l(double mf) {
			max_l = mf;
		}
		void add_pop_vec(pop_info pi) {
			pop_vec.push_back(pi);
		}
		void clear_pop_vec() {
			pop_vec.clear();
		}
		int get_pop_contact() {
			return pop_vec.size();
		}
		double get_pop_vx(int i) {
			return pop_vec[i].vx;
		}
		double get_pop_vy(int i) {
			return pop_vec[i].vy;
		}
		double get_pop_vz(int i) {
			return pop_vec[i].vz;
		}
		double get_pop_fx(int i) {
			return pop_vec[i].fx;
		}
		double get_pop_fy(int i) {
			return pop_vec[i].fy;
		}
		double get_pop_fz(int i) {
			return pop_vec[i].fz;
		}
		double get_pop_nz(int i) {
			return pop_vec[i].nz;
		}
		double get_mf1() {
			return mf1;
		}
		double get_mf2() {
			return mf2;
		}
		double get_mf3() {
			return mf3;
		}
		double get_of1() {
			return of1;
		}
		double get_of2() {
			return of2;
		}
		double get_of3() {
			return of3;
		}
		double get_max_l() {
			return max_l;
		}
		my6Vec get_noise() {
			my6Vec noise;
			noise.x = dx;
			noise.y = dy;
			noise.z = dz;
			noise.nx = dnx;
			noise.ny = dny;
			noise.nz = dnz;
			return noise;
		}
		// Return the state of the cell
		double get_x() { return x; }
		double get_y() { return y; }
		double get_z() { return z; }
		double get_nx() { return nx; }
		double get_ny() { return ny; }
		double get_nz() { return nz; }
		double get_m1x() { return m1x; }
		double get_m1y() { return m1y; }
		double get_m1z() { return m1z; }
		double get_m2x() { return m2x; }
		double get_m2y() { return m2y; }
		double get_m2z() { return m2z; }
		double get_l() { return l; }
		double get_Lf() { return Lfi; }
		double get_ai() { return ai; }
		double get_tborn() { return tborn; }
		string get_lin() { return lineage; }
		int get_fpop() { return fpop; }
		double get_current_nz() { return current_nz; }
		int get_current_sa() { return current_sa; }
		double get_sf_axis(int dof) {
			double t0, t1, t2, a0, a1, a2, an;

			a0 = cross(0, 0, 1, nx, ny, nz, 0);
			a1 = cross(0, 0, 1, nx, ny, nz, 1);
			a2 = cross(0, 0, 1, nx, ny, nz, 2);

			an = sqrt(a0*a0 + a1*a1 + a2*a2);

			if (dof == 0) {
				return a0/an;
			}
			else if (dof == 1) {
				return a1/an;
			}
			else if (dof == 2) {
				return a2/an;
			} else {
				cout << "SF AXIS ERROR" << endl;
				return 0;
			}
		}
		int get_confine() { return confine; }

};

// workhourse function in fix_wall_body_polyhedron_agent.cpp
// returns the reaction force and torque of cell-surface interaction
my6Vec cell_surface_gforce(Cell* cell_1, double kn, double A);

// see if the cell is in contact with z=0 plane
int cell_contact(Cell* cell_1);

// see if the cell is verticalized
int cell_vertical(Cell* cell_1);

// returns the damping force and torque of cell-surface interaction
my6Vec compute_surface_damping_force(Cell* cell_1, double* v, double* omega, double nu_1);

// returns equivalent contact area density at relative position rr
double contact_area_density(Cell* cell_1, double rr);

// kernel function in fix_wall_body_polyhedron_agent.cpp
void contact_forces_new(int ibody, Cell* cell, double* v, double* omega, double** f, double** torque, double kn, double cn, double ct, int shift_flag);

// keep verticalized cells on the surface
void keep_verticalized_cell(int ibody, Cell *cell, double **f, double **torque);

// shift a vector by given flag
void shift_vector(double* v, int n, int flag);

// gaussian quadrature of n points
double gaussian_quadrature(double* f, double* weights, int n);