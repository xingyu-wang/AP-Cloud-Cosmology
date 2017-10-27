#include <vector>
#include <ctime>
#include <math.h>
#include <omp.h>
#include "solver.h"
#include "viewer.h"
#include "mover.h"
#include "parameter.h"
#include "particle.h"

#ifndef CONTROLLER
class Controller {

private:

	Physical_particle *physical_particle;
	Computational_particle *computational_particle;
	Solver *solver;
	Viewer *viewer;
	Mover *mover;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	int BC;
	int WriteInterval;
	int Octree_Depth;
	int lm;
	double GFD_Constant;
	double Start_A;
	double End_A;
	double Delta_A;
	double Current_A;

	double Constant_For_Density_Estimator;
	double Constant_For_Potential_Solver;
	double Constant_For_Field_Solver;
	double Constant_For_Mover_Position;
	double Constant_For_Mover_Velocity;

private:
	double Get_Constant_For_Density_Estimator(double t);
	double Get_Constant_For_Potential_Solver(double t);
	double Get_Constant_For_Field_Solver(double t);
	double Get_Constant_For_Mover_Position(double t);
	double Get_Constant_For_Mover_Velocity(double t);

public:
	Controller(Physical_particle *physical_particle1, Computational_particle *computational_particle1, double Start_A1, double End_A1, double Delta_A1, Solver *solver1, Viewer *viewer1, Mover *mover1, int Write_Interval1, int Boundary_Condition1, double x_min1, double x_max1, double y_min1, double y_max1, double z_min1, double z_max1, int Octree_Depth1, double GFD_Constant1, int lm1);

public:
	virtual int solve();

};
#define CONTROLLER
#endif
