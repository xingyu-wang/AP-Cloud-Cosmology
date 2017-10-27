#include <vector>
#include <lapacke.h>
#include <cmath>
#include <iostream>
#include <ctime>
#include <omp.h>
#include "particle.h"
#include "parameter.h"
#include "petscksp.h"

#ifndef SOLVER

class Solver {

private:
	Computational_particle *cp;
	double Constant_For_Density_Estimator;
	double Constant_For_Potential_Solver;
	double Constant_For_Field_Solver;
	int BC;
	int lm;
	int ln;
	int dirichlet_point;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;

  int argc;
  char **args;

	Mat A1,A2;
	Vec rhs,rho,phi,bv,vgather;
  VecScatter ctx;

	KSP ksp1,ksp2;
	PC pc1,pc2;
	PetscScalar *temp_array;
		

	int Build_Linear_System();
	int Solve_Density();
	int Solve_Potential();
	int Differentiate();
	int CleanPetscData();
	double periodic_wrap(double x, int dimension);
	double get_squared_distance(double dx,double dy,double dz);
	int searchneighboroctant(int i, std::vector<double> *distance, std::vector<int> *neigh);

public:
	Solver(Computational_particle *computational_particle1,int argc_o, char *args_o[]);
	~Solver();
	int Set_Boundary_Condition(int bc){BC=bc;return 0;}
	int Set_NumberOfGFDNeighbors(int n){lm=n;return 0;}
	int Set_Constant_For_Density_Estimator(double c){Constant_For_Density_Estimator=c;return 0;}
	int Set_Constant_For_Potential_Solver(double c){Constant_For_Potential_Solver=c;return 0;}
	int Set_Constant_For_Field_Solver(double c){Constant_For_Field_Solver=c;return 0;}
	int Set_BoundingBox(double x1, double x2, double y1, double y2, double z1, double z2){x_min=x1;x_max=x2;y_min=y1;y_max=y2;z_min=z1;z_max=z2;return 0;}

	int solve();

};

void aquickSort(std::vector<double> &x, std::vector<int> &y, int l, int r);
int apartition(std::vector<double> &a,std::vector<int> &b, int l, int r);

#define SOLVER
#endif
