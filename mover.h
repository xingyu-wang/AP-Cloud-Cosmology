#include "particle.h"
#include "parameter.h"
#include "octree.h"
#include <math.h>
#include <omp.h>
#include <iostream>

#ifndef MOVER

class Mover {
private:
	Physical_particle *physical_particle;
	Computational_particle *computational_particle;
	Octree *octree;
	int BC;
	int Octree_Depth;
	double GFD_Constant;
	double Constant_For_Mover_Position;
	double Constant_For_Mover_Velocity;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	int maxnc;
	double periodic_wrap(double x,int dimension);


public:
	Mover(Physical_particle *physical_particle1, Computational_particle *computational_particle1);
	~Mover();

	int Set_Constant_For_Mover_Position(double c){Constant_For_Mover_Position=c;return 0;}
	int Set_Constant_For_Mover_Velocity(double c){Constant_For_Mover_Velocity=c;return 0;}
	int Set_Constant_For_GFD_Constant(double c){GFD_Constant=c;return 0;}
	int Set_Octree_Depth(int c){Octree_Depth=c;return 0;}	
	int Set_BC(int c){BC=c;return 0;}	
	int Set_BoundingBox(double x1, double x2, double y1, double y2, double z1, double z2){x_min=x1;x_max=x2;y_min=y1;y_max=y2;z_min=z1;z_max=z2;return 0;};
	int Set_maxnc(int c){maxnc=c;return 0;};
	int select_computational_particles();
	double move_position(double dt);
	double move_velocity(double dt);

};

#define MOVER
#endif
