#include <vector>


#ifndef PARTICLE
struct Physical_particle {

  int size;	//number of physical particles
  double *x;	//Position
  double *y;
  double *z;
  double *v_x;	//Velocity
  double *v_y;
  double *v_z;
  int *ID;
  double *a_x;	//Acceleration
  double *a_y;
  double *a_z;
  double *mass;	//Mass
	double uniform_mass;

};

struct Computational_particle {

  int size;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> CellLength;
  std::vector<int> IndexOfFirstPhysicalParticle;//Note this index is the particle index in octree, different from the ID of physical particle. This index will be used in interpolation from computational particle to physical particle.
  std::vector<int> NumberOfPhysicalParticles;
  std::vector<int> Type;
	std::vector<int>	FirstNeigh;//Location of first neighbor in neigh24
	std::vector<int>	NumberNeigh;//Number of neighbors	
	std::vector<int>	Neigh24;
  std::vector<double> Density;
  std::vector<double> Potential;
  std::vector<double> p_x;
  std::vector<double> p_y;
  std::vector<double> p_z;
  std::vector<double> p_xx;
  std::vector<double> p_yy;
  std::vector<double> p_zz;
  std::vector<double> p_xy;
  std::vector<double> p_xz;
  std::vector<double> p_yz;

	int nc;
	int nv;
	int nb;


};

int clear(Computational_particle *computational_particle);
int remove(Computational_particle *computational_particle);
#define PARTICLE
#endif
