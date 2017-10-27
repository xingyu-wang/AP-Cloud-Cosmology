#include "particle.h"

int clear(Computational_particle *computational_particle)
{
	computational_particle->x.clear();
	computational_particle->y.clear();
	computational_particle->z.clear();
	computational_particle->CellLength.clear();
	computational_particle->IndexOfFirstPhysicalParticle.clear();
	computational_particle->NumberOfPhysicalParticles.clear();
	computational_particle->Type.clear();
	computational_particle->FirstNeigh.clear();
	computational_particle->NumberNeigh.clear();
	computational_particle->Neigh24.clear();
	computational_particle->Density.clear();
	computational_particle->Potential.clear();
	computational_particle->p_x.clear();
	computational_particle->p_y.clear();
	computational_particle->p_z.clear();
	computational_particle->p_xx.clear();
	computational_particle->p_yy.clear();
	computational_particle->p_zz.clear();
	computational_particle->p_xy.clear();
	computational_particle->p_xz.clear();
	computational_particle->p_yz.clear();
}

int remove(Computational_particle *computational_particle)
{
	computational_particle->x.~vector();
	computational_particle->y.~vector();
	computational_particle->z.~vector();
	computational_particle->CellLength.~vector();
	computational_particle->IndexOfFirstPhysicalParticle.~vector();
	computational_particle->NumberOfPhysicalParticles.~vector();
	computational_particle->Type.~vector();
	computational_particle->FirstNeigh.~vector();
	computational_particle->NumberNeigh.~vector();
	computational_particle->Neigh24.~vector();
	computational_particle->Density.~vector();
	computational_particle->Potential.~vector();
	computational_particle->p_x.~vector();
	computational_particle->p_y.~vector();
	computational_particle->p_z.~vector();
	computational_particle->p_xx.~vector();
	computational_particle->p_yy.~vector();
	computational_particle->p_zz.~vector();
	computational_particle->p_xy.~vector();
	computational_particle->p_xz.~vector();
	computational_particle->p_yz.~vector();
}

