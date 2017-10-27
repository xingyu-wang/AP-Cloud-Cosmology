#include "mover.h"
#include "apcloud_macro.h"

//#define OUTPUT_TIME_MOVER

Mover::Mover(Physical_particle *physical_particle1, Computational_particle *computational_particle1)
: physical_particle(physical_particle1), computational_particle(computational_particle1)
{
octree=NULL;
}

Mover::~Mover()
{
octree->~Octree();
}

int Mover::select_computational_particles()
{
	int number=-1;
	int MAXNC=1000000;
#ifdef OUTPUT_TIME_MOVER
	double time1,time2;
#endif
	if(octree)
		octree->~Octree();
#ifdef OUTPUT_TIME_MOVER
  time1 = omp_get_wtime();
	time2=time1;
#endif
	octree = new Octree(physical_particle,1*Octree_Depth,x_min,x_max,y_min,y_max,z_min,z_max,1.0*BC);
#ifdef OUTPUT_TIME_MOVER
	time1= omp_get_wtime()-time1;
	std::cout << "Octree Build Time : " << time1 << std::endl;
	time1= omp_get_wtime();
#endif
	while((number<0)&&(maxnc<MAXNC))
	{
		number=octree->select_computational_particles_2to1_3d(computational_particle, 1.0*GFD_Constant, 2, maxnc, 0,1*Octree_Depth,BC);
		if(number<0)	maxnc=maxnc*2;
	}
	if(number<0)	printf("Not enough memory for computational particles.\n");
#ifdef OUTPUT_TIME_MOVER
	time1= omp_get_wtime()-time1;
	std::cout << "Computational Paticles Selection time : " << time1 << std::endl;
	time2= omp_get_wtime()-time2;
	std::cout << "Mover time : " << time2 << std::endl;	
#endif
	return number;
}

double Mover::move_velocity(double dt)
{
#ifdef OUTPUT_TIME_MOVER
	double time1;
  time1 = omp_get_wtime();
#endif
	octree->interpolation_from_node_to_particle_2to1_3d(computational_particle, physical_particle);
#ifdef OUTPUT_TIME_MOVER
	time1= omp_get_wtime()-time1;
	std::cout << "Interpolation Time : " << time1 << std::endl;
#endif
	int i;
	double MAX_CHANGE=0.0;
#ifdef OUTPUT_TIME_MOVER
  time1 = omp_get_wtime();
#endif
#pragma omp parallel
{
	double max_change=0.0,change;
	#pragma omp for
	for(i=0;i<physical_particle->size;i++)
	{
		physical_particle->v_x[i]=physical_particle->v_x[i]+physical_particle->a_x[i]*Constant_For_Mover_Velocity*dt;
		physical_particle->v_y[i]=physical_particle->v_y[i]+physical_particle->a_y[i]*Constant_For_Mover_Velocity*dt;
		physical_particle->v_z[i]=physical_particle->v_z[i]+physical_particle->a_z[i]*Constant_For_Mover_Velocity*dt;
		change=(physical_particle->a_x[i]*physical_particle->a_x[i]+physical_particle->a_y[i]*physical_particle->a_y[i]+physical_particle->a_z[i]*physical_particle->a_z[i])*Constant_For_Mover_Velocity*Constant_For_Mover_Velocity*dt*dt;
		if(change>max_change)
			max_change=change;
	}
	#pragma omp critical
	if(max_change>MAX_CHANGE)
		MAX_CHANGE=max_change;
}
#ifdef OUTPUT_TIME_MOVER
	time1= omp_get_wtime()-time1;
	std::cout << "Velocity Update Time : " << time1 << std::endl;
#endif
	return sqrt(MAX_CHANGE);
}

double Mover::move_position(double dt)
{
#ifdef OUTPUT_TIME_MOVER
	double time1;
  time1 = omp_get_wtime();
#endif
	int i;
	double MAX_CHANGE=0.0;
#pragma omp parallel
{
	double max_change=0.0,change;
	#pragma omp for
	for(i=0;i<physical_particle->size;i++)
	{
		physical_particle->x[i]=periodic_wrap(physical_particle->x[i]+physical_particle->v_x[i]*Constant_For_Mover_Position*dt,0);
		physical_particle->y[i]=periodic_wrap(physical_particle->y[i]+physical_particle->v_y[i]*Constant_For_Mover_Position*dt,1);
		physical_particle->z[i]=periodic_wrap(physical_particle->z[i]+physical_particle->v_z[i]*Constant_For_Mover_Position*dt,2);
		change=(physical_particle->v_x[i]*physical_particle->v_x[i]+physical_particle->v_y[i]*physical_particle->v_y[i]+physical_particle->v_z[i]*physical_particle->v_z[i])*Constant_For_Mover_Position*Constant_For_Mover_Position*dt*dt;
		if(change>max_change)
			max_change=change;
	}
	#pragma omp critical
	if(max_change>MAX_CHANGE)
		MAX_CHANGE=max_change;
}
#ifdef OUTPUT_TIME_MOVER
	time1= omp_get_wtime()-time1;
	std::cout << "Position Update Time : " << time1 << std::endl;
#endif
	return sqrt(MAX_CHANGE);
}

double Mover::periodic_wrap(double x,int dimension)
{
  if(BC==1)
  {
		if(dimension==0)//x
		{
   	while(x >= x_max)
     	x -= x_max-x_min;

   	while(x < x_min)
     	x += x_max-x_min;
		}
		if(dimension==1)//y
		{
   	while(x >= y_max)
     	x -= y_max-y_min;

   	while(x < y_min)
     	x += y_max-y_min;
		}
		if(dimension==2)//z
		{
   	while(x >= z_max)
     	x -= z_max-z_min;

   	while(x < z_min)
     	x += z_max-z_min;
		}
  }
  return x;
}
