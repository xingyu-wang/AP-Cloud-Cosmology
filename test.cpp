#include <iostream>
#include <vector>
#include "petscksp.h"
#include "parameter.h"
#include "particle.h"
#include "initialization.h"
#include "controller.h"
#include "solver.h"
#include "viewer.h"
#include "mover.h"
#include "apcloud_macro.h"
#include <omp.h>

//#define _PARALLEL

int main(int argc, char *argv[])
{
  if(argc < 2)
    {
	  fprintf(stdout, "\nParameters are missing.\n");
	  fprintf(stdout, "Call with <ParameterFile>\n\n");
    }



PetscInitialize(&argc,&argv,(char*)0,NULL);

#ifdef _PARALLEL
	int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);	
#else
	int rank=0;
#endif

	if(rank==0)
	{
#ifdef _OPENMP
  	std::cout << "OMP_NUM_THREADS : " << omp_get_max_threads() << std::endl;
#endif
		fprintf(stdout, "\nReading parameters...\n");
	}

	read_parameterfile(argv[1]);

	if(rank==0)
		fprintf(stdout, "\nDone with read parameters.\n");

Physical_particle physical_particle;

if(rank==0)
{
	fprintf(stdout, "\nReading initial conditions...\n");
	initialization(&physical_particle);
	fprintf(stdout, "\nDone with read initial conditions.\n");
	printf("Number of physical particles=%d\n",physical_particle.size);
}

Computational_particle computational_particle;

Solver solver(&computational_particle, argc, argv);

Mover mover(&physical_particle,&computational_particle);

Viewer viewer(&physical_particle,&computational_particle);

Controller controller(&physical_particle,&computational_particle,1.0/(1+Redshift),1.0/(1+EndRedShift), 1.0*DeltaA, &solver, &viewer, &mover, 10, 1*Periodic, 1.0*BoundingBox_min_x, 1.0*BoundingBox_max_x, 1.0*BoundingBox_min_y, 1.0*BoundingBox_max_y, 1.0*BoundingBox_min_z, 1.0*BoundingBox_max_z, 1*MaxOctreeLevel, 1.0*CellSizeParameter, 1*NumberOfNeighbors);

MPI_Barrier(MPI_COMM_WORLD);

controller.solve();

MPI_Barrier(MPI_COMM_WORLD);

if(rank==0)
{
	printf("\nExecution Finished!\n");
}

PetscFinalize();
}
