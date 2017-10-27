#include "apcloud_macro.h"
#include "controller.h"

double check_radius(Physical_particle *pp, double x0, double y0, double z0);

//#define _PARALLEL

Controller::Controller(Physical_particle *Physical_Particle1, Computational_particle *Computational_Particle1, double Start_A1, double End_A1, double Delta_A1, Solver *solver1, Viewer *viewer1, Mover *mover1, int Write_Interval1, int Boundary_Condition1, double x_min1, double x_max1, double y_min1, double y_max1, double z_min1, double z_max1, int Octree_Depth1, double GFD_Constant1, int lm1)
: physical_particle(Physical_Particle1), computational_particle(Computational_Particle1), solver(solver1), viewer(viewer1), mover(mover1), Start_A(Start_A1), End_A(End_A1), Delta_A(Delta_A1), WriteInterval(Write_Interval1), BC(Boundary_Condition1), x_min(x_min1), x_max(x_max1), y_min(y_min1), y_max(y_max1), z_min(z_min1), z_max(z_max1), Octree_Depth(Octree_Depth1), GFD_Constant(GFD_Constant1), lm(lm1)
{
}

double Controller::Get_Constant_For_Density_Estimator(double a)
{
#ifdef COSMOLOGY_UNIT
	return 8.0*a*a*a*PI_VALUE*GravitionalConstant/3.0/100000.0/100000.0*MPC_CM*MPC_CM/Omega*UnitMass_in_g/UnitLength_in_cm/UnitLength_in_cm/UnitLength_in_cm*MassPerParticle;
#else
	return 1.0;
#endif
}

double Controller::Get_Constant_For_Potential_Solver(double a)
{
#ifdef COSMOLOGY_UNIT
	return 3.0/2.0*Omega/a;
#else
	return 1000.0;
#endif
}

double Controller::Get_Constant_For_Field_Solver(double a)
{
	return -1.0;
}

double Controller::Get_Constant_For_Mover_Position(double a)
{
#ifdef COSMOLOGY_UNIT
	return 1.0/sqrt(Omega/a+OmegaLambda*a*a);
#else
	return 1.0;
#endif
}

double Controller::Get_Constant_For_Mover_Velocity(double a)
{
#ifdef COSMOLOGY_UNIT
	return 1.0/sqrt(Omega/a+OmegaLambda*a*a)/a/a;
#else
	return 1.0;
#endif
}

int Controller::solve()
{
#ifdef _PARALLEL
	int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
#else
	int rank=0;
#endif
	int errorflag,ite;

	int dynamic_time_step=0;
#ifdef ADAPTIVE_TIME_STEP
  dynamic_time_step=1;
#endif
	double dynamic_time_step_size=0.05*Delta_A;
	std::vector<double> dynamic_time;
//	std::vector<double> dynamic_radius;

	double t=Start_A;
	double max_pchange=1.0,max_vchange,delta_t;
	double t_Output=0.0,t_Computational=0.0,t_Solve=0.0,t_Move=0.0,t_Total=0.0,t_SO=0.0,t_SC=0.0,t_SS=0.0,t_SM=0.0,t_ST=0.0;
	mover->Set_BoundingBox(x_min,x_max,y_min,y_max,z_min,z_max);
	mover->Set_maxnc(10000);
	mover->Set_Octree_Depth(1*Octree_Depth);
	mover->Set_Constant_For_GFD_Constant(1.0*GFD_Constant);
	mover->Set_BC(1*BC);
	mover->Set_Constant_For_Mover_Position(Get_Constant_For_Mover_Position(t));
	mover->Set_Constant_For_Mover_Velocity(Get_Constant_For_Mover_Velocity(t));
	if(rank==0)
	{
		#ifdef COSMOLOGY_UNIT
			printf("System of units in Gadget is used.\n");
//			printf("%e\n",8.0*t*t*t*PI_VALUE*GRAVITATIONAL_CONSTANT/3.0/100000.0/100000.0/Omega*UnitMass_in_g/UnitLength_in_cm/UnitLength_in_cm/UnitLength_in_cm*MassPerParticle);
		#endif
		viewer->View_Physical_Particles_Limit(t, x_min, y_min, z_min, x_max, y_max, z_max);
		errorflag=mover->select_computational_particles();
	}

//	return 0;

	MPI_Bcast(&(computational_particle->size), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&(computational_particle->nc), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&(computational_particle->nv), 1, MPI_INT, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&(computational_particle->nb), 1, MPI_INT, 0, PETSC_COMM_WORLD);

	solver->Set_BoundingBox(x_min,x_max,y_min,y_max,z_min,z_max);
	solver->Set_Boundary_Condition(BC);
	solver->Set_NumberOfGFDNeighbors(lm);
	solver->Set_Constant_For_Density_Estimator(Get_Constant_For_Density_Estimator(t));
	solver->Set_Constant_For_Potential_Solver(Get_Constant_For_Potential_Solver(t));
	solver->Set_Constant_For_Field_Solver(Get_Constant_For_Field_Solver(t));

	if(rank==0)
		max_pchange=mover->move_position(0.5*Delta_A);
	if(rank==0)
		printf("Maximum position change from %f to %f is: %f\n",t,t+0.5*Delta_A,max_pchange);

	t+=0.5*Delta_A;

	ite=0;
	delta_t=Delta_A;
	

	while(t<End_A)
	{
		if(rank==0)
		{
			t_Total=omp_get_wtime();
			if(ite%WriteInterval==0)
			{
				t_Output=omp_get_wtime();
				printf("Output result... ");
				viewer->View_Physical_Particles(t);
				viewer->View_Computational_Particles(t);
				t_Output=omp_get_wtime()-t_Output;
				t_SO=t_SO+t_Output;
				printf("time=%f\n",t_Output);	
			}

			mover->Set_Constant_For_Mover_Velocity(Get_Constant_For_Mover_Velocity(t));
			solver->Set_Constant_For_Density_Estimator(Get_Constant_For_Density_Estimator(t));
			solver->Set_Constant_For_Potential_Solver(Get_Constant_For_Potential_Solver(t));
			solver->Set_Constant_For_Field_Solver(Get_Constant_For_Field_Solver(t));

			printf("Build octree and select computational particles... ");
			t_Computational=omp_get_wtime();
			errorflag=mover->select_computational_particles();
			t_Computational=omp_get_wtime()-t_Computational;
			t_SC=t_SC+t_Computational;
			printf("time=%f\n",t_Computational);
			printf("Number of computational particles=%d\n",errorflag);

			if(dynamic_time_step)
			{
//				delta_t=dynamic_time_step_size*delta_t/max_pchange;

				double max_density=0.0;
				for(int i=0;i<computational_particle->size;i++)
				{
					double temp_density=computational_particle->NumberOfPhysicalParticles[i]/computational_particle->CellLength[i]/computational_particle->CellLength[i]/computational_particle->CellLength[i];
					if(max_density<temp_density)
						max_density=temp_density;
				}
				double min_interparticle_distance=1.0/pow(max_density,1.0/2);
				delta_t=dynamic_time_step_size*min_interparticle_distance;

				if(delta_t>Delta_A)
					delta_t=Delta_A;
			}
			else
				delta_t=Delta_A;
			printf("Time step size=%f\n",delta_t);

			if((t+delta_t)>End_A)
				delta_t=End_A-t;
			printf("Estimate density, solve potential and gravity field... \n");
			t_Solve=omp_get_wtime();


		}

		MPI_Bcast(&(computational_particle->size), 1, MPI_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&(computational_particle->nc), 1, MPI_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&(computational_particle->nv), 1, MPI_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&(computational_particle->nb), 1, MPI_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&t, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&delta_t, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

		mover->Set_Constant_For_Mover_Velocity(Get_Constant_For_Mover_Velocity(t));
		mover->Set_Constant_For_Mover_Position(Get_Constant_For_Mover_Position(t+delta_t/2));
		solver->Set_Constant_For_Density_Estimator(Get_Constant_For_Density_Estimator(t));
		solver->Set_Constant_For_Potential_Solver(Get_Constant_For_Potential_Solver(t));
		solver->Set_Constant_For_Field_Solver(Get_Constant_For_Field_Solver(t));

		errorflag=solver->solve();
		if(rank==0)
		{
			t_Solve=omp_get_wtime()-t_Solve;
			t_SS=t_SS+t_Solve;
			printf("time=%f\n",t_Solve);		
			printf("Update velocity and position of particles... ");
			t_Move=omp_get_wtime();
			max_vchange=mover->move_velocity(delta_t);	
			max_pchange=mover->move_position(delta_t);
			t_Move=omp_get_wtime()-t_Move;
			t_SM=t_SM+t_Move;
			printf("time=%f\n",t_Move);	
			printf("Maximum velocity change from %f to %f is: %e\n",t,t+delta_t,max_vchange);	
			printf("Maximum position change from %f to %f is: %e\n",t,t+delta_t,max_pchange);

			t_Total=omp_get_wtime()-t_Total;
			printf("Running time for %dth step=%f\n\n",ite,t_Total);	

			t_ST=t_ST+t_Total;
			t+=delta_t;

//			printf("\nTime Radius\n");
//			printf("%f %f\n",t,check_radius(physical_particle, 100.0, 200.0, 300.0));

			dynamic_time.push_back(t);
//			dynamic_radius.push_back(check_radius(physical_particle, 100.0, 200.0, 300.0));

			ite=ite+1;
		}
		MPI_Bcast(&t, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
	}
	if(rank==0)
	{
		viewer->View_Physical_Particles(t);
		viewer->View_Computational_Particles(t);

		int np;
		MPI_Comm_size(PETSC_COMM_WORLD,&np);
		std::cout << "NP : " << np <<std::endl;	
		std::cout << "OMP_NUM_THREADS : " << omp_get_max_threads() << std::endl;
		printf("Total running time: %f\n",t_ST);
		printf("Running time for output: %f\n",t_SO);
		printf("Running time for computational particles select: %f\n",t_SC);
		printf("Running time for density estimation, potential solver and differentiation: %f\n",t_SS);
		printf("Running time for particle position and velocity update: %f\n",t_SM);

		for(int i=0;i<dynamic_time.size();i++)
			printf("%f\n",dynamic_time[i]);
//			printf("%f %f\n",dynamic_time[i],dynamic_radius[i]);
	}

	return errorflag;
}

double check_radius(Physical_particle *pp, double x0, double y0, double z0)
{
	double max_radius=0,r;
	for (int i=0; i<pp->size; i++)
	{
		r=(pp->x[i]-x0)*(pp->x[i]-x0)+(pp->y[i]-y0)*(pp->y[i]-y0)+(pp->z[i]-z0)*(pp->z[i]-z0);
		if (r>max_radius)
			max_radius=r;
	}
	max_radius=sqrt(max_radius);
	return max_radius;
}
