#include "solver.h"
#include "apcloud_macro.h"

//#define _PARALLEL

//#define OUTPUT_TIME_SOLVER

Solver::Solver(Computational_particle *computational_particle1,int argc_o, char *args_o[])
: cp(computational_particle1),argc(argc_o),args(args_o)
{
BC=1;
Constant_For_Density_Estimator=1.0;
Constant_For_Potential_Solver=1.0;
Constant_For_Field_Solver=1.0;
lm=17;
ln=9;
dirichlet_point=1;
}

Solver::~Solver()
{
}

int Solver::Build_Linear_System()
{
	PetscErrorCode ierr;
	int i;

#ifdef _PARALLEL
	int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);	
#else
	int rank=0;
#endif

	ierr = VecCreate(PETSC_COMM_WORLD,&rhs);CHKERRQ(ierr);
	ierr = VecSetSizes(rhs,PETSC_DECIDE,cp->nc+cp->nv);CHKERRQ(ierr);
	ierr = VecSetFromOptions(rhs);CHKERRQ(ierr);
	ierr = VecDuplicate(rhs,&rho);CHKERRQ(ierr);
	ierr = VecDuplicate(rhs,&phi);CHKERRQ(ierr);
	if(BC-1)
	{
		double p1;
		ierr = VecCreate(PETSC_COMM_WORLD,&bv);CHKERRQ(ierr);
		ierr = VecSetSizes(bv,PETSC_DECIDE,cp->nb);CHKERRQ(ierr);
		ierr = VecSetFromOptions(bv);CHKERRQ(ierr);
		p1=0.0;
		for(i=0;i<cp->nb;i++)
			ierr=VecSetValues(bv,1,&i,&p1,INSERT_VALUES);CHKERRQ(ierr);
	}

	ierr = MatCreate(PETSC_COMM_WORLD,&A1);CHKERRQ(ierr);
	ierr = MatSetSizes(A1,PETSC_DECIDE,PETSC_DECIDE,cp->nc+cp->nv,cp->nc+cp->nv);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A1);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A1,lm+1,NULL,lm+1,NULL);CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(A1,lm+1,NULL);CHKERRQ(ierr);
	ierr = MatSetUp(A1);CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&A2);CHKERRQ(ierr);
	ierr = MatSetSizes(A2,PETSC_DECIDE,PETSC_DECIDE,cp->nc+cp->nv,cp->nc+cp->nv);CHKERRQ(ierr);
	ierr = MatSetFromOptions(A2);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A2,lm+1,NULL,lm+1,NULL);CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(A2,lm+1,NULL);CHKERRQ(ierr);
	ierr = MatSetUp(A2);CHKERRQ(ierr);

	if(rank==0)
	{
		int j,k,info;
		double w,p1,p2,p3,p4,maxdis,gamma=0.7,weight;
		std::vector<int> neigh(16*lm);
		std::vector<double> distance(16*lm);
		double B[lm*ln],b[lm*lm];
		for (i=0;i<(cp->nc+cp->nv);i++)
		{
			searchneighboroctant(i, &distance, &neigh);

			maxdis=0;
			for (j=0;j<lm;j++)
				if((distance)[j]>maxdis)
					maxdis=(distance)[j];
			maxdis=maxdis*gamma;
			for (j=0;j<lm;j++)
			{
					/*exponential weight fucntion*/
				weight=(exp(-(distance)[j]*(distance)[j]/maxdis/maxdis)-exp(-4))/(1-exp(-4));
		
				p1=periodic_wrap(cp->x[(neigh)[j]]-cp->x[i],0)/cp->CellLength[i];/*hc[i] is a pre conditioning factor*/
				p2=periodic_wrap(cp->y[(neigh)[j]]-cp->y[i],1)/cp->CellLength[i];
				p3=periodic_wrap(cp->z[(neigh)[j]]-cp->z[i],2)/cp->CellLength[i];
				B[j*ln]=p1;
				B[j*ln+1]=p2;
				B[j*ln+2]=p3;
				B[j*ln+3]=p1*p1/2.0;
				B[j*ln+4]=p2*p2/2.0;
				B[j*ln+5]=p3*p3/2.0;
				B[j*ln+6]=p1*p2;
				B[j*ln+7]=p1*p3;
				B[j*ln+8]=p2*p3;
				for(k=0;k<ln;k++)
					B[j*ln+k]=B[j*ln+k]*weight;
				for(k=0;k<lm;k++)
					b[j*lm+k]=0.0;
				b[j*lm+j]=weight;
			}
			/*	Solving LSQ by QR	*/
			info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N', lm, ln, lm, B, ln,b, lm );
			if( info > 0 ) {
				      printf( "The diagonal element %i of the triangular factor ", info );
				      printf( "of B is zero, so that B does not have full rank;\n" );
				      printf( "the least squares solution could not be computed.\n" );
				printf("%d %f %f %f\n",i,cp->x[i],cp->y[i],cp->z[i]);
				for(j=0;j<lm;j++)
					printf("%d %f %f %f\n",(neigh)[j],cp->x[(neigh)[j]],cp->y[(neigh)[j]],cp->z[(neigh)[j]]);
				      exit( 1 );
			}		
			p2=0.0;
			for (j=0;j<lm;j++){
				p1=(b[3*lm+j]+b[4*lm+j]+b[5*lm+j])/24.0;

				p2=p2-p1;
				if((neigh[j])<(cp->nc+cp->nv))
					ierr = MatSetValues(A1,1,&i,1,&(neigh[j]),&p1,INSERT_VALUES);CHKERRQ(ierr);
			}

			p2=p2+1;
			ierr = MatSetValues(A1,1,&i,1,&i,&p2,INSERT_VALUES);CHKERRQ(ierr);		
			/* Right hand side given by Monte Carlo integration*/
			p1=1.0*cp->NumberOfPhysicalParticles[i]/cp->CellLength[i]/cp->CellLength[i]/cp->CellLength[i]*Constant_For_Density_Estimator;

	//	p1=double_np[i]/n/hc[i]/hc[i];
	/*	p4=0;
			for(int j1=-20;j1<=20;j1++)
				for(int j2=-20;j2<=20;j2++)
					for(int j3=-20;j3<20;j3++)
					{
						p1=xc[i]+j1*hc[i]/40;
						p2=yc[i]+j2*hc[i]/40;
						p3=zc[i]+j2*hc[i]/40;
						w=sqrt(p1*p1+p2*p2+p3*p3);
						p1=a1*exp(-w*w/sigma1/sigma1/2.0)+a2*exp(-(w-r0)*(w-r0)/sigma2/sigma2/2.0);
						p1=p1/scale;
						p4=p4+p1;
				}
			p1=p4/41/41/41;*/
			ierr = VecSetValues(rhs,1,&i,&p1,INSERT_VALUES);CHKERRQ(ierr);
	/* Matrix for phi*/
			if((i==0)&&(BC==1)&&dirichlet_point)//Dirichlet point
			{
				p1=1.0;
				p3=0.0;
				ierr = MatSetValues(A2,1,&i,1,&i,&p1,INSERT_VALUES);CHKERRQ(ierr);
			}	
			else
			{
				p2=0.0;
				p3=0.0;
				for (j=0;j<lm;j++){
					p1=(b[3*lm+j]+b[4*lm+j]+b[5*lm+j])/cp->CellLength[i]/cp->CellLength[i];
					p2=p2-p1;
					if(neigh[j]>cp->nc+cp->nv)	printf("%d %f %f %f %d\n",i,cp->x[i],cp->y[i],cp->z[i],neigh[j]);
					if(cp->Type[neigh[j]]<1){
						ierr = MatSetValues(A2,1,&i,1,&(neigh[j]),&p1,INSERT_VALUES);CHKERRQ(ierr);
					}
					else{/*Move Dirichlet boundary condition to the right hand side.*/
						k=neigh[j]-cp->nc-cp->nv;
						ierr = VecGetValues(bv, 1, &k, &p4);CHKERRQ(ierr);
						p3=p3-p4*p1;
					}
				}
				ierr = MatSetValues(A2,1,&i,1,&i,&p2,INSERT_VALUES);CHKERRQ(ierr);
			}
			cp->p_x[i]=p3;	
		}
	}
	ierr = MatAssemblyBegin(A1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A1,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(A2,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A2,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	return 0;
}

int Solver::Solve_Density()
{
	PetscErrorCode ierr;
	int i;

#ifdef _PARALLEL
	int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);	
#else
	int rank=0;
#endif

  VecAssemblyBegin(rhs);
  VecAssemblyEnd(rhs);
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp1);CHKERRQ(ierr);
  ierr = KSPSetType(ksp1,KSPGMRES);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp1,A1,A1);CHKERRQ(ierr);

	ierr = KSPGetPC(ksp1,&pc1);CHKERRQ(ierr);
  ierr = PCSetType(pc1,PCJACOBI);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp1);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp1);CHKERRQ(ierr);
	ierr = KSPSolve(ksp1,rhs,rho);CHKERRQ(ierr);

  ierr = VecScatterCreateToZero(rho,&ctx,&vgather);CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx,rho,vgather,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx,rho,vgather,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

	if(rank==0)
	{
		ierr = VecGetArray(vgather,&temp_array);CHKERRQ(ierr);
		for(i=0;i<cp->size;i++)
			cp->Density[i]=temp_array[i];	
		ierr = VecRestoreArray(vgather,&temp_array);CHKERRQ(ierr);
	}	
	return 0;
}

int Solver::Solve_Potential()
{
	PetscErrorCode ierr;
	int i;
//	MatNullSpace   nullsp;

#ifdef _PARALLEL
	int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);	
#else
	int rank=0;
#endif

	if(rank==0)
	{
		double p1,p2,p3;

		ierr = VecGetArray(vgather,&temp_array);CHKERRQ(ierr);
		if(BC==0)
			for (i=0;i<cp->nc+cp->nv;i++){
				p1=temp_array[i];
				p1=(p1+cp->p_x[i])*Constant_For_Potential_Solver;	
				ierr = VecSetValues(rhs,1,&i,&p1,INSERT_VALUES);CHKERRQ(ierr);
			}
		else
		{
			p2=0;
			p3=0;
			for (i=0;i<cp->nc+cp->nv;i++){
				p1=temp_array[i];
				p2=p2+p1*cp->CellLength[i]*cp->CellLength[i]*cp->CellLength[i];
				p3=p3+cp->CellLength[i]*cp->CellLength[i]*cp->CellLength[i];	
			}
			p2=p2/p3;//density average
			for (i=0;i<cp->nc+cp->nv;i++){
				p1=temp_array[i];
				p1=(p1-p2)*Constant_For_Potential_Solver;	
				ierr = VecSetValues(rhs,1,&i,&p1,INSERT_VALUES);CHKERRQ(ierr);
			}
		}
		ierr = VecRestoreArray(vgather,&temp_array);CHKERRQ(ierr);
	}
  VecAssemblyBegin(rhs);
  VecAssemblyEnd(rhs);

	MPI_Barrier(MPI_COMM_WORLD);

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp2);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp2,A2,A2);CHKERRQ(ierr);
  ierr = KSPSetType(ksp2,KSPGMRES);CHKERRQ(ierr);
//	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nullsp);CHKERRQ(ierr);
//	ierr = MatSetNullSpace(A2, nullsp);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp2,&pc2);CHKERRQ(ierr);
  ierr = PCSetType(pc2,PCJACOBI);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp2);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp2);CHKERRQ(ierr);
	ierr = KSPSolve(ksp2,rhs,phi);CHKERRQ(ierr);

//	ierr = MatNullSpaceRemove(nullsp,phi);CHKERRQ(ierr);
//	ierr = MatNullSpaceDestroy(&nullsp);CHKERRQ(ierr);

  ierr = VecScatterCreateToZero(phi,&ctx,&vgather);CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx,phi,vgather,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx,phi,vgather,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

	if(rank==0)
	{
		ierr = VecGetArray(vgather,&temp_array);CHKERRQ(ierr);
		for(i=0;i<cp->size;i++)
			cp->Potential[i]=temp_array[i];	
		ierr = VecRestoreArray(vgather,&temp_array);CHKERRQ(ierr);
	}	
	return 0;
}

int Solver::Differentiate()
{
	PetscErrorCode ierr;
	int i,j,k,info;
	double maxdis,gamma=0.7,weight;
	double p1,p2,p3,p4;
	double B[lm*ln],b[lm*lm];
	std::vector<int> neigh(16*lm);
	std::vector<double> distance(16*lm);
	ierr = VecGetArray(vgather,&temp_array);CHKERRQ(ierr);
	for(i=0;i<cp->nc+cp->nv;i++)
	{
		searchneighboroctant(i,&distance,&neigh);
		maxdis=0;
		for (j=0;j<lm;j++)
			if(distance[j]>maxdis)
				maxdis=distance[j];
		maxdis=maxdis*gamma;
			p4=temp_array[i];
		for (j=0;j<lm;j++)
		{
			weight=(exp(-distance[j]*distance[j]/maxdis/maxdis)-exp(-4))/(1-exp(-4));
			p1=periodic_wrap(cp->x[neigh[j]]-cp->x[i],0)/cp->CellLength[i];/*hc[i] is a pre conditioning factor*/
			p2=periodic_wrap(cp->y[neigh[j]]-cp->y[i],1)/cp->CellLength[i];
			p3=periodic_wrap(cp->z[neigh[j]]-cp->z[i],2)/cp->CellLength[i];
			B[j*ln]=p1;
			B[j*ln+1]=p2;
			B[j*ln+2]=p3;
			B[j*ln+3]=p1*p1/2.0;
			B[j*ln+4]=p2*p2/2.0;
			B[j*ln+5]=p3*p3/2.0;
			B[j*ln+6]=p1*p2;
			B[j*ln+7]=p1*p3;
			B[j*ln+8]=p2*p3;

			for(k=0;k<ln;k++)
				B[j*ln+k]=B[j*ln+k]*weight;
			k=neigh[j];
			if(k<(cp->nc+cp->nv)){
				p1=temp_array[k];
			}
			else{
				k=k-cp->nc-cp->nv;
				ierr = VecGetValues(bv, 1, &k, &p1);CHKERRQ(ierr);
			}
			b[j]=weight*(p1-p4)*Constant_For_Field_Solver;
		}
		info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N', lm, ln, 1, B, ln,b, 1 );
		if( info > 0 ) {
		        printf( "The diagonal element %i of the triangular factor ", info );
		        printf( "of B is zero, so that B does not have full rank;\n" );
		        printf( "the least squares solution could not be computed.\n" );
		        exit( 1 );
		}
		cp->p_x[i]=b[0]/cp->CellLength[i];/*phi_x*/
		cp->p_y[i]=b[1]/cp->CellLength[i];/*phi_y*/
		cp->p_z[i]=b[2]/cp->CellLength[i];/*phi_z*/

		cp->p_xx[i]=b[3]/cp->CellLength[i]/cp->CellLength[i];/*phi_xx*/
		cp->p_yy[i]=b[4]/cp->CellLength[i]/cp->CellLength[i];/*phi_yy*/
		cp->p_zz[i]=b[5]/cp->CellLength[i]/cp->CellLength[i];/*phi_zz*/
		cp->p_xy[i]=b[6]/cp->CellLength[i]/cp->CellLength[i];/*phi_xy*/
		cp->p_xz[i]=b[7]/cp->CellLength[i]/cp->CellLength[i];/*phi_xz*/
		cp->p_yz[i]=b[8]/cp->CellLength[i]/cp->CellLength[i];/*phi_yz*/	
	}
	ierr = VecRestoreArray(vgather,&temp_array);CHKERRQ(ierr);	
	return 0;
}

int Solver::CleanPetscData()
{
//	PetscFinalize();
	PetscErrorCode ierr;
	ierr = MatDestroy(&A1);CHKERRQ(ierr);
	ierr = MatDestroy(&A2);CHKERRQ(ierr);
	ierr = VecDestroy(&rhs);CHKERRQ(ierr);
	ierr = VecDestroy(&phi);CHKERRQ(ierr);
	ierr = VecDestroy(&rho);CHKERRQ(ierr);
	ierr = VecDestroy(&vgather);CHKERRQ(ierr);
	if(BC==0)
		ierr = VecDestroy(&bv);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp1);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp2);CHKERRQ(ierr);
  ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
	return 0;
}

int Solver::solve()
{
//  PetscInitialize(&argc,&args,(char*)0,NULL);
	int errorflag=0;
#ifdef OUTPUT_TIME_SOLVER
	double time;
#endif
#ifdef _PARALLEL
	int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);	
#else
	int rank=0;
#endif
	if(rank==0)
	{
#ifdef OUTPUT_TIME_SOLVER
		time=omp_get_wtime();
		printf("Building linear system...\n");
#endif
	}
		errorflag=Build_Linear_System();
	if(rank==0)
	{
#ifdef OUTPUT_TIME_SOLVER
		time=omp_get_wtime()-time;
		printf("Time to build linear system is %f\n",time);
#endif
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(errorflag)	return errorflag;
//	printf("Linear system build done!\n");
#ifdef OUTPUT_TIME_SOLVER
	time=omp_get_wtime();
#endif
	errorflag=Solve_Density();
#ifdef OUTPUT_TIME_SOLVER
	time=omp_get_wtime()-time;
	if(rank==0)
	{
		printf("Time to solve linear system for rho is %f\n",time);
	}
#endif
	MPI_Barrier(MPI_COMM_WORLD);
	if(errorflag)	return errorflag;
#ifdef OUTPUT_TIME_SOLVER
//	printf("Density solved!\n");
	time=omp_get_wtime();
#endif
	errorflag=Solve_Potential();
#ifdef OUTPUT_TIME_SOLVER
	time=omp_get_wtime()-time;
#endif
	if(rank==0)
	{
#ifdef OUTPUT_TIME_SOLVER
		printf("Time to solve linear system for phi is %f\n",time);
#endif
		if(errorflag)	return errorflag;
#ifdef OUTPUT_TIME_SOLVER
	//	printf("Potential solved!\n");
		time=omp_get_wtime();
#endif
		errorflag=Differentiate();
#ifdef OUTPUT_TIME_SOLVER
		time=omp_get_wtime()-time;
		printf("Time to differentiate is %f\n",time);
#endif
		if(errorflag)	return errorflag;
	}
//	printf("Differentiation done!\n");	
	CleanPetscData();
	return errorflag;
}

int Solver::searchneighboroctant(int i0, std::vector<double> *distance, std::vector<int> *neigh)
{
	int i,j,k,l,flag,nq,m,ite;
	double penalty;
	int nn[8];
	nq=10*lm;
	std::vector<int> temp_neigh(8*nq);
	std::vector<double> temp_dis(8*nq);
	double alpha=0.01,beta=0.001,gamma=0.0001;
	double A[9];
	int oct;
	A[0]=cos(alpha)*cos(gamma)-sin(alpha)*sin(beta)*sin(gamma);
	A[1]=-sin(alpha)*cos(beta);
	A[2]=-cos(alpha)*sin(gamma)-sin(alpha)*sin(beta)*cos(gamma);
	A[3]=cos(alpha)*sin(beta)*sin(gamma)+sin(alpha)*cos(gamma);
	A[4]=cos(alpha)*cos(beta);
	A[5]=cos(alpha)*sin(beta)*cos(gamma)-sin(alpha)*sin(gamma);
	A[6]=cos(beta)*sin(gamma);
	A[7]=-sin(beta);
	A[8]=cos(beta)*cos(gamma);

	double xt,yt,zt;
	for(i=0;i<8;i++)
		nn[i]=0;

	for(i=cp->FirstNeigh[i0];i<cp->FirstNeigh[i0]+cp->NumberNeigh[i0];i++)
	{
		m=cp->Neigh24[i];
		xt=A[0]*periodic_wrap(cp->x[m]-cp->x[i0],0)+A[1]*periodic_wrap(cp->y[m]-cp->y[i0],1)+A[2]*periodic_wrap(cp->z[m]-cp->z[i0],2);
		yt=A[3]*periodic_wrap(cp->x[m]-cp->x[i0],0)+A[4]*periodic_wrap(cp->y[m]-cp->y[i0],1)+A[5]*periodic_wrap(cp->z[m]-cp->z[i0],2);
		zt=A[6]*periodic_wrap(cp->x[m]-cp->x[i0],0)+A[7]*periodic_wrap(cp->y[m]-cp->y[i0],1)+A[8]*periodic_wrap(cp->z[m]-cp->z[i0],2);

		oct=(xt>0)+2*(yt>0)+4*(zt>0);
		temp_neigh[oct*nq+nn[oct]]=m;
		temp_dis[oct*nq+nn[oct]]=get_squared_distance(cp->x[m]-cp->x[i0],cp->y[m]-cp->y[i0],cp->z[m]-cp->z[i0]);
		nn[oct]=nn[oct]+1;
	}
	
	ite=0;
	for(l=nn[0]+nn[1]+nn[2]+nn[3]+nn[4]+nn[5]+nn[6]+nn[7];l<lm;l=nn[0]+nn[1]+nn[2]+nn[3]+nn[4]+nn[5]+nn[6]+nn[7])
	{
		k=0;
		for(i=0;i<8;i++)
		{
			for(j=0;j<nn[i];j++)
			{
				(*neigh)[k]=temp_neigh[i*nq+j];
				k++;
			}
		}
		for(i=0;i<l;i++)
		{
			for(j=cp->FirstNeigh[(*neigh)[i]];j<cp->FirstNeigh[(*neigh)[i]]+cp->NumberNeigh[(*neigh)[i]];j++)
			{
				m=cp->Neigh24[j];
				if(m>cp->nc+cp->nv+cp->nb)
				{
					printf("Invalid neighbour!\n");
					printf("%d %f %f %f %d %d %d %d\n",(*neigh)[i],cp->x[(*neigh)[i]],cp->y[(*neigh)[i]],cp->z[(*neigh)[i]],cp->FirstNeigh[(*neigh)[i]],cp->NumberNeigh[(*neigh)[i]],j,m);
					for(k=cp->FirstNeigh[(*neigh)[i]];k<cp->FirstNeigh[(*neigh)[i]]+cp->NumberNeigh[(*neigh)[i]];k++)
						printf("%d %f %f %f\n",cp->Neigh24[k],cp->x[cp->Neigh24[k]],cp->y[cp->Neigh24[k]],cp->z[cp->Neigh24[k]]);
					printf("\n");
					exit(0);
				}
				xt=A[0]*periodic_wrap(cp->x[m]-cp->x[i0],0)+A[1]*periodic_wrap(cp->y[m]-cp->y[i0],1)+A[2]*periodic_wrap(cp->z[m]-cp->z[i0],2);
				yt=A[3]*periodic_wrap(cp->x[m]-cp->x[i0],0)+A[4]*periodic_wrap(cp->y[m]-cp->y[i0],1)+A[5]*periodic_wrap(cp->z[m]-cp->z[i0],2);
				zt=A[6]*periodic_wrap(cp->x[m]-cp->x[i0],0)+A[7]*periodic_wrap(cp->y[m]-cp->y[i0],1)+A[8]*periodic_wrap(cp->z[m]-cp->z[i0],2);
				oct=(xt>0)+2*(yt>0)+4*(zt>0);

				flag=0;
				for(k=0;k<nn[oct];k++)
					if(m==temp_neigh[oct*nq+k])
						flag=1;
			
				if((flag==0)&&(m-i0))//this is a new neighbor
				{
					temp_neigh[oct*nq+nn[oct]]=m;
					temp_dis[oct*nq+nn[oct]]=get_squared_distance(cp->x[m]-cp->x[i0],cp->y[m]-cp->y[i0],cp->z[m]-cp->z[i0]);
					nn[oct]=nn[oct]+1;		
					if(nn[oct]>nq)	
					{	
						printf("Too many neighbours in %dth octant for %dth particle!\n",oct,i0);	
						exit(EXIT_FAILURE);	
					}
				}
			}
		}
		ite=ite+1;
		if(ite>5)
		{
			printf("Not enough neighbours.\n");
			printf("%d %f %f %f %d\n",i0,cp->x[i0],cp->y[i0],cp->z[i0],l);
			printf("\n");
			exit (EXIT_FAILURE);
		}
	}

	k=0;
	for(i=0;i<8;i++)
	{
		penalty=1.0;
		aquickSort(temp_dis,temp_neigh,k,k+nn[i]-1);
		for(j=0;j<nn[i];j++)
		{
			(*neigh)[k]=temp_neigh[i*nq+j];
			(*distance)[k]=penalty*temp_dis[i*nq+j];
			penalty=penalty*10000.0;
			k=k+1;
		}
	}

	aquickSort(*distance,*neigh,0,l-1);
	for(i=0;i<lm;i++)
		(*distance)[i]=sqrt(get_squared_distance(cp->x[(*neigh)[i]]-cp->x[i0],cp->y[(*neigh)[i]]-cp->y[i0],cp->z[(*neigh)[i]]-cp->z[i0]));
}



double Solver::periodic_wrap(double x, int dimension)
{
	double length;
  if(BC==1)
  {
		if(dimension==0)//x
			length=x_max-x_min;
		if(dimension==1)//y
			length=y_max-y_min;
		if(dimension==2)//z
			length=z_max-z_min;

   	while(x >= (0.5*length))
     	x -= length;

   	while(x < -(0.5*length))
     	x += length;
  }
  return x;
}

double Solver::get_squared_distance(double dx,double dy,double dz)
{
return periodic_wrap(dx,0)*periodic_wrap(dx,0)+periodic_wrap(dy,1)*periodic_wrap(dy,1)+periodic_wrap(dz,2)*periodic_wrap(dz,2);
}

void aquickSort(std::vector<double> &x, std::vector<int> &y, int l, int r)
{
   int j;

   if( l < r )
   {
   	// divide and conquer
        j = apartition( x, y, l, r);
       aquickSort( x, y, l, j-1);
       aquickSort( x, y, j+1, r);
   }
}



int apartition(std::vector<double> &a,std::vector<int> &b, int l, int r) 
{
   int i, j;
   double pivot, t;
   pivot = a[l];
   i = l; j = r+1;
   while( 1)
   {
   	do ++i; while( a[i] <= pivot && i <= r );
   	do --j; while( a[j] > pivot );
   	if( i >= j ) break;
   	t = a[i]; a[i] = a[j]; a[j] = t;
   	t = b[i]; b[i] = b[j]; b[j] = t;
   }
   t = a[l]; a[l] = a[j]; a[j] = t;
   t = b[l]; b[l] = b[j]; b[j] = t;
   return j;
}
