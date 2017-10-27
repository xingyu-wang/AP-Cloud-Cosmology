EXEC   = apcloud

PETSC_DIR=/home/xingyu/Downloads/petsc-3.5.1
OPENMPI_DIR=/usr/lib/openmpi
LAPACK_DIR=/home/xingyu/Downloads/lapack-3.5.0

CC       =  /usr/bin/mpic++
MPICHLIB = #-L/usr/mpi/gcc/openmpi-1.2.2-1/lib64 -lmpi

OBJS   = test.o parameter.o initialization.o controller.o solver.o mover.o octree.o viewer.o particle.o

INCL	= -I ./ -I ${PETSC_DIR}/include -I ${PETSC_DIR}/arch-linux2-c-debug/include  -I${OPENMPI_DIR}/include -I${OPENMPI_DIR}/include/openmpi ${LAPACK_DIR}/liblapacke.a ${LAPACK_DIR}/liblapack.a ${LAPACK_DIR}/librefblas.a

LFLAGS = -L ./ -Wl,-rpath,${PETSC_DIR}/arch-linux2-c-debug/lib -L ${PETSC_DIR}/arch-linux2-c-debug/lib  -Wl,-rpath,${OPENMPI_DIR}/lib -L${OPENMPI_DIR}/lib

OPTIONS= -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -O -fopenmp

LIBS   =  -lpetsc -llapack -lblas -lpthread -lm  -lgfortran $(MPICHLIB) -lgomp

OPTIMIZE =   -O3 -Wno-reorder

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)



apcloud: test.cpp parameter.cpp particle.cpp initialization.cpp controller.cpp solver.cpp mover.cpp octree.cpp viewer.cpp parameter.h particle.h initialization.h controller.h solver.h mover.h viewer.h octree.h apcloud_macro.h Makefile
	$(CC) $(CFLAGS) $(LFLAGS) test.cpp parameter.cpp particle.cpp initialization.cpp controller.cpp solver.cpp mover.cpp octree.cpp viewer.cpp $(INCL) $(LIBS) -o ./${EXEC}

run: 
	mpirun -np 2 -x OMP_NUM_THREADS ./${EXEC} ./parameter.param -ksp_gmres_restart 200 -ksp_max_it 1000 -ksp_atol 1.e-60 > ./../cosmology_result/cosmology_0922_adaptive.txt
ball: 
	mpirun -np 2 -x OMP_NUM_THREADS ./${EXEC} ./ball.param -ksp_gmres_restart 200 -ksp_max_it 1000 > ./../cosmology_result/ball_0922.txt

.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)
