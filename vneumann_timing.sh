#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/cosmology

export Number_P=1
export OMP_NUM_THREADS=1
make -f vneumann_makefile ball

export Number_P=2
export OMP_NUM_THREADS=1
make -f vneumann_makefile ball

export Number_P=4
export OMP_NUM_THREADS=1
make -f vneumann_makefile ball

export Number_P=8
export OMP_NUM_THREADS=1
make -f vneumann_makefile ball

export Number_P=16
export OMP_NUM_THREADS=1
make -f vneumann_makefile ball

export Number_P=1
export OMP_NUM_THREADS=2
make -f vneumann_makefile ball

export Number_P=1
export OMP_NUM_THREADS=4
make -f vneumann_makefile ball

export Number_P=1
export OMP_NUM_THREADS=8
make -f vneumann_makefile ball

export Number_P=2
export OMP_NUM_THREADS=16
make -f vneumann_makefile ball

