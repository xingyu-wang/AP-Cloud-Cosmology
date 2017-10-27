#!/bin/bash
#PBS -l nodes=1:ppn=24

cd /home/xingyu/cosmology

export Number_P=1
export OMP_NUM_THREADS=4

make -f vneumann_makefile ball
