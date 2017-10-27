#!/bin/bash
#PBS -l nodes=1:towel:ppn=12

cd /nfs/user02/xywang/cosmology/cosmology

export OMP_NUM_THREADS=2

make -f vogon_makefile run
