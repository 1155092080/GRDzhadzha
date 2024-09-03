#!/bin/bash

#SBATCH -p astro3_short
#SBATCH -t 11:59:00
#SBATCH -N 1
#SBATCH -n 128
#SBATCH --ntasks-per-node 128
#SBATCH -J v01a05N48_M1_reflect
#SBATCH -o /lustre/astro/chen/rongzi/project/GRDzhadzha/Examples/MimickerBHScalar/v01a05N48_M1_reflect/outextract.txt
#SBATCH -e /lustre/astro/chen/rongzi/project/GRDzhadzha/Examples/MimickerBHScalar/v01a05N48_M1_reflect/errextract.txt
export HDF5_USE_FILE_LOCKING=FALSE
mpirun -n 128 /lustre/astro/chen/rongzi/project/GRDzhadzha/Examples/MimickerBHScalar/Main_MimickerBHScalar3d.Linux.64.mpicxx.gfortran.DEBUG.OPT.MPI.OPENMPCC.ex /lustre/astro/chen/rongzi/project/GRDzhadzha/Examples/MimickerBHScalar/params1_r_reflect.txt
