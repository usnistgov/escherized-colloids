#!/bin/bash

# This is a very basic installation of LAMMPS using OpenMPI and GNU compilers - see https://docs.lammps.org/Install.html for detailed instructions

module load openmpi/4.0.5;
module load cmake/3.20.0;
git clone -b stable https://github.com/lammps/lammps.git lammps;
cd lammps;
mkdir build; cd build;
cmake ../cmake -DBUILD_MPI=yes -DBUILD_OMP=yes -DLAMMPS_MACHINE=mpi -DCMAKE_CXX_COMPILER=g++ -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPKG_OPENMP=yes -DPKG_RIGID=yes -DPKG_MOLECULE=yes -DLAMMPS_EXCEPTIONS=yes;
make -j8;
