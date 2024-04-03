#! /bin/bash

# Change the paths to your 

module load StdEnv/2020
module load gcc/9.3.0
module load openblas/0.3.17
module load netcdf-fortran/4.5.2 # HDF5/1.10.6
module load boost

# Point ARMA_INCLUDES to the include directory of your Armadillo installation 
export ARMA_INCLUDES=/globalhome/kck540/HPC/OpenWQ-Projects/armadillo-10.3.0/include
# Point ARMA_LIB to the libarmadillo.so.10 file of you Armadillo installation
export ARMA_LIB=/globalhome/kck540/HPC/OpenWQ-Projects/armadillo-10.3.0/build/libarmadillo.so.10


cmake -S. -B_build -DCMAKE_BUILD_TYPE=debug
cmake --build _build -j 4