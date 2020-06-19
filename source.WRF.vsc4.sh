#!/bin/bash

# User specific aliases and functions
module purge
module load intel-mpi/2019.3 intel/19.1.0 netcdf-fortran/4.4.5-intel-19.0.5.281-qye4cqn zlib/1.2.11-intel-19.1.0.166-hs6m2qh  hdf5/1.10.5-intel-19.0.5.281-qyzojtm netcdf/4.7.0-intel-19.0.5.281-75t52g6 py-numpy/1.16.4-intel-19.0.5.281-l42jh75
# netcdf/4.7.0-intel-19.0.5.281-cgrpqof netcdf-fortran/4.4.5-intel-19.0.5.281-qye4cqn zlib/1.2.11-intel-19.1.0.166-hs6m2qh  hdf5/1.10.5-intel-19.0.5.281-qyzojtm
#module load gcc/4.8.5 hdf5/1.10.5-gcc-9.1.0-vm6bmqc netcdf/4.7.0-gcc-9.1.0-vtvq5rz netcdf-fortran/4.4.5-gcc-9.1.0-rrbh4ue py-numpy/1.17.4-gcc-9.1.0-sppplc4 zlib/1.2.11-gcc-9.1.0-ocd3sak mpi/openmpi-x86_64

setenv NETCDF /opt/sw/spack-0.12.1/opt/spack/linux-centos7-x86_64/intel-19.0.5.281/netcdf-4.7.0-75t52g6tfbjvhhnfh64cd6nvp2n7kewl 
