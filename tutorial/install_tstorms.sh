#!/bin/bash
# This script will clone, build, and install TSTORMS.

# It is assumed that the following dependencies are installed:
# - Fortran compiler (ifort is assumed, others require modification to the TSTORMS
#                     build scripts)
# - NetCDF with Fortran bindings

# Skip installation if TSTORMS is already built
if [ -f TSTORMS/tstorms_driver/tstorms_driver.exe ]; then
    echo "TSTORMS is already installed, skipping installation."
else
    # Clone and checkout specific commit TCTrack has been tested against
    git clone https://github.com/Cambridge-ICCS/TSTORMS.git
    cd TSTORMS/

    # Build using Make
    cd tstorms_driver/
    make FC=gfortran
    cd ../trajectory_analysis/
    make FC=gfortran

    # return to tutorial directory
    cd ../../
fi
