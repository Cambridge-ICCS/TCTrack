# This script will clone, build, and install TSTORMS.

# It is assumed that the following dependencies are installed:
# - Fortran compiler (ifort is assumed, others require modification to the TSTORMS
#                     build scripts)
# - NetCDF with Fortran bindings

# Clone and checkout specific commit TCTrack has been tested against
git clone https://github.com/Cambridge-ICCS/TSTORMS.git
cd TSTORMS/

# Build using Make
cd tstorms_driver/
make
cd ../trajectory_analysis/
make

# return to tutorial directory
cd ../../
