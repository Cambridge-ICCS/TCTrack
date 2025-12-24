# This script will clone, build, and install Tempest Extremes.

# It is assumed that the following dependencies are installed:
# - CMake
# - C++ compiler
# - NetCDF with C++ bindings

# Clone and checkout specific commit TCTrack has been tested against
git clone https://github.com/ClimateGlobalChange/tempestextremes.git
cd tempestextremes/
git checkout 5feb3a04d29fd62a1f13fa9c0b85daeefcbecd6f

# Build using CMake
mkdir build
cmake -B build/ -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=OFF .
cmake --build build/

# Add Tempest extremes binaries to the path
export PATH="$PATH:$PWD/build/bin"

# return to tutorial directory
cd ../
