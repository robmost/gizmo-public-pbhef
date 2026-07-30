#!/bin/bash            # this line only there to enable syntax highlighting in this file
#
# Base configuration used by the build matrix in .github/workflows/build.yml. Each matrix entry appends
# its own flags to a copy of this file. Keep it close to the sbcluster test problem, so the matrix
# exercises a configuration we actually run.

BOX_PERIODIC
HYDRO_MESHLESS_FINITE_MASS
PMGRID=64
PM_HIRES_REGION_CLIPPING=1000
ADAPTIVE_GRAVSOFT_FORGAS
EOS_GAMMA=(5.0/3.0)
USE_FFTW3
DOUBLEPRECISION_FFTW
