
/*! \file includeDef.h
list of included libraries for the TPSC program
created by Dominic Bergeron
*/


#ifdef USE_MPI
#include "fftw3-mpi.h"
#include "mpi.h"
#endif

#include "fftw3.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <complex>
#include <string>
#include <cstring>
#include <ctime>
#include <limits>
#include <sys/stat.h>
#include <omp.h>
//#include <gsl/gsl_sf_fermi_dirac.h>
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_math.h>

#ifndef PI
#define PI acos((double)-1.0)
#endif

//#ifndef USE_MPI
//#define RELOAD_DATA 1
//#endif

#define EPSILON numeric_limits<double>::epsilon()

using namespace std;

typedef complex<double> dcomplex;

#ifndef I
#define I dcomplex(0.0,1.0)
#endif
