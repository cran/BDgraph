#ifndef UTIL_H
#define UTIL_H

#include <Rconfig.h>

#if defined(_OPENMP) && __GNUG__ && defined(__linux__)
#define __PARALLEL__ true
#else
#define __PARALLEL__ false
#endif

#if __PARALLEL__
#define __NO_PARALLEL__ false
#else
#define __NO_PARALLEL__ true
#endif

#include <R.h>
#include <Rmath.h>

#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <R_ext/Complex.h>
#include <R_ext/Print.h>
#include <R_ext/Arith.h>     // for the special values like NA, NaN  

#include <R_ext/Visibility.h>

#if __PARALLEL__
#include <omp.h>
#endif

#include <sstream>
#include <string>            // std::string, std::to_string
#include <vector>            // for using vector

#include <math.h>            // isinf, sqrt
#include <limits>            // for std::numeric_limits<double>::max()
#include <algorithm>         // for transform function
#include <functional>        // for transform function
#include <climits>

using namespace std;

#endif













