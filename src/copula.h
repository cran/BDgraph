#ifndef copula_H
#define copula_H

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include "matrix.h"

using namespace std;

extern "C" {
	void getMean( double Z[], double K[], double *muij, double *sigma, int *i, int *j, int *n, int *p );

	void getBounds( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copula( double Z[], double K[], int R[], int *n, int *p );
	 
	void getBoundsNA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copulaNA( double Z[], double K[], int R[], int *n, int *p );

	void getDs( double K[], double Z[], int R[], double D[], double Ds[], int *gcgm, int *n, int *p );

	void getTs( double Ds[], double Ts[], int *p );
}

#endif
