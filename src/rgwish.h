#ifndef rgwish_H
#define rgwish_H

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
	void rwish( double Ts[], double K[], int *b, int *p );

	void rgwish( int G[], double Ts[], double K[], int *b, int *p, double *threshold );

	void rgwish_sigma( int G[], int size_node[], double Ts[], double K[], double sigma[], int *b_star, int *p, double *threshold,
					double sigma_start[], double inv_C[], double beta_star[], double sigma_i[], 
					vector<double> &sigma_start_N_i, vector<double> &sigma_N_i, vector<int> &N_i );

	void log_exp_mc( int G[], int nu[], int *b, double H[], int *check_H, int *mc, int *p, double f_T[] );
}

#endif
