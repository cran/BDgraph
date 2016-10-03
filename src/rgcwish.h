#ifndef rgcwish_H
#define rgcwish_H

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <R_ext/Complex.h>
#include "MyLapack.h"
#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include "matrix.h"

using namespace std;

extern "C" {
	void rcwish( double Ls[], Rcomplex *K, int *b, int *p );

	void rgcwish( int G[], double Ls[], Rcomplex *K, int *b, int *p, double *threshold );
	
	void rgcwish_sigma( int G[], int size_node[], double Ls[], Rcomplex *K, double r_sigma[], double i_sigma[], Rcomplex *csigma, Rcomplex *Ind, int *b_star, int *p, double *threshold,
					double r_sigma_start[], double i_sigma_start[], double X[], double Y[], double r_beta_star[], double i_beta_star[], double joint[], double r_sigma_i[], 
					double i_sigma_i[], vector<double> &r_sigma_start_N_i, vector<double> &r_sigma_start_N_i_2, vector<double> &i_sigma_start_N_i, vector<double> &r_sigma_N_i, vector<double> &r_sigma_N_i_2,
					vector<double> &i_sigma_N_i, vector<int> &N_i, vector<double> &IR, vector<double> &Inv_R );
}

#endif
