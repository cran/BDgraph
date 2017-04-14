// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2016 Mohammadi A. and Wit C. E.
//
//     This file is part of BDgraph package.
//
//     BDgraph is free software: you can redistribute it and/or modify it under 
//     the terms of the GNU General Public License as published by the Free 
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.
//
//     Maintainer:
//     Abdolreza Mohammadi: a.mohammadi@rug.nl or a.mohammadi@uvt.nl
// ----------------------------------------------------------------------------|
  
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
	void get_mean( double Z[], double K[], double *mu_ij, double *sigma, int *i, int *j, int *n, int *p );

	void get_bounds( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copula( double Z[], double K[], int R[], int *n, int *p );
	 
	void get_bounds_NA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copula_NA( double Z[], double K[], int R[], int *n, int *p );

	void get_Ds( double K[], double Z[], int R[], double D[], double Ds[], double S[], int *gcgm, int *n, int *p );

	void get_Ts( double Ds[], double Ts[], double inv_Ds[], double copy_Ds[], int *p );
}

#endif
