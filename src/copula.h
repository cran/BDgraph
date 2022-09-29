// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2012 - 2020  Reza Mohammadi                                                   |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  
#ifndef copula_H
#define copula_H

#include "matrix.h"

extern "C" {
	void get_mean( double Z[], double K[], double *mu_ij, double *sigma, int *i, int *j, int *n, int *p );

	void get_bounds( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copula( double Z[], double K[], int R[], int not_continuous[], int *n, int *p );
	
	void copula_dw( double Z[], double K[], int Y[], double lower_bounds[], double upper_bounds[], int *n, int *p );
	    
    void copula_dw_NA( double Z[], double K[], int Y[], double lower_bounds[], double upper_bounds[], int *n, int *p );

    void get_Ds_dw( double K[], double Z[], int Y[], double lower_bounds[], double upper_bounds[], double D[], double Ds[], double S[], int *gcgm, int *n, int *p );
	    
	void get_bounds_NA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n );

	void copula_NA( double Z[], double K[], int R[], int not_continuous[], int *n, int *p );

	void get_Ds( double K[], double Z[], int R[], int not_continuous[], double D[], double Ds[], double S[], int *gcgm, int *n, int *p );

	void get_Ts( double Ds[], double Ts[], double inv_Ds[], double copy_Ds[], int *p );
	
    void update_tu( double data[], double K[], double tu[], double mu[], double *nu, int *n, int *p );

    void get_Ds_tgm( double data[], double D[], double mu[], double tu[], double Ds[], double S[], int *n, int *p );
    
    void update_mu( double data[], double mu[], double tu[], int *n, int *p );

}

#endif
