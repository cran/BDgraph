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

#include "util.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

extern "C" {
	void omp_set_num_cores( int *cores, int *verbose_core ) 
	{
	    #ifdef _OPENMP
	        omp_set_num_threads( *cores );
	    #else
	        if( *verbose_core == 1 ) Rprintf( "  This OS does not support multi-threading for the BDgraph package  \n" ); 
	    #endif
	}
}
