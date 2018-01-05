#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {
	void omp_set_num_cores( int *cores ) 
	{
		#ifdef _OPENMP
			omp_set_num_threads( *cores );
		#endif
	}
}
