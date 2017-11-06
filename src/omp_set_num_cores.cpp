#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

extern "C" {
	void omp_set_num_cores( int *cores ) 
	{
		#ifdef SUPPORT_OPENMP
			omp_set_num_threads( *cores );
		#endif
	}
}
