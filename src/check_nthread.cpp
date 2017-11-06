#include "util.h"

extern "C" {
  
	void check_nthread( int *R_nthread )
	{
		if( !__PARALLEL__ )
			*R_nthread = 1;
	}
  
}
