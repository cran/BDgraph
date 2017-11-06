#include "util.h"

extern "C" {
  
	void check_os( int *R_os )
	{
		if( __PARALLEL__)
		{
			*R_os = 1;
		}else{
			*R_os = 0;
		}
	}
}
