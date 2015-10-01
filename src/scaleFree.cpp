#include <R.h>
#include <Rmath.h>
#include <vector>        // for using vector

extern "C" {
	// For generating scale-free graphs: matrix G (p x p) is an adjacency matrix
	void scaleFree( int *G, int *p )
	{
		int i, j, tmp, total, dim = *p, p0 = 2;
		double randomValue;
		std::vector<int> size_a( dim ); 

		for( i = 0; i < p0 - 1; i++ )
		{
			G[i * dim + i + 1]   = 1;
			G[(i + 1) * dim + i] = 1;
		}
			
		for( i = 0; i < p0; i++ ) size_a[i] = 2;
		
		for( i = p0; i < dim; i++ ) size_a[i] = 0;
		
		total = 2 * p0;
		
		GetRNGstate();
		for( i = p0; i < dim; i++ )
		{
		   randomValue = (double) total * runif( 0, 1 );
		   
			tmp = 0;
			j   = 0;
			
			while ( tmp < randomValue && j < i ) 
				tmp += size_a[j++];
			
			j--;
			G[i * dim + j] = 1;
			G[j * dim + i] = 1;
			total += 2;
			size_a[j]++;
			size_a[i]++;
		}
		PutRNGstate();
	}
} // exturn "C"


