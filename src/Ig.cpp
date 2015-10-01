#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <vector>        // for using vector
#include <sstream>

using namespace std;

extern "C" {

// Part of function for calculating I.g = function( G, b = 3, D = diag( ncol(G) ), mc = 100 ) out logIg
void log_exp_mc( int G[], int nu[], int *b, double H[], int *check_H, int *mc, int *p, double f_T[] )
{
	int iter, i, j, ij, h, r, mc_iter = *mc, dim = *p, pxp = dim * dim, b_c = *b;
	double sumPsi, sumPsiH, sumPsiHi, sumPsiHj;
	vector<double> psi( pxp );      //building vector
	
	GetRNGstate();
	if( *check_H == 1 )
	{
		for( iter = 0; iter < mc_iter; iter++ ) //for ( iter in 1 : mc )
		{
			for( i = 0; i < dim; i++ )
				psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );

			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					if ( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 );
					psi[i * dim + j] = 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if ( G[ij] == 0 )
					{
						psi[ij] = 0.0;

						if ( i > 0 )  //if ( i > 1 )
						{
							sumPsi = 0.0;
							//sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] )
							for( h = 0; h < ( i - 1 ); h++ )
								sumPsi += ( psi[ i * dim + h ] * psi[ j * dim + h ] ); 
							
							//psi[i, j] <- - sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] ) / psi[i, i]
							psi[ij] = - sumPsi / psi[ i * dim + i ];
						}
						
						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[iter] += ( psi[ij] * psi[ij] ); 
					}
				}
		} 
	}
	else
	{
		for( iter = 0; iter < mc_iter; iter++ ) //for ( iter in 1 : mc )
		{
			for( i = 0; i < dim; i++ )
				psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );

			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					if ( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 );
					psi[i * dim + j] = 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if ( G[ij] == 0 )
					{
						//psi[i, j] = - sum( psi[ i, i : ( j - 1 ) ] * H[ i : ( j - 1 ), j ] )
						sumPsiH = 0.0;
						for( h = i; h < j; h++ )
							sumPsiH += ( psi[ h * dim + i ] * H[ j * dim + h ] ); 
						psi[ij] = - sumPsiH;
						
						if ( i > 0 )  //if ( i > 1 )
							for ( r = 0; r < i; r++ ) //for ( r in 1 : ( i - 1 ) )
							{
								//sum( psi[ r, r : i ] * H[ r : i, i ] )
								sumPsiHi = 0.0;
								for( h = r; h < i + 1; h++  )
									sumPsiHi += ( psi[ h * dim + r ] * H[ i * dim + h ] );
									
								//sum( psi[ r, r : j ] * H[ r : j, j ] ) )
								sumPsiHj = 0.0;
								for( h = r; h < j + 1; h++  )
									sumPsiHj += ( psi[ h * dim + r ] * H[ j * dim + h ] );
								
								//psi[i, j] <- psi[i, j] - ( ( sum( psi[ r, r : i ] * H[ r : i, i ] ) ) * ( sum( psi[ r, r : j ] * H[ r : j, j ] ) ) ) / ( psi[i, i] )
								psi[ij] -= ( sumPsiHi * sumPsiHj ) / psi[i * dim + i];
							}

						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[iter] += ( psi[ij] * psi[ij] ); 
					}
				}
		}
	}
	PutRNGstate();	
} 
     
} // End extern "C"
