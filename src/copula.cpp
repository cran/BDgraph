#include "copula.h"
   
// for copula function
void getMean( double Z[], double K[], double *muij, double *sigma, int *i, int *j, int *n, int *p )
{
	int k, dim = *p, number = *n, row = *i, col = *j;
	double mu_ij = 0.0;
	
	for( k = 0; k < col; k++ ) 
		mu_ij += Z[k * number + row] * K[col * dim + k];

	for( k = col + 1; k < dim; k++ ) 
		mu_ij += Z[k * number + row] * K[col * dim + k];

	*muij = - mu_ij * *sigma;
}

// for copula function
void getBounds( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n )
{
	int kj, ij, row = *i, col = *j;
	double lowb = -1e308, upperb = +1e308;

	for( int k = 0, number = *n; k < number; k++ )
	{
		kj = col * number + k;
		ij = col * number + row;
		
     // if( R[k, j] < R[i, j] ) lb = max( Z[ k, j], lb )
	 // if( R[k, j] > R[i, j] ) ub = min( Z[ k, j], ub )										
		if( R[kj] < R[ij] ) 
			lowb = max( Z[kj], lowb );	
		else if( R[kj] > R[ij] ) 
			upperb = min( Z[kj], upperb );
	}

	*lb = lowb;
	*ub = upperb;	
}
 
// copula part
void copula( double Z[], double K[], int R[], int *n, int *p )
{
	int number = *n, dim = *p;
	double sigma, sdj, muij;
	
	double lb, ub;
	double runifValue, pnormlb, pnormub;
	
	for( int j = 0; j < dim; j++ )
	{   
		sigma = 1 / K[j * dim + j];
		sdj   = sqrt( sigma );
		
		// interval and sampling
		for( int i = 0; i < number; i++ )
		{
			getMean( Z, K, &muij, &sigma, &i, &j, &number, &dim );
			
			getBounds( Z, R, &lb, &ub, &i, &j, &number );
			
			// runifValue = runif( 1, pnorm( lowerBound, muij, sdj ), pnorm( upperBound, muij, sdj ) )
			// Z[i,j]     = qnorm( runifValue, muij, sdj )									
			GetRNGstate();
			pnormlb           = pnorm( lb, muij, sdj, TRUE, FALSE );
			pnormub           = pnorm( ub, muij, sdj, TRUE, FALSE );
			runifValue        = runif( pnormlb, pnormub );
			Z[j * number + i] = qnorm( runifValue, muij, sdj, TRUE, FALSE );
			PutRNGstate();				
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// for copula function with missing data 
void getBoundsNA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n )
{
	int kj, ij, row = *i, col = *j;
	double lowb = -1e308, upperb = +1e308;

	for( int k = 0, number = *n; k < number; k++ )
	{
		kj = col * number + k;
		ij = col * number + row;
		
		if( R[kj] != 0 )
		{
			// if( R[k, j] < R[i, j] ) lb = max( Z[ k, j], lb )
			// if( R[k, j] > R[i, j] ) ub = min( Z[ k, j], ub )										
			if( R[kj] < R[ij] ) 
				lowb = max( Z[kj], lowb );	
			else if( R[kj] > R[ij] ) 
				upperb = min( Z[kj], upperb );
		}
	}
	
	*lb = lowb;
	*ub = upperb;		
}
 
// copula part
void copulaNA( double Z[], double K[], int R[], int *n, int *p )
{
	int ij, number = *n, dim = *p;
	double sigma, sdj, muij, lb, ub, runifValue, pnormlb, pnormub;
	
	for( int j = 0; j < dim; j++ )
	{   
		sigma = 1 / K[j * dim + j];
		sdj   = sqrt( sigma );
		
		// interval and sampling
		for( int i = 0; i < number; i++ )
		{
			getMean( Z, K, &muij, &sigma, &i, &j, &number, &dim );
			
			ij = j * number + i;
			GetRNGstate();
			if( R[ij] != 0 )
			{
				getBoundsNA( Z, R, &lb, &ub, &i, &j, &number );
				
				pnormlb    = pnorm( lb, muij, sdj, TRUE, FALSE );
				pnormub    = pnorm( ub, muij, sdj, TRUE, FALSE );
				runifValue = runif( pnormlb, pnormub );
				Z[ij]      = qnorm( runifValue, muij, sdj, TRUE, FALSE );
			} 
			else 
				Z[ij] = rnorm( muij, sdj );
			PutRNGstate();				
		}
	}
}
    
// for bdmcmcCopulaDmh function
void getDs( double K[], double Z[], int R[], double D[], double Ds[], int *gcgm, int *n, int *p )
{
	int gcgmCheck = *gcgm, dim = *p, pxp = dim * dim;

	if( gcgmCheck == 0 )
		copula( Z, K, R, n, &dim );
	else
		copulaNA( Z, K, R, n, &dim );
	
	vector<double> S( pxp ); 
	// S <- t(Z) %*% Z
	// Here, I'm using Ds instead of S, for saving memory
	double alpha = 1.0, beta  = 0.0;
	char transA = 'T', transB = 'N';
	// LAPACK function to compute  C := alpha * A * B + beta * C
	//        DGEMM ( TRANSA,  TRANSB, M, N, K,  ALPHA, A,LDA,B, LDB,BETA, C, LDC )																				
	F77_NAME(dgemm)( &transA, &transB, &dim, &dim, n, &alpha, Z, n, Z, n, &beta, &S[0], &dim );		
	// Ds = D + S
	for( int i = 0; i < pxp ; i++ ) Ds[i] = D[i] + S[i];		
}

// for bdmcmcCopulaDmh function
void getTs( double Ds[], double Ts[], int *p )
{
	int dim = *p, pxp = dim * dim;
	vector<double> invDs( pxp ); 
	vector<double> copyDs( pxp ); 

	//~ for( int i = 0; i < pxp; i++ ) copyDs[i] = Ds[i]; 	
	memcpy( &copyDs[0], Ds, sizeof( double ) * pxp );
	
	inverse( &copyDs[0], &invDs[0], &dim );	

	cholesky( &invDs[0], Ts, &dim );	
}
    
