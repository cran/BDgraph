#include "rgwish.h"

// sampling from Wishart distribution
// Ts = chol( solve( Ds ) )
void rwish( double Ts[], double K[], int *b, int *p )
{
	int i, j, dim = *p, pxp = dim * dim, bK = *b;
	vector<double> psi( pxp, 0.0 ); 

	// ---- Sample values in Psi matrix ---
    GetRNGstate();
	for( i = 0; i < dim; i++ )
		psi[i * dim + i] = sqrt( rchisq( bK + dim - i - 1 ) );

	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
			psi[j * dim + i] = rnorm( 0, 1 );
	PutRNGstate();
	// ------------------------------------

    // C = psi %*% Ts   I used   psi = psi %*% Ts
	double alpha = 1.0, beta  = 0.0;
	char transT  = 'T', transN = 'N', side = 'R', upper = 'U';																	
	// dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
	F77_NAME(dtrmm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, Ts, &dim, &psi[0], &dim );

	// K = t(C) %*% C 
	// LAPACK function to compute  C := alpha * A * B + beta * C																				
	F77_NAME(dgemm)( &transT, &transN, &dim, &dim, &dim, &alpha, &psi[0], &dim, &psi[0], &dim, &beta, K, &dim );
}

// G is adjacency matrix which has zero in its diagonal
// threshold = 1e-8
void rgwish( int G[], double Ts[], double K[], int *b, int *p, double *threshold )
{
	char transN = 'N', uplo = 'U'; 
	double alpha = 1.0, beta  = 0.0;

	int info, i, j, l, size_i, one = 1, dim = *p, pxp = dim * dim;	
	double temp, threshold_C = *threshold;
	
	rwish( Ts, K, b, &dim );
	
	vector<double> sigma( pxp ); 
	inverse( K, &sigma[0], &dim );
	
	vector<double> W( sigma ); 
	vector<double> W_last( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 

	vector<double> sigma_N_i( dim );   // For dynamic memory used
	vector<int> N_i( dim );            // For dynamic memory used
	vector<double> W_N_i( pxp );       // For dynamic memory used

	double max_diff = 1.0;	
	while ( max_diff > threshold_C )
	{
		memcpy( &W_last[0], &W[0], sizeof( double ) * pxp );
		
		for( i = 0; i < dim; i++ )
		{
			// Count  size of note
			size_i = 0;
			for( j = 0; j < dim; j++ ) size_i += G[j * dim + i];

			if( size_i > 0 )
			{
				// Record size of node and initialize zero in beta_star for next steps
				sigma_N_i.resize( size_i );  // vector<double> sigma_N_i( size_i );
				N_i.resize( size_i );        // vector<int> N_i( size_i );
				
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					if( G[j * dim + i] )
					{
						sigma_N_i[l] = sigma[i * dim + j]; // sigma_N_i[j] = sigma[i * dim + N_i[j]];
						N_i[l++]     = j;
					}
					else
						beta_star[j] = 0.0; // for( j = 0; j < *p; j++ ) beta_star[j] = 0.0;
				}
				// -------------------------------------------------------------
				
				W_N_i.resize( size_i * size_i ); //vector<double> W_N_i( size_i * size_i );
				subMatrix( &W[0], &W_N_i[0], &N_i[0], &size_i, &dim );
					
				// A * X = B   for   sigma_N_i := (W_N_i)^{-1} * sigma_N_i
				F77_NAME(dposv)( &uplo, &size_i, &one, &W_N_i[0], &size_i, &sigma_N_i[0], &size_i, &info );

				for( j = 0; j < size_i; j++ ) beta_star[N_i[j]] = sigma_N_i[j];
				
				// multiplyMatrix( W, &beta_star[0], &sigma_i[0], &dim, &one, &dim );	// sigma_i = W * beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &W[0], &dim, &beta_star[0], &dim, &beta, &sigma_i[0], &dim );
				
				for( j = 0; j < i; j++ )
				{
					W[j * dim + i] = sigma_i[j];
					W[i * dim + j] = sigma_i[j];
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					W[j * dim + i] = sigma_i[j];
					W[i * dim + j] = sigma_i[j];
				}
			} 
			else 
			{
				for( j = 0; j < i; j++ )
				{
					W[j * dim + i] = 0.0;
					W[i * dim + j] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					W[j * dim + i] = 0.0;
					W[i * dim + j] = 0.0;
				}
			} 
		}

		max_diff = fabs( static_cast<double>( W[0] - W_last[0] ) );
		for( i = 1; i < pxp; i++ )
		{
			temp = fabs( static_cast<double>( W[i] - W_last[i] ) );
			if( temp > max_diff ) max_diff = temp; 
		}		
	}

	inverse( &W[0], K, &dim );
}
     
// rgwish ONLY for inside of MCMC algorithm
void rgwish_sigma( int G[], double Ts[], double K[], double sigma[], int *bstar, int *p, double *threshold,
					double sigma_start[], double invC[], double beta_star[], double sigma_i[], 
					vector<double> &sigma_start_N_i, vector<double> &sigma_N_i, vector<int> &N_i )
{
	char uplo = 'U';
	int i, j, ij, l, size_i, info, one = 1, dim = *p, pxp = dim * dim, bK = *bstar;	
	double temp, threshold_C = *threshold;
// STEP 1: sampling from wishart distributions
	//~ vector<double> sigma_start( pxp ); 
	// ---- Sample values in Psi matrix ---
    GetRNGstate();
	for( i = 0; i < dim; i++ )
		sigma_start[i * dim + i] = sqrt( rchisq( bK + dim - i - 1 ) );

	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			sigma_start[j * dim + i] = rnorm( 0, 1 );
			sigma_start[i * dim + j] = 0.0;
		}
	PutRNGstate();
	// ------------------------------------

    // C = psi %*% Ts   I used   psi = psi %*% Ts
	double alpha = 1.0; 
	char transT  = 'T', transN = 'N', side = 'R', upper = 'U';																	
	// dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
	F77_NAME(dtrmm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, Ts, &dim, &sigma_start[0], &dim );

	//~ vector<double> invC( pxp ); 
	side = 'L';
	// creating an identity matrix
	for( i = 0; i < dim; i++ )
		for( j = 0; j < dim; j++ )
			invC[j * dim + i] = (i == j);	
	// op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
	F77_NAME(dtrsm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, &sigma_start[0], &dim, &invC[0], &dim );

	// sigma_start <- invC %*% t(invC)   
	double beta  = 0.0;
	// LAPACK function to compute  C := alpha * A * B + beta * C																				
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &dim, &alpha, &invC[0], &dim, &invC[0], &dim, &beta, &sigma_start[0], &dim );

	memcpy( sigma, &sigma_start[0], sizeof( double ) * pxp ); 

	//~ vector<double> beta_star( dim ); 
	//~ vector<double> sigma_i( dim ); 
//~ 
	//~ vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	//~ vector<int> N_i( dim );            // For dynamic memory used
	//~ vector<double> sigma_N_i( pxp );       // For dynamic memory used
	
	double max_diff = 1.0;	
	while ( max_diff > threshold_C )
	{
		max_diff = 0.0;
		
		for( i = 0; i < dim; i++ )
		{
			// Count  size of note
			size_i = 0;
			for( j = 0; j < dim; j++ )  size_i += G[j * dim + i];

			if( size_i > 0 )
			{
				// Record size of node and initialize zero in beta_star for next steps
				sigma_start_N_i.resize( size_i );  // vector<double> sigma_start_N_i( size_i );
				N_i.resize( size_i );        // vector<int> N_i( size_i );
				
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					ij = j * dim + i;
					if( G[ij] )
					{
						sigma_start_N_i[l] = sigma_start[ij]; // sigma_start_N_i[j] = sigma_start[i * dim + N_i[j]];
						N_i[l++]     = j;
					}
					else
						beta_star[j] = 0.0; // for( j = 0; j < *p; j++ ) beta_star[j] = 0.0;
				}
				// -------------------------------------------------------------
				
				sigma_N_i.resize( size_i * size_i ); //vector<double> sigma_N_i( size_i * size_i );
				subMatrix( sigma, &sigma_N_i[0], &N_i[0], &size_i, &dim );
					
				// A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
				F77_NAME(dposv)( &uplo, &size_i, &one, &sigma_N_i[0], &size_i, &sigma_start_N_i[0], &size_i, &info );

				for( j = 0; j < size_i; j++ ) beta_star[N_i[j]] = sigma_start_N_i[j];
	
				// multiplyMatrix( sigma, &beta_star[0], &sigma_i[0], &dim, &one, &dim );	// sigma_i = sigma * beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, sigma, &dim, &beta_star[0], &dim, &beta, &sigma_i[0], &dim );

				for( j = 0; j < i; j++ )
				{
					ij   = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] - sigma_i[j] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij]          = sigma_i[j];
					sigma[i * dim + j] = sigma_i[j];
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					ij   = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] - sigma_i[j] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij]          = sigma_i[j];
					sigma[i * dim + j] = sigma_i[j];
				}
			} 
			else 
			{
				for( j = 0; j < i; j++ )
				{
					ij   = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij]          = 0.0;
					sigma[i * dim + j] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					ij = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij]          = 0.0;
					sigma[i * dim + j] = 0.0;				
				}
			} 
		}
	}
	
	memcpy( &sigma_start[0], sigma, sizeof( double ) * pxp );	 	
	
	inverse( &sigma_start[0], K, &dim );
}
        
