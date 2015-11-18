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

	int info, i, j, l, size_node, one = 1, dim = *p, pxp = dim * dim;	
	double temp, threshold_C = *threshold;
	
	rwish( Ts, K, b, &dim );
	
	vector<double> sigma_start( pxp ); 
	inverse( K, &sigma_start[0], &dim );
	
	vector<double> sigma( sigma_start ); 
	vector<double> sigma_last( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_start_i( dim ); 

	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<int> N_i( dim );            // For dynamic memory used
	vector<double> sigma_N_i( pxp );       // For dynamic memory used

	double max_diff = 1.0;	
	while ( max_diff > threshold_C )
	{
		memcpy( &sigma_last[0], &sigma[0], sizeof( double ) * pxp );
		
		for( i = 0; i < dim; i++ )
		{
			// Count  size of note
			size_node = 0;
			for( j = 0; j < dim; j++ ) size_node += G[j * dim + i];

			if( size_node > 0 )
			{
				// Record size of node and initialize zero in beta_star for next steps
				sigma_start_N_i.resize( size_node );  // vector<double> sigma_start_N_i( size_node );
				N_i.resize( size_node );        // vector<int> N_i( size_node );
				
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					if( G[j * dim + i] )
					{
						sigma_start_N_i[l] = sigma_start[i * dim + j]; // sigma_start_N_i[j] = sigma_start[i * dim + N_i[j]];
						N_i[l++]     = j;
					}
					else
						beta_star[j] = 0.0; // for( j = 0; j < *p; j++ ) beta_star[j] = 0.0;
				}
				// -------------------------------------------------------------
				
				sigma_N_i.resize( size_node * size_node ); //vector<double> sigma_N_i( size_node * size_node );
				sub_matrix( &sigma[0], &sigma_N_i[0], &N_i[0], &size_node, &dim );
					
				// A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
				F77_NAME(dposv)( &uplo, &size_node, &one, &sigma_N_i[0], &size_node, &sigma_start_N_i[0], &size_node, &info );

				for( j = 0; j < size_node; j++ ) beta_star[N_i[j]] = sigma_start_N_i[j];
				
				// multiply_matrix( sigma, &beta_star[0], &sigma_start_i[0], &dim, &one, &dim );	// sigma_start_i = sigma * beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &sigma[0], &dim, &beta_star[0], &dim, &beta, &sigma_start_i[0], &dim );
				
				for( j = 0; j < i; j++ )
				{
					sigma[j * dim + i] = sigma_start_i[j];
					sigma[i * dim + j] = sigma_start_i[j];
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					sigma[j * dim + i] = sigma_start_i[j];
					sigma[i * dim + j] = sigma_start_i[j];
				}
			} 
			else 
			{
				for( j = 0; j < i; j++ )
				{
					sigma[j * dim + i] = 0.0;
					sigma[i * dim + j] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					sigma[j * dim + i] = 0.0;
					sigma[i * dim + j] = 0.0;
				}
			} 
		}

		max_diff = fabs( static_cast<double>( sigma[0] - sigma_last[0] ) );
		for( i = 1; i < pxp; i++ )
		{
			temp = fabs( static_cast<double>( sigma[i] - sigma_last[i] ) );
			if( temp > max_diff ) max_diff = temp; 
		}		
	}

	inverse( &sigma[0], K, &dim );
}
     
// rgwish ONLY for inside of MCMC algorithm
void rgwish_sigma( int G[], int size_node[], double Ts[], double K[], double sigma[], int *b_star, int *p, double *threshold,
					double sigma_start[], double inv_C[], double beta_star[], double sigma_i[], 
					vector<double> &sigma_start_N_i, vector<double> &sigma_N_i, vector<int> &N_i )
{
	char uplo = 'U';
	int i, j, ij, ip, l, size_node_i, info, one = 1, dim = *p, pxp = dim * dim, bK = *b_star;	
	double temp, threshold_C = *threshold;
// STEP 1: sampling from wishart distributions
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

	side = 'L';
	// creating an identity matrix
	for( i = 0; i < dim; i++ )
		for( j = 0; j < dim; j++ )
			inv_C[j * dim + i] = ( i == j );	
	// op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
	F77_NAME(dtrsm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, &sigma_start[0], &dim, &inv_C[0], &dim );

	// sigma_start <- inv_C %*% t(inv_C)   
	double beta  = 0.0;
	// LAPACK function to compute  C := alpha * A * B + beta * C																				
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &dim, &alpha, &inv_C[0], &dim, &inv_C[0], &dim, &beta, &sigma_start[0], &dim );

	memcpy( sigma, &sigma_start[0], sizeof( double ) * pxp ); 
	
	double max_diff = 1.0;	
	while ( max_diff > threshold_C )
	{
		max_diff = 0.0;
		
		for( i = 0; i < dim; i++ )
		{
			ip = i * dim;

			size_node_i = size_node[i];
			if( size_node_i > 0 )
			{
				// Record node size and initialize zero in beta_star
				sigma_start_N_i.resize( size_node_i );  // vector<double> sigma_start_N_i( size_node_i );
				N_i.resize( size_node_i );              // vector<int> N_i( size_node_i );
				
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					ij = ip + j;
					if( G[ij] )
					{
						sigma_start_N_i[l] = sigma_start[ij]; // sigma_start_N_i[j] = sigma_start[i * dim + N_i[j]];
						N_i[l++] = j;
					}
					else
						beta_star[j] = 0.0; // for( j = 0; j < *p; j++ ) beta_star[j] = 0.0;
				}
				// -------------------------------------------------------------
				
				sigma_N_i.resize( size_node_i * size_node_i ); //vector<double> sigma_N_i( size_node_i * size_node_i );
				sub_matrix( sigma, &sigma_N_i[0], &N_i[0], &size_node_i, &dim );
					
				// A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
				F77_NAME(dposv)( &uplo, &size_node_i, &one, &sigma_N_i[0], &size_node_i, &sigma_start_N_i[0], &size_node_i, &info );

				for( j = 0; j < size_node_i; j++ ) beta_star[N_i[j]] = sigma_start_N_i[j];
	
				// multiply_matrix( sigma, &beta_star[0], &sigma_i[0], &dim, &one, &dim );	// sigma_i = sigma * beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, sigma, &dim, &beta_star[0], &dim, &beta, &sigma_i[0], &dim );
				
				//~ for( j = 0; j < i; j++ ) sigma[ip + j] = sigma_i[j];	
				memcpy( sigma + ip, sigma_i,  sizeof( double ) * i );	
				
				for( j = 0; j < i; j++ )
				{
					ij   = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] - sigma_i[j] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij] = sigma_i[j];
				}
				
				//~ for( j = i + 1; j < dim; j++ ) sigma[ip + j] = sigma_i[j];	
				memcpy( sigma + ip + i + 1, sigma_i + i + 1, sizeof( double ) * ( dim - i - 1 ) );	

				for( j = i + 1; j < dim; j++ )
				{
					ij   = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] - sigma_i[j] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij] = sigma_i[j];
				}
			} 
			else 
			{				
				for( j = 0; j < i; j++ )
				{
					ij   = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij]     = 0.0;
					sigma[ip + j] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					ij = j * dim + i;
					temp = fabs( static_cast<double>( sigma[ij] ) );
					if( temp > max_diff ) max_diff = temp; 					

					sigma[ij]     = 0.0;
					sigma[ip + j] = 0.0;				
				}
			} 
		}
	}
	
	memcpy( &sigma_start[0], sigma, sizeof( double ) * pxp );	 	
	
	inverse( &sigma_start[0], K, &dim );
}
        
