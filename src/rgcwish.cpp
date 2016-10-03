#include "rgcwish.h"

// sampling from COMPLEX Wishart distribution
// Ls = t( chol( solve( Ds ) ) )
void rcwish( double Ls[], Rcomplex *K, int *b, int *p )
{
	int i, j, dim = *p, bK = *b, n = bK + dim, p2 = 2 * dim, pxn = dim * n, p2xn = 2 * pxn, pxp = dim * dim, info;
	int i2p, ip;
	double *joint = new double[p2xn];
	double *X = new double[pxn];
	double *Y = new double[pxn];
	vector<double> r_K( pxp );
	vector<double> i_K( pxp );
	Rcomplex *csigma = new Rcomplex[pxp];
	Rcomplex *Ind = new Rcomplex[pxp];

	// ---- Sample values in Joint matrix ---
	GetRNGstate();
	for( j = 0; j < n; j++ )
		for( i = 0; i < p2; i++ )
		{
			joint[j * p2 + i] = rnorm( 0, 1 );
		}
	PutRNGstate();
	// ------------------------------------

	double alpha = 1.0, malpha = -1.0; 
	char transT  = 'T', transN = 'N', side = 'L', upper = 'U', lower = 'L';																	
	// dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
	// C = Ls %*% joint   I used   joint = Ls %*% joint
	F77_NAME(dtrmm)( &side, &lower, &transN, &transN, &p2, &n, &alpha, Ls, &p2, &joint[0], &p2 );
	for ( i = 0; i < n; i++ )
	{
		i2p = i * p2;
		ip = i * dim;
		memcpy(X + ip, &joint[i2p], sizeof(double) * dim);
		memcpy(Y + ip, &joint[i2p + dim], sizeof(double) * dim);
	}
	double beta  = 0.0;
	// The real part of K_start = X %*% t(X) + Y %*% t(Y)
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &X[0], &dim, &X[0], &dim, &beta, &r_K[0], &dim );
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &Y[0], &dim, &Y[0], &dim, &alpha, &r_K[0], &dim );

	// The imaginary part of K_start = Y %*% t(X) - X %*% t(Y)
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &Y[0], &dim, &X[0], &dim, &beta, &i_K[0], &dim );
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &malpha, &X[0], &dim, &Y[0], &dim, &alpha, &i_K[0], &dim );
	
	for (j = 0; j < pxp; j++)
	{
		K[j].r = r_K[j];
		K[j].i = i_K[j];
	}
	delete[] csigma, Ind, joint, X, Y;
}

// sampling from COMPLEX G-Wishart distribution
void rgcwish( int G[], double Ls[], Rcomplex *K, int *b, int *p, double *threshold )
{
	int i, j, ij, ji, i2p, ip, l, size_node_i, info, one = 1, dim = *p, p2 = 2*dim, pxp = dim * dim, bK = *b, n = bK + dim, pxn = dim * n, p2xn = 2 * pxn;	
	double temp, threshold_C = *threshold, done = 1.0, dmone = -1.0;
	double alpha = 1.0, malpha = -1.0, beta  = 0.0; 
	char transT  = 'T', transN = 'N', side = 'L', upper = 'U', lower = 'L';	
	// STEP 1: sampling from complex wishart distributions
	double *joint = new double[p2xn];
	double *X = new double[pxn];
	double *Y = new double[pxn];

	vector<double> r_sigma_start(pxp);
	vector<double> i_sigma_start(pxp);
	Rcomplex *csigma = new Rcomplex[pxp];
	Rcomplex *Ind = new Rcomplex[pxp];

	// ---- Sample values in Joint matrix ---
	GetRNGstate();
	for( j = 0; j < n; j++ )
		for( i = 0; i < p2; i++ )
		{
			joint[j * p2 + i] = rnorm( 0, 1 );
		}
	PutRNGstate();
	// ------------------------------------
																
	// dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
	// C = Ls %*% joint   I used   joint = Ls %*% joint
	F77_NAME(dtrmm)( &side, &lower, &transN, &transN, &p2, &n, &alpha, Ls, &p2, &joint[0], &p2 );
	for ( i = 0; i < n; i++ )
	{
		i2p = i * p2;
		ip = i * dim;
		memcpy(X + ip, &joint[i2p], sizeof(double) * dim);
		memcpy(Y + ip, &joint[i2p + dim], sizeof(double) * dim);
	}
	// The real part of K_start = X %*% t(X) + Y %*% t(Y)
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &X[0], &dim, &X[0], &dim, &beta, &r_sigma_start[0], &dim );
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &Y[0], &dim, &Y[0], &dim, &alpha, &r_sigma_start[0], &dim );

	// The imaginary part of K_start = Y %*% t(X) - X %*% t(Y)
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &Y[0], &dim, &X[0], &dim, &beta, &i_sigma_start[0], &dim );
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &malpha, &X[0], &dim, &Y[0], &dim, &alpha, &i_sigma_start[0], &dim );

	for ( j = 0; j < dim; j++ )
	  for ( int k = 0; k < dim; k++ ){
	    csigma[j*dim+k].r = r_sigma_start[j*dim+k];
	    csigma[j*dim+k].i = i_sigma_start[j*dim+k];
  	  }
	for ( j = 0; j < dim; j++ )
          for ( int k = 0; k < dim; k++ ){
            Ind[j*dim+k].i = 0;
            Ind[j*dim+k].r = (j == k);
          }
	F77_NAME(zpotrf)( &upper, &dim, csigma, &dim, &info ); // Remark: csigma will be changed
	zpotrs( &upper, &dim, &dim, csigma, &dim, Ind, &dim, &info ); //sigma = inv(K)
	
	for ( j = 0; j < pxp; j++ ){
	  r_sigma_start[j] = Ind[j].r;
	  i_sigma_start[j] = Ind[j].i;
	}
	
	vector<double> r_sigma( r_sigma_start );
	vector<double> i_sigma( i_sigma_start );	 
	vector<double> r_beta_star( dim );
	vector<double> i_beta_star( dim );	 
	vector<double> r_sigma_i( dim );	 
	vector<double> i_sigma_i( dim );

	// For dynamic memory used
	vector<double> r_sigma_start_N_i( dim );	
	vector<double> r_sigma_start_N_i_2( dim );	
	vector<double> i_sigma_start_N_i( dim );	
	vector<int> N_i( dim );	
	vector<double> r_sigma_N_i( dim );	
	vector<double> r_sigma_N_i_2( dim );
	vector<double> i_sigma_N_i( dim );	
	vector<double> Inv_R( pxp );	
	vector<double> IR( pxp );	

	
	double max_diff = 1.0, r_diff, i_diff;	
	while ( max_diff > threshold_C )
	{
		max_diff = 0.0; // The maximum difference between two adjacent sigma

		for( i = 0; i < dim; i++ )
		{
			ip = i * dim;
			// Count size of node
			size_node_i = 0;
			for( j = 0; j < dim; j++ ) size_node_i += G[j * dim + i];

			if( size_node_i > 0 )
			{
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					ji = ip + j;
					if( G[ji] )
					{
						r_sigma_start_N_i[l] = r_sigma_start[ji];
						r_sigma_start_N_i_2[l] = r_sigma_start[ji];
						i_sigma_start_N_i[l] = i_sigma_start[ji];
						N_i[l++] = j;
					}
					else
					{
						r_beta_star[j] = 0.0; 
						i_beta_star[j] = 0.0; 
					}
				}
				// -------------------------------------------------------------
				
				sub_matrix( &r_sigma[0], &r_sigma_N_i[0], &N_i[0], &size_node_i, &dim );
				sub_matrix( &i_sigma[0], &i_sigma_N_i[0], &N_i[0], &size_node_i, &dim );
				
				// Inv_R = solve(r_sigma_N_i)
				for (int s = 0; s < size_node_i*size_node_i; s++)
				  r_sigma_N_i_2[s] = r_sigma_N_i[s];
				inverse( &r_sigma_N_i_2[0], &Inv_R[0], &size_node_i );			
				// IR = i_sigma_N_i %*% Inv_R
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &size_node_i, &size_node_i, &alpha, &i_sigma_N_i[0], &size_node_i, &Inv_R[0], &size_node_i, &beta, &IR[0], &size_node_i);
				// A = IR %*% i_sigma_N_i + r_sigma_N_i = r_sigma_N_i
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &size_node_i, &size_node_i, &alpha, &IR[0], &size_node_i, &i_sigma_N_i[0], &size_node_i, &alpha, &r_sigma_N_i[0], &size_node_i);
				// B = IR %*% i_sigma_start_N_i + r_sigma_start_N_i	= r_sigma_start_N_i			
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &one, &size_node_i, &alpha, &IR[0], &size_node_i, &i_sigma_start_N_i[0], &size_node_i, &alpha, &r_sigma_start_N_i[0], &size_node_i);				
				// C = IR %*% r_sigma_start_N_i_2 - i_sigma_start_N_i = i_sigma_start_N_i
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &one, &size_node_i, &alpha, &IR[0], &size_node_i, &r_sigma_start_N_i_2[0], &size_node_i, &dmone, &i_sigma_start_N_i[0], &size_node_i);			
				
				// r_sigma_start_N_i = solve(A) %*% B; i_sigma_start_N_i = solve(A) %*% C
				for (int s = 0; s < size_node_i*size_node_i; s++)
				  r_sigma_N_i_2[s] = r_sigma_N_i[s];
				F77_NAME(dposv)( &upper, &size_node_i, &one, &r_sigma_N_i[0], &size_node_i, &r_sigma_start_N_i[0], &size_node_i, &info );
				F77_NAME(dposv)( &upper, &size_node_i, &one, &r_sigma_N_i_2[0], &size_node_i, &i_sigma_start_N_i[0], &size_node_i, &info );	

				for( j = 0; j < size_node_i; j++ ) 
				{
					r_beta_star[N_i[j]] = r_sigma_start_N_i[j];
					i_beta_star[N_i[j]] = i_sigma_start_N_i[j];
				}
				
				// r_sigma_i = r_sigma %*% r_beta_star + i_sigma %*% i_beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &i_sigma[0], &dim, &i_beta_star[0], &dim, &beta, &r_sigma_i[0], &dim );
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &r_sigma[0], &dim, &r_beta_star[0], &dim, &done, &r_sigma_i[0], &dim );

				// i_sigma_i = i_sigma %*% r_beta_star - r_sigma %*% i_beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &r_sigma[0], &dim, &i_beta_star[0], &dim, &beta, &i_sigma_i[0], &dim );
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &i_sigma[0], &dim, &r_beta_star[0], &dim, &dmone, &i_sigma_i[0], &dim );
					
				// Update the first i elements
				for( j = 0; j < i; j++ )
				{
					ij = j * dim + i;
					ji = i * dim + j;
					r_diff = r_sigma[ij] - r_sigma_i[j];
					i_diff = i_sigma[ij] + i_sigma_i[j];

					temp = sqrt( r_diff * r_diff + i_diff * i_diff );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij] = r_sigma_i[j];
					i_sigma[ij] = - i_sigma_i[j];
					r_sigma[ji] = r_sigma_i[j];
					i_sigma[ji] = i_sigma_i[j];
				}
				

				for( j = i + 1; j < dim; j++ )
				{
					ij = j * dim + i;
					ji = i * dim + j;					
					r_diff = r_sigma[ij] - r_sigma_i[j];
					i_diff = i_sigma[ij] + i_sigma_i[j];

					temp = sqrt( r_diff * r_diff + i_diff * i_diff );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij] = r_sigma_i[j];
					i_sigma[ij] = - i_sigma_i[j];
					r_sigma[ji] = r_sigma_i[j];
					i_sigma[ji] = i_sigma_i[j];					
				}
			} 
			else 
			{				
				for( j = 0; j < i; j++ )
				{
					ij   = j * dim + i;
					temp = sqrt( r_sigma[ij] * r_sigma[ij] + i_sigma[ij] * i_sigma[ij] );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij]     = 0.0;
					i_sigma[ij]     = 0.0;
					r_sigma[ip + j] = 0.0;
					i_sigma[ip + j] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					ij = j * dim + i;
					temp = sqrt( r_sigma[ij] * r_sigma[ij] + i_sigma[ij] * i_sigma[ij] );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij]     = 0.0;
					i_sigma[ij]     = 0.0;
					r_sigma[ip + j] = 0.0;				
					i_sigma[ip + j] = 0.0;				
				}
			} 
		}
	}
	
	memcpy( &r_sigma_start[0], &r_sigma[0], sizeof( double ) * pxp );	 	
	memcpy( &i_sigma_start[0], &i_sigma[0], sizeof( double ) * pxp );	 	
	
	// K = solve(sigma)
	for (int j = 0; j < pxp; j++){
	  csigma[j].r = r_sigma_start[j];
	  csigma[j].i = i_sigma_start[j];
	}
	for (int j = 0; j < dim; j++)
		for (int k = 0; k < dim; k++){
		K[j*dim+k].i = 0;
		K[j*dim+k].r = (j == k);
		}

	F77_NAME(zpotrf)(&upper, &dim, csigma, &dim, &info);
	zpotrs(&upper, &dim, &dim, csigma, &dim, K, &dim, &info );
	delete[] csigma, Ind, joint, X, Y;
}

// rgwish ONLY for inside of MCMC algorithm
// Ls is the cholesky of the covariance matrix of (X;Y) -- sigma = Ls %*% Ls*
// Input (value matters) -- G, size_node, Ls, b_star, p, threshold
// Output -- K, r_sigma, i_sigma
void rgcwish_sigma( int G[], int size_node[], double Ls[], Rcomplex *K, double r_sigma[], double i_sigma[], Rcomplex *csigma, Rcomplex *Ind, int *b_star, int *p, double *threshold,
					double r_sigma_start[], double i_sigma_start[], double X[], double Y[], double r_beta_star[], double i_beta_star[], double joint[], double r_sigma_i[], 
					double i_sigma_i[], vector<double> &r_sigma_start_N_i, vector<double> &r_sigma_start_N_i_2, vector<double> &i_sigma_start_N_i, vector<double> &r_sigma_N_i, vector<double> &r_sigma_N_i_2,
					vector<double> &i_sigma_N_i, vector<int> &N_i, vector<double> &IR, vector<double> &Inv_R )
{
	int i, j, ij, ji, i2p, ip, l, size_node_i, info, one = 1, dim = *p, p2 = 2*dim, pxp = dim * dim, bK = *b_star, n = bK + dim, p2xn = 2*dim*n;	
	double temp, threshold_C = *threshold, done = 1.0, dmone = -1.0;
	// STEP 1: sampling from complex wishart distributions
	// ---- Sample values in Joint matrix ---
	GetRNGstate();
	for( j = 0; j < n; j++ )
		for( i = 0; i < p2; i++ )
		{
			joint[j * p2 + i] = rnorm( 0, 1 );
		}
	PutRNGstate();
	// ------------------------------------

	double alpha = 1.0, malpha = -1.0; 
	char transT  = 'T', transN = 'N', side = 'L', upper = 'U', lower = 'L';																	
	// dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
	// C = Ls %*% joint   I used   joint = Ls %*% joint
	F77_NAME(dtrmm)( &side, &lower, &transN, &transN, &p2, &n, &alpha, Ls, &p2, &joint[0], &p2 );
	for (i = 0; i < n; i++)
	{
		i2p = i * p2;
		ip = i * dim;
		memcpy(X + ip, &joint[i2p], sizeof(double) * dim);
		memcpy(Y + ip, &joint[i2p + dim], sizeof(double) * dim);
	}
	double beta  = 0.0;
	// The real part of K_start = X %*% t(X) + Y %*% t(Y)
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &X[0], &dim, &X[0], &dim, &beta, &r_sigma_start[0], &dim );
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &Y[0], &dim, &Y[0], &dim, &alpha, &r_sigma_start[0], &dim );

	// The imaginary part of K_start = Y %*% t(X) - X %*% t(Y)
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &alpha, &Y[0], &dim, &X[0], &dim, &beta, &i_sigma_start[0], &dim );
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &n, &malpha, &X[0], &dim, &Y[0], &dim, &alpha, &i_sigma_start[0], &dim );

	for (j = 0; j < dim; j++)
	  for (int k = 0; k < dim; k++){
	    csigma[j*dim+k].r = r_sigma_start[j*dim+k];
	    csigma[j*dim+k].i = i_sigma_start[j*dim+k];
  	  }
	for (j = 0; j < dim; j++)
          for (int k = 0; k < dim; k++){
            Ind[j*dim+k].i = 0;
            Ind[j*dim+k].r = (j == k);
          }
	F77_NAME(zpotrf)(&upper, &dim, csigma, &dim, &info); // Remark: csigma will be changed
	zpotrs(&upper, &dim, &dim, csigma, &dim, Ind, &dim, &info ); //sigma = inv(K)
	
	for (j = 0; j < pxp; j++){
	  r_sigma_start[j] = Ind[j].r;
	  i_sigma_start[j] = Ind[j].i;
	}
	
	memcpy( r_sigma, &r_sigma_start[0], sizeof( double ) * pxp );
	memcpy( i_sigma, &i_sigma_start[0], sizeof( double ) * pxp );	 
	
	double max_diff = 1.0, r_diff, i_diff;	
	while ( max_diff > threshold_C )
	{
		max_diff = 0.0; // The maximum difference between two adjacent sigma
		
		for( i = 0; i < dim; i++ )
		{
			ip = i * dim;

			size_node_i = size_node[i];
			if( size_node_i > 0 )
			{
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					ji = ip + j;
					if( G[ji] )
					{
						r_sigma_start_N_i[l] = r_sigma_start[ji];
						r_sigma_start_N_i_2[l] = r_sigma_start[ji];
						i_sigma_start_N_i[l] = i_sigma_start[ji];
						N_i[l++] = j;
					}
					else
					{
						r_beta_star[j] = 0.0; 
						i_beta_star[j] = 0.0; 
					}
				}
				// -------------------------------------------------------------
				
				sub_matrix( r_sigma, &r_sigma_N_i[0], &N_i[0], &size_node_i, &dim );
				sub_matrix( i_sigma, &i_sigma_N_i[0], &N_i[0], &size_node_i, &dim );
				
				// Inv_R = solve(r_sigma_N_i)
				for (int s = 0; s < size_node_i*size_node_i; s++)
				  r_sigma_N_i_2[s] = r_sigma_N_i[s];
				inverse( &r_sigma_N_i_2[0], &Inv_R[0], &size_node_i );			
				// IR = i_sigma_N_i %*% Inv_R
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &size_node_i, &size_node_i, &alpha, &i_sigma_N_i[0], &size_node_i, &Inv_R[0], &size_node_i, &beta, &IR[0], &size_node_i);
				// A = IR %*% i_sigma_N_i + r_sigma_N_i = r_sigma_N_i
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &size_node_i, &size_node_i, &alpha, &IR[0], &size_node_i, &i_sigma_N_i[0], &size_node_i, &alpha, &r_sigma_N_i[0], &size_node_i);
				// B = IR %*% i_sigma_start_N_i + r_sigma_start_N_i	= r_sigma_start_N_i			
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &one, &size_node_i, &alpha, &IR[0], &size_node_i, &i_sigma_start_N_i[0], &size_node_i, &alpha, &r_sigma_start_N_i[0], &size_node_i);				
				// C = IR %*% r_sigma_start_N_i_2 - i_sigma_start_N_i = i_sigma_start_N_i
				F77_NAME(dgemm)( &transN, &transN, &size_node_i, &one, &size_node_i, &alpha, &IR[0], &size_node_i, &r_sigma_start_N_i_2[0], &size_node_i, &dmone, &i_sigma_start_N_i[0], &size_node_i);			
				
				// r_sigma_start_N_i = solve(A) %*% B; i_sigma_start_N_i = solve(A) %*% C
				for (int s = 0; s < size_node_i*size_node_i; s++)
				  r_sigma_N_i_2[s] = r_sigma_N_i[s];
				F77_NAME(dposv)( &upper, &size_node_i, &one, &r_sigma_N_i[0], &size_node_i, &r_sigma_start_N_i[0], &size_node_i, &info );
				F77_NAME(dposv)( &upper, &size_node_i, &one, &r_sigma_N_i_2[0], &size_node_i, &i_sigma_start_N_i[0], &size_node_i, &info );	

				for( j = 0; j < size_node_i; j++ ) 
				{
					r_beta_star[N_i[j]] = r_sigma_start_N_i[j];
					i_beta_star[N_i[j]] = i_sigma_start_N_i[j];
				}
				
				// r_sigma_i = r_sigma %*% r_beta_star + i_sigma %*% i_beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, i_sigma, &dim, &i_beta_star[0], &dim, &beta, &r_sigma_i[0], &dim );
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, r_sigma, &dim, &r_beta_star[0], &dim, &done, &r_sigma_i[0], &dim );

				// i_sigma_i = i_sigma %*% r_beta_star - r_sigma %*% i_beta_star
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, r_sigma, &dim, &i_beta_star[0], &dim, &beta, &i_sigma_i[0], &dim );
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, i_sigma, &dim, &r_beta_star[0], &dim, &dmone, &i_sigma_i[0], &dim );
				
				// Update the first i elements of sigma[-i,i]
				memcpy( r_sigma + ip, r_sigma_i,  sizeof( double ) * i );	
				memcpy( i_sigma + ip, i_sigma_i,  sizeof( double ) * i );	
				// Update the first i elements of sigma[i,-i]
				for( j = 0; j < i; j++ )
				{
					ij   = j * dim + i;
					r_diff = r_sigma[ij] - r_sigma_i[j];
					i_diff = i_sigma[ij] + i_sigma_i[j];

					temp = sqrt( r_diff * r_diff + i_diff * i_diff );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij] = r_sigma_i[j];
					i_sigma[ij] = - i_sigma_i[j];
				}
				
				memcpy( r_sigma + ip + i + 1, r_sigma_i + i + 1, sizeof( double ) * ( dim - i - 1 ) );	
				memcpy( i_sigma + ip + i + 1, i_sigma_i + i + 1, sizeof( double ) * ( dim - i - 1 ) );	

				for( j = i + 1; j < dim; j++ )
				{
					ij   = j * dim + i;
					r_diff = r_sigma[ij] - r_sigma_i[j];
					i_diff = i_sigma[ij] + i_sigma_i[j];

					temp = sqrt( r_diff * r_diff + i_diff * i_diff );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij] = r_sigma_i[j];
					i_sigma[ij] = - i_sigma_i[j];
				}
			} 
			else 
			{				
				for( j = 0; j < i; j++ )
				{
					ij   = j * dim + i;
					temp = sqrt( r_sigma[ij] * r_sigma[ij] + i_sigma[ij] * i_sigma[ij] );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij]     = 0.0;
					i_sigma[ij]     = 0.0;
					r_sigma[ip + j] = 0.0;
					i_sigma[ip + j] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					ij = j * dim + i;
					temp = sqrt( r_sigma[ij] * r_sigma[ij] + i_sigma[ij] * i_sigma[ij] );
					if( temp > max_diff ) max_diff = temp; 					

					r_sigma[ij]     = 0.0;
					i_sigma[ij]     = 0.0;
					r_sigma[ip + j] = 0.0;				
					i_sigma[ip + j] = 0.0;				
				}
			} 
		}
	}
	
	memcpy( &r_sigma_start[0], r_sigma, sizeof( double ) * pxp );	 	
	memcpy( &i_sigma_start[0], i_sigma, sizeof( double ) * pxp );	 	
	
	// K = solve(sigma)
	for (int j = 0; j < pxp; j++){
	  csigma[j].r = r_sigma_start[j];
	  csigma[j].i = i_sigma_start[j];
	}
	for (int j = 0; j < dim; j++)
		for (int k = 0; k < dim; k++){
		K[j*dim+k].i = 0;
		K[j*dim+k].r = (j == k);
		}

	F77_NAME(zpotrf)(&upper, &dim, csigma, &dim, &info);
	zpotrs(&upper, &dim, &dim, csigma, &dim, K, &dim, &info );
} 
