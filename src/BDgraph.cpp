#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <sstream>
# include <string>       // std::string, std::to_string
#include <vector>        // for using vector

using namespace std;

extern "C" {

// copying square matrix A to matrix copyA, for arrary with one dimention
void copyMatrix( double A[], double copyA[], int *pxp )
{
	for( register int i = 0, dim = *pxp; i < dim; i++ ) copyA[i] = A[i]; 	
}
	
// Takes square matrix A (p x p) and retrieves square submatrix B (p_sub x p_sub), dictated by vector sub
void subMatrix( double A[], double subA[], int sub[], int *p_sub, int *p  )
{
	for( int i = 0, psub = *p_sub, pdim = *p; i < psub; i++ )
		for( register int j = 0; j < psub; j++ )
			subA[j * psub + i] = A[sub[j] * pdim + sub[i]]; 
}

// Takes square matrix A (p x p) and retrieves vector subA which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
void subRowMins( double A[], double subA[], int *sub, int *p )
{
	register int i;
	int dimSub = *sub, pdim = *p;
	
	for( i = 0; i < dimSub; i++ )
		subA[i] = A[i * pdim + dimSub];

	for( i = dimSub + 1; i < pdim; i++ )
		subA[i - 1] = A[i * pdim + dimSub];	
}
   
// Takes square matrix A (p x p) and retrieves submatrix subA(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R
void subRowsMins( double A[], double subA[], int *row, int *col, int *p )
{	
	register int j; 
	int l = 0, pdim = *p, sub0 = *row, sub1 = *col;
	
	for( j = 0; j < sub0; j++ )
	{
		subA[l++] = A[j * pdim + sub0]; 
		subA[l++] = A[j * pdim + sub1]; 
	}
	
	for( j = sub0 + 1; j < sub1; j++ )
	{
		subA[l++] = A[j * pdim + sub0]; 
		subA[l++] = A[j * pdim + sub1]; 
	}

	for( j = sub1 + 1; j < pdim; j++ )
	{
		subA[l++] = A[j * pdim + sub0]; 
		subA[l++] = A[j * pdim + sub1]; 
	}
}

// Takes symmatric matrix A (p x p) and retrieves A_jj, A12(1x(p-1)), A21((p-1)x1), and A22((p-1)x(p-1))
// Like A11=A[j, j], A12=A[j, -j], and A22=A[-j, -j] in R
void subMatrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	register int j;
	int i, pdim = *p, p1 = pdim - 1, psub = *sub;

	for( i = 0; i < psub; i++ )
	{	
		A12[i] = A[i * pdim + psub];

		for( j = 0; j < psub; j++ )
			A22[j * p1 + i] = A[j * pdim + i];

		for( j = psub + 1; j < pdim; j++ )
		{
			A22[(j - 1) * p1 + i] = A[j * pdim + i];
			A22[i * p1 + j - 1]   = A[i * pdim + j];
		}
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		A12[i - 1] = A[i * pdim + psub];
		
		for( j = psub + 1; j < pdim; j++ )
			A22[(j - 1) * p1 + i - 1] = A[j * pdim + i];
	}
}

// Takes square matrix A (p x p) and retrieves A11(2x2), A12(2x(p-2)), and A22((p-2)x(p-2))
// Like A11=A[e, e], A12=A[e, -e], and A22=A[-e, -e] in R
void subMatrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	register int j;
	int i, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;
	A11[0] = A[sub0 * pdim + sub0];
	A11[1] = A[sub0 * pdim + sub1];
	A11[2] = A11[1];                   // for symmetric matrices
	A11[3] = A[sub1 * pdim + sub1];

	for( i = 0; i < sub0; i++ )
	{	
		A12[i * 2]     = A[i * pdim + sub0];
		A12[i * 2 + 1] = A[i * pdim + sub1];
	
		for( j = 0; j < sub0; j++ )
			A22[j * p2 + i] = A[j * pdim + i];

		for( j = sub0 + 1; j < sub1; j++ )
		{
			A22[(j - 1) * p2 + i] = A[j * pdim + i];
			A22[i * p2 + j - 1]   = A[i * pdim + j];
		}
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			A22[(j - 2) * p2 + i] = A[j * pdim + i];
			A22[i * p2 + j - 2]   = A[i * pdim + j];
		}
	}

	for( i = sub0 + 1; i < sub1; i++ )
	{
		A12[(i - 1) * 2] = A[i * pdim + sub0];
		A12[i * 2 - 1]   = A[i * pdim + sub1];
	
		for( j = sub0 + 1; j < sub1; j++ )
			A22[(j - 1) * p2 + i - 1] = A[j * pdim + i];
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			A22[(j - 2) * p2 + i - 1] = A[j * pdim + i];
			A22[(i - 1) * p2 + j - 2] = A[i * pdim + j];
		}
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		A12[(i - 2) * 2]     = A[i * pdim + sub0];
		A12[(i - 2) * 2 + 1] = A[i * pdim + sub1];
		
		for( j = sub1 + 1; j < pdim; j++ )
			A22[(j - 2) * p2 + i - 2] = A[j * pdim + i];
	}
}
   
////////////////////////////////////////////////////////////////////////////////
//  Multiplies (p_i x p_k) matrix by (p_k x p_j) matrix to give (p_i x p_j) matrix
//  C := A %*% B
void multiplyMatrix( double A[], double B[], double C[], int *p_i, int *p_j, int *p_k )
{
	double alpha = 1.0, beta  = 0.0;
	char trans   = 'N';																	
	F77_NAME(dgemm)( &trans, &trans, p_i, p_j, p_k, &alpha, A, p_i, B, p_k, &beta, C, p_i );
}

// inverse function for symmetric positive-definite matrices (p x p)
// ******** WARNING: this function change matrix A **************************
void inverse( double A[], double A_inv[], int *p )
{
	int info, dim = *p;
	char uplo = 'U';

	// creating an identity matrix
	for( int i = 0; i < dim; i++ )
		for(register int j = 0; j < dim; j++ )
			A_inv[j * dim + i] = (i == j);
	
	// LAPACK function: computes solution to A x X = B, where A is symmetric positive definite matrix
	F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info );
}

// inverse function for symmetric (2 x 2)
void inverse2x2( double B[], double B_inv[] )
{
	double detB = B[0] * B[3] - B[1] * B[1];
	B_inv[0]    = B[3] / detB;
	B_inv[1]    = - B[1] / detB;
	B_inv[2]    = B_inv[1];
	B_inv[3]    = B[0] / detB;
}

// sampling from Wishart distribution
// Ts = chol( solve( Ds ) )
void rwish( double Ts[], double K[], int *b, int *p )
{
	int dim = *p, pxp = dim * dim, bK = *b;
	vector<double> psi( pxp ); 

	// ---- Sample values in Psi matrix ---
    GetRNGstate();
	for( int i = 0; i < dim; i++ )
		for( register int j = 0; j < dim; j++ )
			psi[j * dim + i] = (i < j) ? rnorm(0, 1) : ( (i > j) ? 0.0 : sqrt( rchisq( bK + dim - i - 1 ) ) );
	PutRNGstate();
	// ------------------------------------

    // C = psi %*% Ts
    vector<double> C( pxp ); 
	multiplyMatrix( &psi[0], Ts, &C[0], &dim, &dim, &dim );

	// K = t(C) %*% C 
	double alpha = 1.0, beta  = 0.0;
	char transA  = 'T', transB  = 'N';
	// LAPACK function to compute  C := alpha * A * B + beta * C																				
	F77_NAME(dgemm)( &transA, &transB, &dim, &dim, &dim, &alpha, &C[0], &dim, &C[0], &dim, &beta, K, &dim );
}

// G is adjacency matrix which has zero in its diagonal
// threshold = 1e-8
void rgwish( int G[], double Ts[], double K[], int *b, int *p )
{
	register int k;
	int j, l, a, one = 1, dim = *p, pxp = dim * dim;	
	double temp;
	
	rwish( Ts, K, b, &dim );
	
	vector<double> Sigma( pxp ); 
	inverse( K, &Sigma[0], &dim );
	
	// copying  matrix sigma to matrix W	
	vector<double> W( Sigma ); 

	vector<double> W_last( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> ww( dim ); 

	double difference = 1.0;	
	while ( difference > 1e-8 )
	{
		// copying  matrix W to matrix W_last	
		copyMatrix( &W[0], &W_last[0], &pxp ); 	
		
		for( j = 0; j < dim; j++ )
		{
			// Count  size of note
			a = 0;
			for( k = 0; k < dim; k++ ) a += G[k * dim + j];

			if( a > 0 )
			{
				// Record size of node and initialize zero in beta_star for next steps
				vector<double> Sigma_N_j( a );
				vector<int> N_j( a );
				l = 0;
				for( k = 0; k < dim; k++ )
				{
					if( G[k * dim + j] )
					{
						Sigma_N_j[l] = Sigma[j * dim + k]; // Sigma_N_j[k] = Sigma[j * dim + N_j[k]];
						N_j[l++]     = k;
					}
					
					beta_star[k] = 0.0; // for( k = 0; k < *p; k++ ) beta_star[k] = 0.0;
				}
				// -------------------------------------------------------------
				
				vector<double> W_N_j( a * a );
				subMatrix( &W[0], &W_N_j[0], &N_j[0], &a, &dim );
				
				vector<double> W_N_j_inv( a * a );
				inverse( &W_N_j[0], &W_N_j_inv[0], &a );
				
				vector<double> beta_star_hat_j( a );   
				multiplyMatrix( &W_N_j_inv[0], &Sigma_N_j[0], &beta_star_hat_j[0], &a, &one, &a );
				
				for( k = 0; k < a; k++ ) beta_star[N_j[k]] = beta_star_hat_j[k];
				
				multiplyMatrix( &W[0], &beta_star[0], &ww[0], &dim, &one, &dim );
				
				for( k = 0; k < j; k++ )
				{
					W[k * dim + j] = ww[k];
					W[j * dim + k] = ww[k];
				}
				
				for( k = j + 1; k < dim; k++ )
				{
					W[k * dim + j] = ww[k];
					W[j * dim + k] = ww[k];
				}
			} 
			else 
			{
				for( k = 0; k < j; k++ )
				{
					W[k * dim + j] = 0.0;
					W[j * dim + k] = 0.0;
				}
				
				for( k = j + 1; k < dim; k++ )
				{
					W[k * dim + j] = 0.0;
					W[j * dim + k] = 0.0;
				}
			} 
		}

		difference = fabs( W[0] - W_last[0] );
		for( k = 1; k < pxp; k++ )
		{
			temp = fabs( W[k] - W_last[k] );
			if( temp > difference ) difference = temp; 
		}		
	}

	inverse( &W[0], K, &dim );
}
     
// G is adjacency matrix which has zero in its diagonal
// threshold = 1e-8
void rgwish_sigma( int G[], double Ts[], double K[], double W[], int *b, int *p )
{
	register int k;
	int i, j, l, a, one = 1, dim = *p, pxp = dim * dim, bK = *b;	
	double temp;
// STEP 1: sampling from wishart distributions
	vector<double> psi( pxp ); 
	// ---- Sample values in Psi matrix ---
    GetRNGstate();
	for( i = 0; i < dim; i++ )
		for( j = 0; j < dim; j++ )
			psi[j * dim + i] = (i < j) ? rnorm(0, 1) : ( (i > j) ? 0.0 : sqrt( rchisq( bK + dim - i - 1 ) ) );
	PutRNGstate();
	// ------------------------------------

    // C <- psi %*% Ts
    vector<double> C( pxp ); 
	multiplyMatrix( &psi[0], Ts, &C[0], &dim, &dim, &dim );

	vector<double> invC( pxp ); 
	char side = 'L', up = 'U', transA = 'N', diag = 'N';
	double alpha = 1.0;
	// creating an identity matrix
	for( i = 0; i < dim; i++ )
		for( j = 0; j < dim; j++ )
			invC[j * dim + i] = (i == j);
	// op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
	F77_NAME(dtrsm)( &side, &up, &transA, &diag, &dim, &dim, &alpha, &C[0], &dim, &invC[0], &dim );

	vector<double> Sigma( pxp ); 
	// Sigma <- invC %*% t(invC) 
	double beta  = 0.0;
	char transB  = 'T';
	// LAPACK function to compute  C := alpha * A * B + beta * C																				
	F77_NAME(dgemm)( &transA, &transB, &dim, &dim, &dim, &alpha, &invC[0], &dim, &invC[0], &dim, &beta, &Sigma[0], &dim );

	// copying  matrix sigma to matrix W	
	//~ copyMatrix( &Sigma[0], W, &pxp ); 
	for( i = 0; i < pxp; i++ ) W[i] = Sigma[i]; 	 

	vector<double> W_last( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> ww( dim ); 
	
	double difference = 1.0;	
	while ( difference > 1e-8 )
	{
		// copying  matrix W to matrix W_last	
		//~ copyMatrix( W, &W_last[0], &pxp );
		for( i = 0; i < pxp; i++ ) W_last[i] = W[i]; 	 	
		
		for( j = 0; j < dim; j++ )
		{
			// Count  size of note
			a = 0;
			for( k = 0; k < dim; k++ )  a += G[k * dim + j];

			if( a > 0 )
			{
				// Record size of node and initialize zero in beta_star for next steps
				vector<double> Sigma_N_j( a );
				vector<int> N_j( a );
				l = 0;
				for( k = 0; k < dim; k++ )
				{
					if( G[k * dim + j] )
					{
						Sigma_N_j[l] = Sigma[j * dim + k]; // Sigma_N_j[k] = Sigma[j * dim + N_j[k]];
						N_j[l++]     = k;
					}
					
					beta_star[k] = 0.0; // for( k = 0; k < *p; k++ ) beta_star[k] = 0.0;
				}
				// -------------------------------------------------------------
				
				vector<double> W_N_j( a * a );
				subMatrix( W, &W_N_j[0], &N_j[0], &a, &dim );
				
				vector<double> W_N_j_inv( a * a );
				inverse( &W_N_j[0], &W_N_j_inv[0], &a );
				
				vector<double> beta_star_hat_j( a );   
				multiplyMatrix( &W_N_j_inv[0], &Sigma_N_j[0], &beta_star_hat_j[0], &a, &one, &a );
				
				for( k = 0; k < a; k++ ) beta_star[N_j[k]] = beta_star_hat_j[k];
				
				multiplyMatrix( W, &beta_star[0], &ww[0], &dim, &one, &dim );
				
				for( k = 0; k < j; k++ )
				{
					W[k * dim + j] = ww[k];
					W[j * dim + k] = ww[k];
				}
				
				for( k = j + 1; k < dim; k++ )
				{
					W[k * dim + j] = ww[k];
					W[j * dim + k] = ww[k];
				}
			} 
			else 
			{
				for( k = 0; k < j; k++ )
				{
					W[k * dim + j] = 0.0;
					W[j * dim + k] = 0.0;
				}
				
				for( k = j + 1; k < dim; k++ )
				{
					W[k * dim + j] = 0.0;
					W[j * dim + k] = 0.0;
				}
			} 
		}

		difference = fabs( W[0] - W_last[0] );
		for( k = 1; k < pxp; k++ )
		{
			temp = fabs( W[k] - W_last[k] );
			if( temp > difference ) difference = temp; 
		}		
	}
	
	//~ copyMatrix( W, &Sigma[0], &pxp );
	for( i = 0; i < pxp; i++ ) Sigma[i] = W[i]; 	 	
	
	inverse( &Sigma[0], K, &dim );
}
     
////////////////////////////////////////////////////////////////////////////////
// bdmcmc algoirthm with exact value of normalizing constant for D = I_p
////////////////////////////////////////////////////////////////////////////////
void bdmcmcExact( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 char *allGraphs[], double allWeights[], double Ksum[], 
			 char *sampleGraphs[], double graphWeights[], int *sizeSampleG,
			 int lastGraph[], double lastK[],
			 int *b, int *bstar, double Ds[] )
{
	vector<string> allGraphs_C( *iter );
	vector<string> sampleGraphs_C( *iter );
	
	register int col; 
	bool thisOne;

	int selectedEdgei, selectedEdgej, selectedEdgeij, burn_in = *burnin, sizeSampleGraph = *sizeSampleG;
	int row, rowCol, i, j, k, ij, jj, Dsjj, Dsij, counter, nustar, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	double rate, maxRates, sumRates, sumDiag, K022, a11, b1 = *b, sigmaj11;

	//~ vector<double> rates( pxp, 0.0 ); 
	vector<double> sigma( pxp ); 
	copyMatrix( K, lastK, &pxp );
	inverse( lastK, &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> charG( qp ); // char stringG[pp];

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );             // K[j, -j]
	vector<double> sigmaj12( p1 );         // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );      // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 
	vector<double> K12( p2x2 );            // K[e, -e]
	vector<double> sigma11( 4 );           // sigma[e, e]
	vector<double> sigma12( p2x2 );        // sigma[e, -e]
	vector<double> sigma22( p2xp2 );       // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   

	for( int g = 0, iteration = *iter; g < iteration; g++ )
	{
		if( ( g + 1 ) % 1000 == 0 )	Rprintf( " Iteration  %d                 \n", g + 1 ); 
		
// STEP 1: calculating birth-death rates
		// computing birth and death rates
		//~ ratesMatrixExact( K, &sigma[0], &invDsee[0], G, &rates[0], b, Ds, &dim );

		counter  = 0;
		sumRates = 0.0;
		maxRates = -1.7e307; 
		
		for( j = 1; j < dim; j++ )
		{
			jj = j * dim + j;
			Dsjj = Ds[jj];
			
			sigmaj11 = sigma[jj];        // sigma[j, j]  
			subMatrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			for( row = 0; row < p1; row++ )
				for( col = 0; col < p1; col++ )
				{
					rowCol = col * p1 + row;
					Kj22_inv[rowCol] = sigmaj22[rowCol] - sigmaj12[row] * sigmaj12[col] / sigmaj11;
				}
			
			for( i = 0; i < j; i++ )
			{
				ij = j * dim + i;
				Dsij = Ds[ij];

	// For (i,j) = 0 ---------------------------------------------------------------|
				// For (i,j) = 0 
				// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
				//~ K022ij( &K[0], &sigma[0], &K022, &i, &j, &dim );

				subRowMins( K, &Kj12[0], &j, &dim );  // K12 = K[j, -j]  
				Kj12[ i ] = 0.0;                      // K12[1,i] = 0

				// K12 %*% K22_inv
				multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );

				// K121 = K12 %*% solve(K[-e, -e]) %*% t(K12) 
				F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			
	// Finished (i,j) = 0 ----------------------------------------------------------|

	// For (i,j) = 1 ---------------------------------------------------------------|
				// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 
				//~ K121output(  K[],  sigma, K121[],    *i, *j, *p )
				//~ K121output( &K[0], sigma, &K121[0], &i, &j, &dim );
				
				subRowsMins( K, &K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
				
				subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &i, &j, &dim );

				// solve( sigma[e, e] )
				inverse2x2( &sigma11[0], &sigma11_inv[0] );

				// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
				F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma21xsigma11_inv %*% sigma12
				multiplyMatrix( &sigma21xsigma11_inv[0], &sigma12[0], &sigma2112[0], &p2, &p2, &two );

				// solve( K[-e, -e] ) = sigma22 - sigma2112
				for( k = 0; k < p2xp2 ; k++ ) K22_inv[k] = sigma22[k] - sigma2112[k];	
				
				// K12 %*% K22_inv
				multiplyMatrix( &K12[0], &K22_inv[0], &K12xK22_inv[0], &two, &p2, &p2 );
				
				// K121 <- K12 %*% solve(K[-e, -e]) %*% t(K12) 																
				F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
	// Finished (i,j) = 1-----------------------------------------------------------|

				a11 = K[i * dim + i] - K121[0];	
				//~ sumDiagAB( &Dsee[0], &K0_ij[0], &K121[0], &sumDiag );
				//~ sumDiag = Dsii * a11 - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );
				sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

				//	nustar = b + sum( Gf[,i] * Gf[,j] )
				nustar = b1;
				for( k = 0; k < dim; k++ )
					nustar += G[i * dim + k] * G[j * dim + k];

				//~ rates[ij] = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) + ( ( Dsii - Dsij * Dsij / Dsjj ) * a11 - sumDiag ) / 2 );
				rate = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2 );

				if( G[ij] == 0 ) rate = 1.0 / rate;	

				sumRates += rate;        // sumUpperMatrix( &rates[0], &allWeights[g], p );
				
				if( rate > maxRates )    // selectEdge( &rates[0], selectedEdge, p );
				{
					maxRates      = rate; 
					selectedEdgei = i;
					selectedEdgej = j;
				}
				
				charG[counter++] = G[ij] + '0'; // adjToString( G, &allGraphs[g], p );
			}
		}	
		
		allGraphs_C[g]  = std::string( charG.begin(), charG.end() );	
		allWeights[g] = 1 / sumRates;
////////////////////////////////////////////////////////////////////////////////	
		
		if( g > burn_in )
		{
			for( i = 0; i < pxp ; i++ ) Ksum[i] += K[i];	
			
			thisOne = false;
			for( i = 0; i < sizeSampleGraph; i++ )
				if( sampleGraphs_C[i] == allGraphs_C[g] )
				{
					graphWeights[i] += allWeights[g];
					thisOne = true;
					break;
				} 
			
			if( !thisOne || sizeSampleGraph == 0 )
			{
				sampleGraphs_C[sizeSampleGraph] = allGraphs_C[g];
				graphWeights[sizeSampleGraph] = allWeights[g];
				sizeSampleGraph++;				
			} 
		}

		selectedEdgeij    = selectedEdgej * dim + selectedEdgei;
		G[selectedEdgeij] = 1 - G[selectedEdgeij];
		G[selectedEdgei * dim + selectedEdgej] = G[selectedEdgeij];

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim );
	}

	for( i = 0; i < *iter; i++ ) 
	{
		allGraphs_C[i].copy(allGraphs[i], qp, 0);
		allGraphs[i][qp] = '\0';
		sampleGraphs_C[i].copy(sampleGraphs[i], qp, 0);
		sampleGraphs[i][qp] = '\0';
	}
	
	*sizeSampleG = sizeSampleGraph;
	// For last graph and its precision matrix
	for( i = 0; i < pxp; i++ ) 
	{
		lastK[i]     = K[i];	
		lastGraph[i] = G[i];
	}
}
        
////////////////////////////////////////////////////////////////////////////////
// RJMCMC algoirthm with exact value of normalizing constant for D = I_p
////////////////////////////////////////////////////////////////////////////////
void rjmcmcExact( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 char *allGraphs[], double allWeights[], double Ksum[], 
			 char *sampleGraphs[], double graphWeights[], int *sizeSampleG,
			 int lastGraph[], double lastK[],
			 int *b, int *bstar, double Ds[] )
{
	vector<string> allGraphs_C( *iter );
	vector<string> sampleGraphs_C( *iter );

	register int col; 
	bool thisOne;

	int randomEdge, selectedEdgei, selectedEdgej, burn_in = *burnin, sizeSampleGraph = *sizeSampleG;
	int row, rowCol, i, j, k, ij, jj, Dsjj, Dsij, counter, nustar, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	double sumDiag, K022, a11, b1 = *b, sigmaj11;

	vector<double> sigma( pxp ); 
	copyMatrix( K, lastK, &pxp );
	inverse( lastK, &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> charG( qp ); // char stringG[pp];

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );             // K[j, -j]
	vector<double> sigmaj12( p1 );         // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );      // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 
	vector<double> K12( p2x2 );            // K[e, -e]
	vector<double> sigma11( 4 );           // sigma[e, e]
	vector<double> sigma12( p2x2 );        // sigma[e, -e]
	vector<double> sigma22( p2xp2 );       // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   

	double alpha_ij;
	GetRNGstate();
	for( int g = 0, iteration = *iter; g < iteration; g++ )
	{
		if( ( g + 1 ) % 1000 == 0 )	Rprintf( " Iteration  %d                 \n", g + 1 ); 
		
// STEP 1: selecting random edge and calculating alpha

		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		randomEdge = static_cast<int>( runif( 0, 1 ) * qp );

		counter = 0;
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				if( counter == randomEdge )    
				{
					selectedEdgei = i;
					selectedEdgej = j;
				}
				
				charG[counter++] = G[j * dim + i] + '0'; // adjToString( G, &allGraphs[g], p );
			}
		
		// -------- Calculating alpha -------------------------------------|
		ij   = selectedEdgej * dim + selectedEdgei;
		jj   = selectedEdgej * dim + selectedEdgej;
		Dsij = Ds[ij];
		Dsjj = Ds[jj];
		
		sigmaj11 = sigma[jj];        // sigma[j, j]  
		subMatrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &selectedEdgej, &dim );

		// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
		for( row = 0; row < p1; row++ )
			for( col = 0; col < p1; col++ )
			{
				rowCol = col * p1 + row;
				Kj22_inv[rowCol] = sigmaj22[rowCol] - sigmaj12[row] * sigmaj12[col] / sigmaj11;
			}

// For (i,j) = 0 ---------------------------------------------------------------|
		// For (i,j) = 0 
		// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
		//~ K022ij( &K[0], &sigma[0], &K022, &i, &j, &dim );

		subRowMins( K, &Kj12[0], &selectedEdgej, &dim );  // K12 = K[j, -j]  
		Kj12[ selectedEdgei ] = 0.0;                      // K12[1,i] = 0

		// K12 %*% K22_inv
		multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );

		// K121 = K12 %*% solve(K[-e, -e]) %*% t(K12) 
		F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			
// Finished (i,j) = 0 ----------------------------------------------------------|

// For (i,j) = 1 ---------------------------------------------------------------|
		// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 
		//~ K121output(  K[],  sigma, K121[],    *i, *j, *p )
		//~ K121output( &K[0], sigma, &K121[0], &i, &j, &dim );
		
		subRowsMins( K, &K12[0], &selectedEdgei, &selectedEdgej, &dim );  // K12 = K[e, -e]  
		
		subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &selectedEdgei, &selectedEdgej, &dim );

		// solve( sigma[e, e] )
		inverse2x2( &sigma11[0], &sigma11_inv[0] );

		// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
		F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

		// sigma21xsigma11_inv %*% sigma12
		multiplyMatrix( &sigma21xsigma11_inv[0], &sigma12[0], &sigma2112[0], &p2, &p2, &two );

		// solve( K[-e, -e] ) = sigma22 - sigma2112
		for( k = 0; k < p2xp2 ; k++ ) K22_inv[k] = sigma22[k] - sigma2112[k];	
		
		// K12 %*% K22_inv
		multiplyMatrix( &K12[0], &K22_inv[0], &K12xK22_inv[0], &two, &p2, &p2 );
		
		// K121 <- K12 %*% solve(K[-e, -e]) %*% t(K12) 																
		F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
// Finished (i,j) = 1-----------------------------------------------------------|

		a11 = K[selectedEdgei * dim + selectedEdgei] - K121[0];	
		//~ sumDiagAB( &Dsee[0], &K0_ij[0], &K121[0], &sumDiag );
		//~ sumDiag = Dsii * a11 - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );
		sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

		//	nustar = b + sum( Gf[,i] * Gf[,j] )
		nustar = b1;
		for( k = 0; k < dim; k++ )
			nustar += G[selectedEdgei * dim + k] * G[selectedEdgej * dim + k];

		//~ alpha_ij[ij] = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) + ( ( Dsii - Dsij * Dsij / Dsjj ) * a11 - sumDiag ) / 2 );
		//~ alpha_ij = log( sqrt( 2.0 * Dsjj / a11 ) ) + lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2;
		alpha_ij = ( log(2.0) + log(Dsjj) - log(a11) ) / 2 + lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2;

		if( G[ij] == 0 ) alpha_ij = - alpha_ij;	
		// -------- End calculating alpha -----------------------------|
		
		allGraphs_C[g]  = std::string( charG.begin(), charG.end() );	
		allWeights[g] = 1.0;
////////////////////////////////////////////////////////////////////////////////	
		
		if( g > burn_in )
		{
			for( i = 0; i < pxp ; i++ ) Ksum[i] += K[i];	
			
			thisOne = false;
			for( i = 0; i < sizeSampleGraph; i++ )
				if( sampleGraphs_C[i] == allGraphs_C[g] )
				{
					graphWeights[i] += allWeights[g];
					thisOne = true;
					break;
				} 
			
			if( !thisOne || sizeSampleGraph == 0 )
			{
				sampleGraphs_C[sizeSampleGraph] = allGraphs_C[g];
				graphWeights[sizeSampleGraph] = allWeights[g];
				sizeSampleGraph++;				
			} 
		}
  		
		if( log( runif( 0, 1 ) ) < alpha_ij )
		{
			G[ij] = 1 - G[ij];
			G[selectedEdgei * dim + selectedEdgej] = G[ij];
		}

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim );
	}
	PutRNGstate();

	for( i = 0; i < *iter; i++ ) 
	{
		allGraphs_C[i].copy(allGraphs[i], qp, 0);
		allGraphs[i][qp] = '\0';
		sampleGraphs_C[i].copy(sampleGraphs[i], qp, 0);
		sampleGraphs[i][qp] = '\0';
	}
	
	*sizeSampleG = sizeSampleGraph;
	// For last graph and its precision matrix
	for( i = 0; i < pxp; i++ ) 
	{
		lastK[i]     = K[i];	
		lastGraph[i] = G[i];
	}
}
        
////////////////////////////////////////////////////////////////////////////////
// Gaussian copula graphical models
////////////////////////////////////////////////////////////////////////////////
// for copula function
void getMean( double Z[], double K[], double *muij, double *sigma, int *i, int *j, int *n, int *p )
{
	register int k;
	int dim = *p, number = *n, row = *i, col = *j;
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

	for( register int k = 0, number = *n; k < number; k++ )
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
		for( register int i = 0; i < number; i++ )
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

	for( register int k = 0, number = *n; k < number; k++ )
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
		for( register int i = 0; i < number; i++ )
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
	for( register int i = 0; i < pxp ; i++ ) Ds[i] = D[i] + S[i];		
}

// Cholesky decomposition of symmetric positive-definite matrix
// Note: the matrix you pass this function is overwritten with the result.
// A = U' %*% U
void cholesky( double A[], double U[], int *p )
{
	char uplo = 'U';
	register int j;
	int info, dim = *p, i, pxp = dim * dim;
	// copying  matrix A to matrix L	
	for( i = 0; i < pxp; i++ ) U[i] = A[i]; 	
	
	F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info );	

	// dpotrf function doesn't zero out the lower diagonal
	for( i = 0; i < dim; i++ )
		for( j = 0; j < i; j++ )
			U[j * dim + i] = 0.0;
}

// for bdmcmcCopulaDmh function
void getTs( double Ds[], double Ts[], int *p )
{
	int dim = *p, pxp = dim * dim;
	vector<double> invDs( pxp ); 
	vector<double> copyDs( pxp ); 

	//~ copyMatrix( Ds, &copyDs[0], &pxp );
	for( register int i = 0; i < pxp; i++ ) copyDs[i] = Ds[i]; 	
	
	inverse( &copyDs[0], &invDs[0], &dim );	

	cholesky( &invDs[0], Ts, &dim );	
}

////////////////////////////////////////////////////////////////////////////////
// Gaussian copula graphical models 
// Based on reversible jump MCMC algorithm
//********* with NEW idea of exact normalizing constant ***********************
////////////////////////////////////////////////////////////////////////////////
void rjmcmcCopula( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 char *allGraphs[], double allWeights[], double Ksum[], 
			 char *sampleGraphs[], double graphWeights[], int *sizeSampleG,
			 int lastGraph[], double lastK[],
			 int *b, int *bstar, double D[], double Ds[] )
{
	vector<string> allGraphs_C( *iter );
	vector<string> sampleGraphs_C( *iter );	
	
	int randomEdge, counter, selectedEdgei, selectedEdgej, burn_in = *burnin, sizeSampleGraph = *sizeSampleG;
	bool thisOne;

	register int col; 
	double Dsjj, Dsij, sumDiag, K022, a11, b1 = *b, sigmaj11;
	int row, rowCol, i, j, k, ij, jj, nustar, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	//~ vector<double> rates( pxp, 0.0 ); 
	vector<double> sigma( pxp ); 
	copyMatrix( K, &lastK[0], &pxp );
	inverse( &lastK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> charG( qp ); // char stringG[pp];
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	vector<double> Kj12( p1 );             // K[j, -j]
	vector<double> sigmaj12( p1 );         // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );      // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );            // K[e, -e]
	vector<double> sigma11( 4 );           // sigma[e, e]
	vector<double> sigma12( p2x2 );        // sigma[e, -e]
	vector<double> sigma22( p2xp2 );       // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   

	double alpha_ij;
	GetRNGstate();
	for( int g = 0, iteration = *iter; g < iteration; g++ )
	{
		if( ( g + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                  \n", g + 1 );

		getDs( K, Z, R, D, Ds, gcgm, n, &dim );

		getTs( Ds, Ts, &dim );
		
// STEP 1: selecting random edge and calculating alpha

		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		randomEdge = static_cast<int>( runif( 0, 1 ) * qp );

		counter   = 0;
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				if( counter == randomEdge )    
				{
					selectedEdgei = i;
					selectedEdgej = j;
				}
				
				charG[counter++] = G[j * dim + i] + '0'; // adjToString( G, &allGraphs[g], p );
			}
		
		// -------- Calculating alpha -------------------------------------|
		ij   = selectedEdgej * dim + selectedEdgei;
		jj   = selectedEdgej * dim + selectedEdgej;
		Dsij = Ds[ij];
		Dsjj = Ds[jj];
		
		sigmaj11 = sigma[jj];        // sigma[j, j]  
		subMatrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &selectedEdgej, &dim );

		// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
		for( row = 0; row < p1; row++ )
			for( col = 0; col < p1; col++ )
			{
				rowCol = col * p1 + row;
				Kj22_inv[rowCol] = sigmaj22[rowCol] - sigmaj12[row] * sigmaj12[col] / sigmaj11;
			}

// For (i,j) = 0 ---------------------------------------------------------------|
		// For (i,j) = 0 
		// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
		//~ K022ij( &K[0], &sigma[0], &K022, &i, &j, &dim );

		subRowMins( K, &Kj12[0], &selectedEdgej, &dim );  // K12 = K[j, -j]  
		Kj12[ selectedEdgei ] = 0.0;                      // K12[1,i] = 0

		// K12 %*% K22_inv
		multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );

		// K121 = K12 %*% solve(K[-e, -e]) %*% t(K12) 
		F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			
// Finished (i,j) = 0 ----------------------------------------------------------|

// For (i,j) = 1 ---------------------------------------------------------------|
		// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 
		//~ K121output(  K[],  sigma, K121[],    *i, *j, *p )
		//~ K121output( &K[0], sigma, &K121[0], &i, &j, &dim );
		
		subRowsMins( K, &K12[0], &selectedEdgei, &selectedEdgej, &dim );  // K12 = K[e, -e]  
		
		subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &selectedEdgei, &selectedEdgej, &dim );

		// solve( sigma[e, e] )
		inverse2x2( &sigma11[0], &sigma11_inv[0] );

		// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
		F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

		// sigma21xsigma11_inv %*% sigma12
		multiplyMatrix( &sigma21xsigma11_inv[0], &sigma12[0], &sigma2112[0], &p2, &p2, &two );

		// solve( K[-e, -e] ) = sigma22 - sigma2112
		for( k = 0; k < p2xp2 ; k++ ) K22_inv[k] = sigma22[k] - sigma2112[k];	
		
		// K12 %*% K22_inv
		multiplyMatrix( &K12[0], &K22_inv[0], &K12xK22_inv[0], &two, &p2, &p2 );
		
		// K121 <- K12 %*% solve(K[-e, -e]) %*% t(K12) 																
		F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
// Finished (i,j) = 1-----------------------------------------------------------|

		a11 = K[selectedEdgei * dim + selectedEdgei] - K121[0];	
		//~ sumDiagAB( &Dsee[0], &K0_ij[0], &K121[0], &sumDiag );
		//~ sumDiag = Dsii * a11 - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );
		sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

		//	nustar = b + sum( Gf[,i] * Gf[,j] )
		nustar = b1;
		for( k = 0; k < dim; k++ )
			nustar += G[selectedEdgei * dim + k] * G[selectedEdgej * dim + k];

		//~ alpha_ij[ij] = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) + ( ( Dsii - Dsij * Dsij / Dsjj ) * a11 - sumDiag ) / 2 );
		//~ alpha_ij = log( sqrt( 2.0 * Dsjj / a11 ) ) + lgamma( ( nustar + 1 ) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2;
		alpha_ij = ( log(2.0) + log(Dsjj) - log(a11) ) / 2 + lgamma( ( nustar + 1 ) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2;

		if( G[ij] == 0 ) alpha_ij = - alpha_ij;	
		// -------- End calculating alpha -----------------------------|
		
		allGraphs_C[g]  = std::string( charG.begin(), charG.end() );	
		allWeights[g] = 1;
////////////////////////////////////////////////////////////////////////////////	
		
		if( g > burn_in )
		{
			for( i = 0; i < pxp ; i++ ) Ksum[i] += K[i];	
			
			thisOne = false;
			for( i = 0; i < sizeSampleGraph; i++ )
				if( sampleGraphs_C[i] == allGraphs_C[g] )
				{
					graphWeights[i] += allWeights[g];
					thisOne = true;
					break;
				} 
			
			if( !thisOne || sizeSampleGraph == 0 )
			{
				sampleGraphs_C[sizeSampleGraph] = allGraphs_C[g];
				graphWeights[sizeSampleGraph] = allWeights[g];
				sizeSampleGraph++;				
			} 
		}

		if( log( runif( 0, 1 ) ) < alpha_ij )
		{
			G[ij] = 1 - G[ij];
			G[selectedEdgei * dim + selectedEdgej] = G[ij];
		}

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim );
	}
	PutRNGstate();

	for( i = 0; i < *iter; i++ ) 
	{
		allGraphs_C[i].copy(allGraphs[i], qp, 0);
		allGraphs[i][qp] = '\0';
		sampleGraphs_C[i].copy(sampleGraphs[i], qp, 0);
		sampleGraphs[i][qp] = '\0';
	}
	
	*sizeSampleG = sizeSampleGraph;
	// For last graph and its precision matrix
	for( i = 0; i < pxp; i++ ) 
	{
		lastK[i]     = K[i];	
		lastGraph[i] = G[i];
	}
}
    
////////////////////////////////////////////////////////////////////////////////
// Gaussian copula graphical models 
// based on birth-death MCMC algorithm ***********************
////////////////////////////////////////////////////////////////////////////////
void bdmcmcCopula( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 char *allGraphs[], double allWeights[], double Ksum[], 
			 char *sampleGraphs[], double graphWeights[], int *sizeSampleG,
			 int lastGraph[], double lastK[],
			 int *b, int *bstar, double D[], double Ds[] )
{
	vector<string> allGraphs_C( *iter );
	vector<string> sampleGraphs_C( *iter );
	
	int selectedEdgei, selectedEdgej, selectedEdgeij, l, burn_in = *burnin, sizeSampleGraph = *sizeSampleG;
	bool thisOne;

	register int col; 
	double Dsjj, Dsij, sumDiag, K022, a11, b1 = *b, sigmaj11;
	int row, rowCol, i, j, k, ij, jj, nustar, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	//~ vector<double> rates( pxp, 0.0 ); 
	vector<double> sigma( pxp ); 
	copyMatrix( K, &lastK[0], &pxp );
	inverse( &lastK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> charG( qp ); // char stringG[pp];
	double rate, maxRates, sumRates; 
	
	vector<double> K121( 4 ); 

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	vector<double> Kj12( p1 );             // K[j, -j]
	vector<double> sigmaj12( p1 );         // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );      // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K12( p2x2 );            // K[e, -e]
	vector<double> sigma11( 4 );           // sigma[e, e]
	vector<double> sigma12( p2x2 );        // sigma[e, -e]
	vector<double> sigma22( p2xp2 );       // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   

	for( int g = 0, iteration = *iter; g < iteration; g++ )
	{
		if( ( g + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                  \n", g + 1 );

		getDs( K, Z, R, D, Ds, gcgm, n, &dim );

		getTs( Ds, Ts, &dim );
		
		l = 0;
		sumRates = 0.0;
		maxRates = -1.7e307; 

		// computing birth and death rates
		//~ ratesMatrixExactCopula( K, &sigma[0], G, &rates[0], b, Ds, &dim );
		for( j = 1; j < dim; j++ )
		{
			jj   = j * dim + j;
			Dsjj = Ds[jj];

			sigmaj11 = sigma[jj];        // sigma[j, j]  
			subMatrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			for( row = 0; row < p1; row++ )
				for( col = 0; col < p1; col++ )
				{
					rowCol = col * p1 + row;
					Kj22_inv[rowCol] = sigmaj22[rowCol] - sigmaj12[row] * sigmaj12[col] / sigmaj11;
				}
			
			for( i = 0; i < j; i++ )
			{
				ij   = j * dim + i;
				Dsij = Ds[ij];

	// For (i,j) = 0 ---------------------------------------------------------------|
				// For (i,j) = 0 
				// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
				//~ K022ij( &K[0], &sigma[0], &K022, &i, &j, &dim );

				subRowMins( K, &Kj12[0], &j, &dim );  // K12 = K[j, -j]  
				Kj12[ i ] = 0.0;                      // K12[1,i] = 0

				// K12 %*% K22_inv
				multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );

				// K121 = K12 %*% solve(K[-e, -e]) %*% t(K12) 
				F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			
	// Finished (i,j) = 0 ----------------------------------------------------------|

	// For (i,j) = 1 ---------------------------------------------------------------|
				// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 
				//~ K121output(  K[],  sigma, K121[],    *i, *j, *p )
				//~ K121output( &K[0], sigma, &K121[0], &i, &j, &dim );
				
				subRowsMins( K, &K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
				
				subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &i, &j, &dim );

				// solve( sigma[e, e] )
				inverse2x2( &sigma11[0], &sigma11_inv[0] );

				// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
				F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma21xsigma11_inv %*% sigma12
				multiplyMatrix( &sigma21xsigma11_inv[0], &sigma12[0], &sigma2112[0], &p2, &p2, &two );

				// solve( K[-e, -e] ) = sigma22 - sigma2112
				for( k = 0; k < p2xp2 ; k++ ) K22_inv[k] = sigma22[k] - sigma2112[k];	
				
				// K12 %*% K22_inv
				multiplyMatrix( &K12[0], &K22_inv[0], &K12xK22_inv[0], &two, &p2, &p2 );
				
				// K121 <- K12 %*% solve(K[-e, -e]) %*% t(K12) 																
				F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
	// Finished (i,j) = 1-----------------------------------------------------------|

				a11 = K[i * dim + i] - K121[0];	
				//~ sumDiagAB( &Dsee[0], &K0_ij[0], &K121[0], &sumDiag );
				//~ sumDiag = Dsii * a11 - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );
				sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

				//	nustar = b + sum( Gf[,i] * Gf[,j] )
				nustar = b1;
				for( k = 0; k < dim; k++ )
					nustar += G[i * dim + k] * G[j * dim + k];

				//~ rates[ij] = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) + ( ( Dsii - Dsij * Dsij / Dsjj ) * a11 - sumDiag ) / 2 );
				rate = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2 );

				if( G[ij] == 0 ) rate = 1.0 / rate;		

				sumRates += rate;        // sumUpperMatrix( &rates[0], &allWeights[g], p );
				
				if( rate > maxRates )    // selectEdge( &rates[0], selectedEdge, p );
				{
					maxRates      = rate; 
					selectedEdgei = i;
					selectedEdgej = j;
				}
				
				charG[l++] = G[ij] + '0'; // adjToString( G, &allGraphs[g], p );
			}
		}
		
		allGraphs_C[g]  = std::string( charG.begin(), charG.end() );	
		allWeights[g] = 1 / sumRates;
////////////////////////////////////////////////////////////////////////////////	
		
		if( g > burn_in )
		{
			for( i = 0; i < pxp ; i++ ) Ksum[i] += K[i];	
			
			thisOne = false;
			for( i = 0; i < sizeSampleGraph; i++ )
				if( sampleGraphs_C[i] == allGraphs_C[g] )
				{
					graphWeights[i] += allWeights[g];
					thisOne = true;
					break;
				} 
			
			if( !thisOne || sizeSampleGraph == 0 )
			{
				sampleGraphs_C[sizeSampleGraph] = allGraphs_C[g];
				graphWeights[sizeSampleGraph] = allWeights[g];
				sizeSampleGraph++;				
			} 
		}

		selectedEdgeij    = selectedEdgej * dim + selectedEdgei;
		G[selectedEdgeij] = 1 - G[selectedEdgeij];
		G[selectedEdgei * dim + selectedEdgej] = G[selectedEdgeij];

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim );
	}

	for( i = 0; i < *iter; i++ ) 
	{
		allGraphs_C[i].copy(allGraphs[i], qp, 0);
		allGraphs[i][qp] = '\0';
		sampleGraphs_C[i].copy(sampleGraphs[i], qp, 0);
		sampleGraphs[i][qp] = '\0';
	}
	
	*sizeSampleG = sizeSampleGraph;
	// For last graph and its precision matrix
	for( i = 0; i < pxp; i++ ) 
	{
		lastK[i]     = K[i];	
		lastGraph[i] = G[i];
	}
}
    
////////////////////////////////////////////////////////////////////////////////
// For generating scale-free graphs: matrix G (p x p) is an adjacency matrix
void scaleFree( int *G, int *p )
{
    int i, j, tmp, total, dim = *p, p0 = 2;
    double randomValue;
    vector<int> size_a( dim ); 

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


