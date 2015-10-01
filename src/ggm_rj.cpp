#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include "matrix.h"
#include "rgwish.h"

using namespace std;

extern "C" {

////////////////////////////////////////////////////////////////////////////////
// RJMCMC algoirthm with exact value of normalizing constant for D = I_p
////////////////////////////////////////////////////////////////////////////////
void rjmcmcExact( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 int allGraphs[], double allWeights[], double Khat[], 
			 char *sampleGraphs[], double graphWeights[], int *sizeSampleG,
			 int *b, int *bstar, double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin;
	int counterallG = 0;
	string stringG;
	vector<string> sampleGraphs_C( iteration - burn_in );

	bool thisOne;

	int randomEdge, selectedEdgei, selectedEdgej, sizeSampleGraph = *sizeSampleG;
	int row, col, rowCol, i, j, k, ij, jj, Dsjj, Dsij, counter, nustar, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	double sumDiag, K022, a11, b1 = *b, sigmaj11, threshold_C = *threshold;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> charG( qp );              // char stringG[pp];

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
// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> invC( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
// ----------------------------

	double alpha_ij;
	GetRNGstate();
	for( int g = 0; g < iteration; g++ )
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
				
				charG[counter++] = G[j * dim + i] + '0'; 
			}
		
		// -------- Calculating alpha -----------------------------------------|
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
		subRowMins( K, &Kj12[0], &selectedEdgej, &dim );  // K12 = K[j, -j]  
		Kj12[ selectedEdgei ] = 0.0;                      // K12[1,i] = 0

		// K12 %*% K22_inv
		//multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );															
		F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &Kj12[0], &one, &Kj22_inv[0], &p1, &beta, &Kj12xK22_inv[0], &one );

		// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12) 
		F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			

// For (i,j) = 1 ---------------------------------------------------------------|
		subRowsMins( K, &K12[0], &selectedEdgei, &selectedEdgej, &dim );  // K12 = K[e, -e]  
		
		subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &selectedEdgei, &selectedEdgej, &dim );

		// solve( sigma[e, e] )
		inverse2x2( &sigma11[0], &sigma11_inv[0] );

		// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
		F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

		// sigma21xsigma11_inv %*% sigma12
		//multiplyMatrix( &sigma21xsigma11_inv[0], &sigma12[0], &sigma2112[0], &p2, &p2, &two );
		F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &sigma21xsigma11_inv[0], &p2, &sigma12[0], &two, &beta, &sigma2112[0], &p2 );

		// solve( K[-e, -e] ) = sigma22 - sigma2112
		for( k = 0; k < p2xp2 ; k++ ) K22_inv[k] = sigma22[k] - sigma2112[k];	
		
		// K12 %*% K22_inv
		// multiplyMatrix( &K12[0], &K22_inv[0], &K12xK22_inv[0], &two, &p2, &p2 );
		F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &K12[0], &two, &K22_inv[0], &p2, &beta, &K12xK22_inv[0], &two );

		// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e])															
		F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
// Finished (i,j) = 1-----------------------------------------------------------|

		a11     = K[selectedEdgei * dim + selectedEdgei] - K121[0];	
		sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

		//	nustar = b + sum( Gf[,i] * Gf[,j] )
		nustar = b1;
		for( k = 0; k < dim; k++ )
			nustar += G[selectedEdgei * dim + k] * G[selectedEdgej * dim + k];

		alpha_ij = ( log( static_cast<double>(2.0) ) + log( static_cast<double>(Dsjj) ) - log( static_cast<double>(a11) ) ) / 2 + 
		          lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2;

		if( G[ij] == 0 ) alpha_ij = - alpha_ij;	
		// -------- End calculating alpha -------------------------------------|
		
// saving result --------------------------------------------------------------|	
		if( g >= burn_in )
		{
			for( i = 0; i < pxp ; i++ ) Khat[i] += K[i];	

			stringG = string( charG.begin(), charG.end() );	
			
			thisOne = false;
			for( i = 0; i < sizeSampleGraph; i++ )
				if( sampleGraphs_C[i] == stringG )
				{
					graphWeights[i]++;           // += allWeights[counterallG];
					allGraphs[counterallG] = i;
					thisOne = true;
					break;
				} 
			
			if( !thisOne || sizeSampleGraph == 0 )
			{
				sampleGraphs_C[sizeSampleGraph] = stringG;
				graphWeights[sizeSampleGraph]   = allWeights[counterallG];
				allGraphs[counterallG]          = sizeSampleGraph; 
				sizeSampleGraph++;				
			}
			
			counterallG++; 
		}
// Finished saving result -----------------------------------------------------|	
  		
		if( log( static_cast<double>( runif( 0, 1 ) ) ) < alpha_ij )
		{
			G[ij] = 1 - G[ij];
			G[selectedEdgei * dim + selectedEdgej] = G[ij];
		}

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim, &threshold_C, &sigma_start[0], &invC[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}
	PutRNGstate();

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sampleGraphs_C[i].copy(sampleGraphs[i], qp, 0);
		sampleGraphs[i][qp] = '\0';
	}
	
	*sizeSampleG = sizeSampleGraph;
}

////////////////////////////////////////////////////////////////////////////////
// RJMCMC algoirthm with exact value of normalizing constant for D = I_p
////////////////////////////////////////////////////////////////////////////////
void rjmcmcExactPhat( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double Khat[], int phat[],
			 int *b, int *bstar, double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin;

	int randomEdge, selectedEdgei, selectedEdgej;
	int row, col, rowCol, i, j, k, ij, jj, Dsjj, Dsij, counter, nustar, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	double sumDiag, K022, a11, b1 = *b, sigmaj11, threshold_C = *threshold;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;

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
// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> invC( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
// ----------------------------

	double alpha_ij;
	GetRNGstate();
	for( int g = 0; g < iteration; g++ )
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
				
				counter++;
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
		subRowMins( K, &Kj12[0], &selectedEdgej, &dim );  // K12 = K[j, -j]  
		Kj12[ selectedEdgei ] = 0.0;                      // K12[1,i] = 0

		// K12 %*% K22_inv
		//multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );															
		F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &Kj12[0], &one, &Kj22_inv[0], &p1, &beta, &Kj12xK22_inv[0], &one );


		// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
		F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			

// For (i,j) = 1 ---------------------------------------------------------------|
		subRowsMins( K, &K12[0], &selectedEdgei, &selectedEdgej, &dim );  // K12 = K[e, -e]  
		
		subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &selectedEdgei, &selectedEdgej, &dim );

		// solve( sigma[e, e] )
		inverse2x2( &sigma11[0], &sigma11_inv[0] );

		// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
		F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

		// sigma21xsigma11_inv %*% sigma12
		//multiplyMatrix( &sigma21xsigma11_inv[0], &sigma12[0], &sigma2112[0], &p2, &p2, &two );
		F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &sigma21xsigma11_inv[0], &p2, &sigma12[0], &two, &beta, &sigma2112[0], &p2 );

		// solve( K[-e, -e] ) = sigma22 - sigma2112
		for( k = 0; k < p2xp2 ; k++ ) K22_inv[k] = sigma22[k] - sigma2112[k];	
		
		// K12 %*% K22_inv
		// multiplyMatrix( &K12[0], &K22_inv[0], &K12xK22_inv[0], &two, &p2, &p2 );
		F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &K12[0], &two, &K22_inv[0], &p2, &beta, &K12xK22_inv[0], &two );
		
		// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 															
		F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
// Finished (i,j) = 1-----------------------------------------------------------|

		a11     = K[selectedEdgei * dim + selectedEdgei] - K121[0];	
		sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

		//	nustar = b + sum( Gf[,i] * Gf[,j] )
		nustar = b1;
		for( k = 0; k < dim; k++ )
			nustar += G[selectedEdgei * dim + k] * G[selectedEdgej * dim + k];

		alpha_ij = ( log( static_cast<double>(2.0) ) + log( static_cast<double>(Dsjj) ) - log( static_cast<double>(a11) ) ) / 2 + 
		          lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2;

		if( G[ij] == 0 ) alpha_ij = - alpha_ij;	
		// -------- End calculating alpha -----------------------------|
		
// saving result --------------------------------------------------------------|	
		if( g >= burn_in )
			for( i = 0; i < pxp ; i++ )
			{
				Khat[i] += K[i];
				phat[i] += G[i];
			}	
// Finished saving result -----------------------------------------------------|	
  		
		if( log( static_cast<double>( runif( 0, 1 ) ) ) < alpha_ij )
		{
			G[ij] = 1 - G[ij];
			G[selectedEdgei * dim + selectedEdgej] = G[ij];
		}

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim, &threshold_C, &sigma_start[0], &invC[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}
	PutRNGstate();
}
           
} // exturn "C"
