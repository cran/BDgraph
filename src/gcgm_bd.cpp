#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include "matrix.h"
#include "rgwish.h"
#include "copula.h"

using namespace std;

extern "C" {
// Gaussian copula graphical models 
// based on birth-death MCMC algorithm ***********************
 
void bdmcmcCopula( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 int allGraphs[], double allWeights[], double Khat[], 
			 char *sampleGraphs[], double graphWeights[], int *sizeSampleG,
			 int *b, int *bstar, double D[], double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin;
	int counterallG = 0;
	string stringG;
	vector<string> sampleGraphs_C( iteration - burn_in );
	
	int selectedEdgei, selectedEdgej, selectedEdgeij, l, sizeSampleGraph = *sizeSampleG;
	bool thisOne;

	double Dsjj, Dsij, sumDiag, K022, a11, b1 = *b, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, nustar, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> charG( qp );              // char stringG[pp];
	double sumWeights = 0.0, rate, maxRates, sumRates; 
	
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
// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> invC( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
// ----------------------------

	for( int g = 0; g < iteration; g++ )
	{
		if( ( g + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                  \n", g + 1 );

		getDs( K, Z, R, D, Ds, gcgm, n, &dim );

		getTs( Ds, Ts, &dim );
		
		l = 0;
		sumRates = 0.0;
		maxRates = -1.7e307; 

		// computing birth and death rates
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
				subRowMins( K, &Kj12[0], &j, &dim );  // K12 = K[j, -j]  
				Kj12[ i ] = 0.0;                      // K12[1,i] = 0

				// K12 %*% K22_inv
				//multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );															
				F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &Kj12[0], &one, &Kj22_inv[0], &p1, &beta, &Kj12xK22_inv[0], &one );
				
				// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
				F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			

	// For (i,j) = 1 ---------------------------------------------------------------|
				subRowsMins( K, &K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
				
				subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &i, &j, &dim );

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

				a11     = K[i * dim + i] - K121[0];	
				sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

				//	nustar = b + sum( Gf[,i] * Gf[,j] )
				nustar = b1;
				for( k = 0; k < dim; k++ )
					nustar += G[i * dim + k] * G[j * dim + k];

				rate = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2 );

				if( G[ij] == 0 ) rate = 1.0 / rate;		

				sumRates += rate;        // sumUpperMatrix( &rates[0], &allWeights[g], p );
				
				if( rate > maxRates )    // selectEdge( &rates[0], selectedEdge, p );
				{
					maxRates      = rate; 
					selectedEdgei = i;
					selectedEdgej = j;
				}
				
				charG[l++] = G[ij] + '0'; 
			}
		}
		
// saving result --------------------------------------------------------------|	
		if( g >= burn_in )
		{
			for( i = 0; i < pxp ; i++ ) Khat[i] += K[i] / sumRates;	

			stringG = string( charG.begin(), charG.end() );	
			allWeights[counterallG] = 1 / sumRates;
			
			thisOne = false;
			for( i = 0; i < sizeSampleGraph; i++ )
				if( sampleGraphs_C[i] == stringG )
				{
					graphWeights[i] += allWeights[counterallG];
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
			sumWeights += 1 / sumRates;
		}
// Finished saving result -----------------------------------------------------|	

		selectedEdgeij    = selectedEdgej * dim + selectedEdgei;
		G[selectedEdgeij] = 1 - G[selectedEdgeij];
		G[selectedEdgei * dim + selectedEdgej] = G[selectedEdgeij];

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim, &threshold_C, &sigma_start[0], &invC[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sampleGraphs_C[i].copy(sampleGraphs[i], qp, 0);
		sampleGraphs[i][qp] = '\0';
	}
	
	*sizeSampleG = sizeSampleGraph;

	for( i = 0; i < pxp; i++ ) 
		Khat[i] /= sumWeights;
}
    
void bdmcmcCopulaPhat( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double Z[], int R[], int *n, int *gcgm,
			 double Khat[], double phat[],
			 int *b, int *bstar, double D[], double Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin;
	
	int selectedEdgei, selectedEdgej, selectedEdgeij, l;

	double Dsjj, Dsij, sumDiag, K022, a11, b1 = *b, sigmaj11, threshold_C = *threshold;
	int row, col, rowCol, i, j, k, ij, jj, nustar, one = 1, two = 2, dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	double sumWeights = 0.0, rate, maxRates, sumRates; 
	
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
// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> invC( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
// ----------------------------

	for( int g = 0; g < iteration; g++ )
	{
		if( ( g + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                  \n", g + 1 );

		getDs( K, Z, R, D, Ds, gcgm, n, &dim );

		getTs( Ds, Ts, &dim );
		
		l = 0;
		sumRates = 0.0;
		maxRates = -1.7e307; 

		// computing birth and death rates
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
				subRowMins( K, &Kj12[0], &j, &dim );  // K12 = K[j, -j]  
				Kj12[ i ] = 0.0;                      // K12[1,i] = 0

				// K12 %*% K22_inv
				//multiplyMatrix( &Kj12[0], &Kj22_inv[0], &Kj12xK22_inv[0], &one, &p1, &p1 );															
				F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &Kj12[0], &one, &Kj22_inv[0], &p1, &beta, &Kj12xK22_inv[0], &one );
				
				// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
				F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			

	// For (i,j) = 1 ---------------------------------------------------------------|
				subRowsMins( K, &K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
				
				subMatrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &i, &j, &dim );

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

				a11     = K[i * dim + i] - K121[0];	
				sumDiag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

				//	nustar = b + sum( Gf[,i] * Gf[,j] )
				nustar = b1;
				for( k = 0; k < dim; k++ )
					nustar += G[i * dim + k] * G[j * dim + k];

				rate = sqrt( 2.0 * Dsjj / a11 ) * exp( lgamma( (nustar + 1) / 2 ) - lgamma( nustar / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sumDiag ) / 2 );

				if( G[ij] == 0 ) rate = 1.0 / rate;		

				sumRates += rate;        // sumUpperMatrix( &rates[0], &allWeights[g], p );
				
				if( rate > maxRates )    // selectEdge( &rates[0], selectedEdge, p );
				{
					maxRates      = rate; 
					selectedEdgei = i;
					selectedEdgej = j;
				}
			}
		}
		
// saving result --------------------------------------------------------------|	
		if( g >= burn_in )
		{
			for( i = 0; i < pxp ; i++ )
			{
				Khat[i] += K[i] / sumRates;
				phat[i] += double( G[i] ) / sumRates;
			}
			
			sumWeights += 1 / sumRates;
		}	
// Finished saving result -----------------------------------------------------|	

		selectedEdgeij    = selectedEdgej * dim + selectedEdgei;
		G[selectedEdgeij] = 1 - G[selectedEdgeij];
		G[selectedEdgei * dim + selectedEdgej] = G[selectedEdgeij];

		rgwish_sigma( G, Ts, K, &sigma[0], bstar, &dim, &threshold_C, &sigma_start[0], &invC[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}

	for( i = 0; i < pxp; i++ ) 
	{
		phat[i] /= sumWeights;
		Khat[i] /= sumWeights;
	}
}
    
} // exturn "C"


