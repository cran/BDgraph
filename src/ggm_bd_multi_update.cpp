#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <limits>        // for std::numeric_limits<long double>::max()
#include <R.h>
#include <Rmath.h>
#include "matrix.h"
#include "rgwish.h"

using namespace std;

extern "C" {

/*
 * Multiple birth-death MCMC for Gaussian Graphical models  
 * with exact value of normalizing constant for D = I_p 
 * it is for Bayesian model averaging
*/
void bdmcmcExactp_links_multi_update( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double K_hat[], double p_links[],
			 int *b, int *b_star, double Ds[], double *threshold, int *multi_update )
{
	int iteration = *iter, burn_in = *burnin;
	int multi_update_C = *multi_update;

	int selected_edge_i, selected_edge_j, selected_edge_ij;
	int row, col, rowCol, i, j, k, ij, jj, Dsjj, Dsij, counter, nu_star, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	double sum_weights = 0.0, sum_diag, K022, a11, b1 = *b, sigmaj11, threshold_C = *threshold;
	long double rate, sum_rates;

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	// For finding the index of rates 
	vector<long double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 
	vector<double> K12( p2x2 );             // K[e, -e]
	vector<double> sigma11( 4 );            // sigma[e, e]
	vector<double> sigma12( p2x2 );         // sigma[e, -e]
	vector<double> sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ----------------------------

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for multiple birth-death MCMC sampling algorithm -------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
		counter = 0;
		// STEP 1: calculating birth and death rates --------------------------|			
		for( j = 1; j < dim; j++ )
		{
			jj   = j * dim + j;
			Dsjj = Ds[jj];
			
			sigmaj11 = sigma[jj];        // sigma[j, j]  
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			for( row = 0; row < p1; row++ )
				for( col = 0; col < p1; col++ )
				{
					rowCol           = col * p1 + row;
					Kj22_inv[rowCol] = sigmaj22[rowCol] - sigmaj12[row] * sigmaj12[col] / sigmaj11;
				}
			
			for( i = 0; i < j; i++ )
			{
				ij   = j * dim + i;
				Dsij = Ds[ij];

				// For (i,j) = 0 ----------------------------------------------|	
				sub_row_mins( K, &Kj12[0], &j, &dim );    // K12 = K[j, -j]  
				Kj12[ i ] = 0.0;                          // K12[1,i] = 0

				// K12 %*% K22_inv
				F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &Kj12[0], &one, &Kj22_inv[0], &p1, &beta, &Kj12xK22_inv[0], &one );
				
				// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12)
				F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			

				// For (i,j) = 1 ----------------------------------------------|
				sub_rows_mins( K, &K12[0], &i, &j, &dim );   // K12 = K[e, -e]  
				
				sub_matrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &i, &j, &dim );

				// solve( sigma[e, e] )
				inverse_2x2( &sigma11[0], &sigma11_inv[0] );

				// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
				F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma21xsigma11_inv %*% sigma12
				F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &sigma21xsigma11_inv[0], &p2, &sigma12[0], &two, &beta, &sigma2112[0], &p2 );

				// solve( K[-e, -e] ) = sigma22 - sigma2112
				for( k = 0; k < p2xp2 ; k++ ) 
					K22_inv[k] = sigma22[k] - sigma2112[k];	
				
				// K12 %*% K22_inv
				F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &K12[0], &two, &K22_inv[0], &p2, &beta, &K12xK22_inv[0], &two );
				
				// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 														
				F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
				// Finished (i,j) = 1------------------------------------------|

				a11      = K[i * dim + i] - K121[0];	
				sum_diag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

				// nu_star = b + sum( Gf[,i] * Gf[,j] )
				nu_star = b1;
				for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k] * G[j * dim + k];

				rate = ( G[ij] ) 
					? sqrt( Dsjj / a11 ) * exp( lgamma( ( nu_star + 1 ) / 2 ) - lgamma( nu_star / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sum_diag ) / 2 )
					: sqrt( a11 / Dsjj ) * exp( lgamma( nu_star / 2 ) - lgamma( ( nu_star + 1 ) / 2 ) + ( Dsij * Dsij * a11 / Dsjj  + sum_diag ) / 2 );
				
				rates[counter++] = ( R_FINITE( rate ) ) ? rate : max_numeric_limits_ld;
			}
		}	

		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}

		// STEP 2: Sampling from G-Wishart for new graph ----------------------|	
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < pxp ; i++ )
			{
				K_hat_Cpp[i]   += K[i] / sum_rates;
				if( G[i] ) p_links_Cpp[i] += 1.0 / sum_rates;
			}
			
			sum_weights += 1.0 / sum_rates;
		} // End of saving result ---------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------|
	PutRNGstate();

	for( i = 0; i < pxp; i++ )
	{	
		p_links[i] = p_links_Cpp[i] / sum_weights;
		K_hat[i]   = K_hat_Cpp[i] / sum_weights;
	}
}
    
/*
 * Multiple birth-death MCMC for Gaussian Graphical models  
 * with exact value of normalizing constant for D = I_p 
 * it is for maximum a posterior probability estimation (MAP)
*/
void bdmcmcExact_multi_update( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g,
			 int *b, int *b_star, double Ds[], double *threshold, int *multi_update )
{
	int multi_update_C = *multi_update;
	int iteration = *iter, burn_in = *burnin;
	int counterallG = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	bool this_one;

	int selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int row, col, rowCol, i, j, k, ij, jj, Dsjj, Dsij, counter, nu_star, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;

	double sum_weights = 0.0, sum_diag, K022, a11, b1 = *b, sigmaj11, threshold_C = *threshold;
	long double rate, sum_rates;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	

	// Counting size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<long double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> Kj22_inv( p1xp1 ); 
	vector<double> Kj12xK22_inv( p1 ); 
	vector<double> K12( p2x2 );             // K[e, -e]
	vector<double> sigma11( 4 );            // sigma[e, e]
	vector<double> sigma12( p2x2 );         // sigma[e, -e]
	vector<double> sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2xp2 ); 
	vector<double> K22_inv( p2xp2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// ---- for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// ----------------------------

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for multiple birth-death MCMC sampling algorithm -------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
		counter = 0;
		// STEP 1: calculating birth and death rates
		for( j = 1; j < dim; j++ )
		{
			jj   = j * dim + j;
			Dsjj = Ds[jj];
			
			sigmaj11 = sigma[jj];        // sigma[j, j]  
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

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

				// For (i,j) = 0 ----------------------------------------------|	
				sub_row_mins( K, &Kj12[0], &j, &dim );   // K12 = K[j, -j]  
				Kj12[ i ] = 0.0;                         // K12[1,i] = 0

				// Kj12 %*% Kj22_inv
				F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &Kj12[0], &one, &Kj22_inv[0], &p1, &beta, &Kj12xK22_inv[0], &one );
				
				// K022  <- Kj12 %*% solve( K0[-j, -j] ) %*% t(Kj12)
				F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1, &alpha, &Kj12xK22_inv[0], &one, &Kj12[0], &one, &beta, &K022, &one );			

				// For (i,j) = 1 ----------------------------------------------|
				sub_rows_mins( K, &K12[0], &i, &j, &dim );   // K12 = K[e, -e]  
				
				sub_matrices( &sigma[0], &sigma11[0], &sigma12[0], &sigma22[0], &i, &j, &dim );

				// solve( sigma[e, e] )
				inverse_2x2( &sigma11[0], &sigma11_inv[0] );

				// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
				F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &sigma12[0], &two, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma21xsigma11_inv %*% sigma12
				F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &sigma21xsigma11_inv[0], &p2, &sigma12[0], &two, &beta, &sigma2112[0], &p2 );

				// solve( K[-e, -e] ) = sigma22 - sigma2112
				for( k = 0; k < p2xp2 ; k++ ) 
					K22_inv[k] = sigma22[k] - sigma2112[k];	
				
				// K12 %*% K22_inv
				F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &K12[0], &two, &K22_inv[0], &p2, &beta, &K12xK22_inv[0], &two );
				
				// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 														
				F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K12[0], &two, &beta, &K121[0], &two );		
				// Finished (i,j) = 1------------------------------------------|

				a11      = K[i * dim + i] - K121[0];	
				sum_diag = - Dsij * K121[1] - Dsij * K121[2] + Dsjj * ( K022 - K121[3] );

				// nu_star = b + sum( Gf[,i] * Gf[,j] )
				nu_star = b1;
				for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k] * G[j * dim + k];   

				rate = ( G[ij] ) 
					? sqrt( Dsjj / a11 ) * exp( lgamma( ( nu_star + 1 ) / 2 ) - lgamma( nu_star / 2 ) - ( Dsij * Dsij * a11 / Dsjj  + sum_diag ) / 2 )
					: sqrt( a11 / Dsjj ) * exp( lgamma( nu_star / 2 ) - lgamma( ( nu_star + 1 ) / 2 ) + ( Dsij * Dsij * a11 / Dsjj  + sum_diag ) / 2 );
				
				rates[counter] = ( R_FINITE( rate ) ) ? rate : max_numeric_limits_ld;
				
				char_g[counter] = G[ij] + '0'; 
				counter++; 
			}
		}	
			
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}

		// STEP 2: Sampling from G-Wishart for new graph ----------------------|
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &threshold_C, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			for( i = 0; i < pxp ; i++ ) K_hat[i] += K[i] / sum_rates;	

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[counterallG] = 1.0 / sum_rates;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[counterallG];
					all_graphs[counterallG] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[counterallG];
				all_graphs[counterallG]          = size_sample_graph; 
				size_sample_graph++;				
			}
			
			counterallG++; 
			sum_weights += 1.0 / sum_rates;
		} // End of saving result ---------------------------------------------|	
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	
	for( i = 0; i < pxp; i++ ) K_hat[i] /= sum_weights;		
}
              
} // End of exturn "C"
