#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <math.h>        // isinf, sqrt
#include <limits>        // for numeric_limits<long double>::max()
#include <R.h>
#include <Rmath.h>
#include "matrix.h"
#include "rgcwish.h"
using namespace std;

extern "C" {
/*
 * birth-death MCMC for time series with Gaussian Graphical models  
 * for D = I_p 
 * it is for Bayesian model averaging
*/
void bdmcmc_for_multi_dim( int *iter, int *burnin, int G[], double Ls[], double r_K[], double i_K[], int *p, 
			   int *freq, double r_sigma[], double i_sigma[], double r_K_hat[], double i_K_hat[], double p_links[],
			   int *b, int *b_star, double r_Ds[], double i_Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin, b1 = 0, T = *freq;

	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij;
	int row, col, rowCol, i, j, t, k, ij, jj, counter, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;
	int pxpxT = pxp * T, txpxp, n = dim, max_b_star = 0;
	bool useful;

	double r_Dsjj, i_Dsjj, r_Dsij, i_Dsij, sum_weights = 0.0, r_sum_diag, r_K022, i_K022, r_a11, i_a11, r_sigmaj11, i_sigmaj11, threshold_C = *threshold;
	long double mod_Dsjj, mod_a11, coef, r_temp, nu_star, rate, sum_rates, G_prior, common_factor = 1.0;
	
	for( i = 0; i < T; i++ )
	{
		b1 = b1 + b[i];
		if( b_star[i] > max_b_star )
			max_b_star = b_star[i];
	}
	
	n  = n + max_b_star;
	b1 = b1 / T;
	// for a matrix A, A[i,j] -> A[(j-1)*p+(i-1)]
	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> r_K_hat_Cpp( pxpxT, 0.0 ); 
	vector<double> i_K_hat_Cpp( pxpxT, 0.0 ); 			

	int qp = dim * ( dim - 1 ) / 2;

	double alpha = 1.0, beta = 0.0, dmone = -1.0;
	char transT = 'T', transN = 'N';																	

	// Counting size of nodes
	int ip;
	vector<int> size_node( dim, 0 );      // degrees of vertex
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<long double> log_rates( qp );
	vector<long double> rates( qp );      // the rates of all edges
	vector<int> index_rates_row( qp );    // 0,0,1,0,1,2,...
	vector<int> index_rates_col( qp );    // 1,2,2,3,3,3,...
	vector<bool> is_infinite(qp, 0);
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> r_K121( 4 ); 
	vector<double> i_K121( 4 ); 
	vector<double> r_Kj12( p1 );              // K[j, -j]
	vector<double> i_Kj12( p1 );              
	vector<double> r_sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> i_sigmaj12( p1 );           
	vector<double> r_sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> i_sigmaj22( p1xp1 );       
	vector<double> r_Kj22_inv( p1xp1 ); 
	vector<double> i_Kj22_inv( p1xp1 ); 
	vector<double> i12xi22_j( p1 ); 
	vector<double> r12xi22_j( p1 );
	vector<double> i12xi22( p2x2 ); 
	vector<double> r12xi22( p2x2 );	 
	vector<double> r21xr11( p2x2 ); 
	vector<double> i21xr11( p2x2 ); 
	vector<double> r_K12( p2x2 );             // K[e, -e]
	vector<double> i_K12( p2x2 );             
	vector<double> r_sigma11( 4 );            // sigma[e, e]
	vector<double> i_sigma11( 4 );            
	vector<double> r_sigma12( p2x2 );         // sigma[e, -e]
	vector<double> i_sigma12( p2x2 );         
	vector<double> r_sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> i_sigma22( p2xp2 );        
	vector<double> r_sigma11_inv( 4 ); 
	vector<double> i_sigma11_inv( 4 ); 
	vector<double> r_sigma2112( p2xp2 ); 
	vector<double> i_sigma2112( p2xp2 ); 
	vector<double> r_K22_inv( p2xp2 ); 
	vector<double> i_K22_inv( p2xp2 ); 
	vector<double> r_K12xK22_inv( p2x2 );   
	vector<double> i_K12xK22_inv( p2x2 );   
	// ---- for rgcwish_sigma 
	Rcomplex *csigma = new Rcomplex[pxp];
	Rcomplex *Ind = new Rcomplex[pxp];
	Rcomplex *K = new Rcomplex[pxp];
	vector<double> X(dim*n);
	vector<double> Y(dim*n);
	vector<double> joint(2*dim*n);
	vector<double> IR(pxp);
	vector<double> Inv_R(pxp);
	vector<double> r_sigma_start( pxp ); 
	vector<double> i_sigma_start( pxp ); 
	vector<double> r_beta_star( dim ); 
	vector<double> i_beta_star( dim ); 
	vector<double> r_sigma_i( dim ); 
	vector<double> i_sigma_i( dim ); 
	vector<double> r_sigma_start_N_i( dim );     // For dynamic memory used
	vector<double> r_sigma_start_N_i_2( dim );   // For dynamic memory used
	vector<double> i_sigma_start_N_i( dim );     // For dynamic memory used
	vector<double> r_sigma_N_i( pxp );           // For dynamic memory used
	vector<double> r_sigma_N_i_2( pxp );         // For dynamic memory used
	vector<double> i_sigma_N_i( pxp );           // For dynamic memory used
	vector<int> N_i( dim );                      // For dynamic memory used
	// ----------------------------

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 

		// STEP 1: calculating birth and death rates --------------------------|			
		counter = 0;
		useful = FALSE;
		for( j = 1; j < dim; j++ )  // The terms related to G in rate
			for( i = 0; i < j; i++ )
			{
				nu_star = b1 + 0.0;
				for( k = 0; k < dim; k++ ) // nu_star = b + sum( Gf[,i] * Gf[,j] )
					nu_star += G[i * dim + k] * G[j * dim + k];   

				G_prior = lgamma( 0.5 * ( nu_star + 1 ) ) - lgamma( 0.5 * nu_star );
				log_rates[counter++] = ( G[j * dim + i] ) ? G_prior : - G_prior;
			}

		for( j = 0; j < qp; j++ ) 
			is_infinite[j] = 0;   // Record if the rate of edge j is infinite

		for( t = 0; t < T; t++ ) // The loop for the frequencies
		{
			txpxp   = t*pxp;
			counter = 0;		
			
			for( j = 1; j < dim; j++ )
			{		
				jj     = j * dim + j + txpxp;
				r_Dsjj = r_Ds[jj];
				i_Dsjj = i_Ds[jj];
				
				r_sigmaj11 = r_sigma[jj];        // sigma[j, j]  
				i_sigmaj11 = i_sigma[jj]; 
				sub_matrices1( &r_sigma[txpxp], &r_sigmaj12[0], &r_sigmaj22[0], &j, &dim );
				Hsub_matrices1( &i_sigma[txpxp], &i_sigmaj12[0], &i_sigmaj22[0], &j, &dim );

				// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
				for( row = 0; row < p1; row++ )       
					for( col = 0; col < p1; col++ )
					{
						rowCol             = col * p1 + row;
						r_Kj22_inv[rowCol] = r_sigmaj22[rowCol] - (r_sigmaj11*(r_sigmaj12[row]*r_sigmaj12[col] + i_sigmaj12[row]*i_sigmaj12[col]) + i_sigmaj11*(r_sigmaj12[row]*i_sigmaj12[col] - i_sigmaj12[row]*r_sigmaj12[col]))/(r_sigmaj11*r_sigmaj11 + i_sigmaj11*i_sigmaj11);
						i_Kj22_inv[rowCol] = i_sigmaj22[rowCol] - (r_sigmaj11*(r_sigmaj12[row]*i_sigmaj12[col] - i_sigmaj12[row]*r_sigmaj12[col]) - i_sigmaj11*(r_sigmaj12[row]*r_sigmaj12[col] + i_sigmaj12[row]*i_sigmaj12[col]))/(r_sigmaj11*r_sigmaj11 + i_sigmaj11*i_sigmaj11);
					}		

				for( i = 0; i < j; i++ )
				{
					ij = j * dim + i;
					if( is_infinite[counter] ) // If the rate of edge counter has already been infinite, then we consider the
					{						   // final rate of it as infinite.
						counter++;
						continue;
					}
					ij    += txpxp;
					r_Dsij = r_Ds[ij];
					i_Dsij = i_Ds[ij];
					
					// For (i,j) = 0 ----------------------------------------------|	
					sub_row_mins( &r_K[txpxp], &r_Kj12[0], &j, &dim );    // K12 = K[j, -j]  
					Hsub_row_mins( &i_K[txpxp], &i_Kj12[0], &j, &dim );   // K12 = K[j, -j]  
					r_Kj12[ i ] = 0.0;                         // K12[i] = 0
					i_Kj12[ i ] = 0.0;                         // K12[i] = 0

					// Kj12xK22_inv = Kj12 %*% Kj22_inv
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &i_Kj12[0], &one, &i_Kj22_inv[0], &p1, &beta, &i12xi22_j[0], &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &r_Kj12[0], &one, &r_Kj22_inv[0], &p1, &dmone, &i12xi22_j[0], &one );				
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &r_Kj12[0], &one, &i_Kj22_inv[0], &p1, &beta, &r12xi22_j[0], &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &i_Kj12[0], &one, &r_Kj22_inv[0], &p1, &alpha, &r12xi22_j[0], &one );				
					// K022  <- Kj12 %*% solve( K0[-j, -j] ) %*% t(Kj12) = c
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &r12xi22_j[0], &one, &i_Kj12[0], &p1, &beta, &r_K022, &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &i12xi22_j[0], &one, &r_Kj12[0], &p1, &alpha, &r_K022, &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &i12xi22_j[0], &one, &i_Kj12[0], &p1, &beta, &i_K022, &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &r12xi22_j[0], &one, &r_Kj12[0], &p1, &dmone, &i_K022, &one );
					// For (i,j) = 1 ----------------------------------------------|
					sub_rows_mins( &r_K[txpxp], &r_K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
					Hsub_rows_mins( &i_K[txpxp], &i_K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
					
					sub_matrices( &r_sigma[txpxp], &r_sigma11[0], &r_sigma12[0], &r_sigma22[0], &i, &j, &dim ); //r_sigma[e,e], r_sigma[e,-e], r_sigma[-e,-e]
					Hsub_matrices( &i_sigma[txpxp], &i_sigma11[0], &i_sigma12[0], &i_sigma22[0], &i, &j, &dim ); //i_sigma[e,e], i_sigma[e,-e], i_sigma[-e,-e]

					// solve( sigma[e, e] )
					cinverse_2x2( &r_sigma11[0], &i_sigma11[0], &r_sigma11_inv[0], &i_sigma11_inv[0] );

					// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &r_sigma12[0], &two, &r_sigma11_inv[0], &two, &beta, &r21xr11[0], &p2 );
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &i_sigma12[0], &two, &i_sigma11_inv[0], &two, &alpha, &r21xr11[0], &p2 );
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &i_sigma12[0], &two, &r_sigma11_inv[0], &two, &beta, &i21xr11[0], &p2 );
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &r_sigma12[0], &two, &i_sigma11_inv[0], &two, &dmone, &i21xr11[0], &p2 );				
					// sigma2112 = sigma21xsigma11_inv %*% sigma12 = sigma[-e,e] %*% solve(sigma[e,e]) %*% sigma[e,-e]
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &i21xr11[0], &p2, &i_sigma12[0], &two, &beta, &r_sigma2112[0], &p2 );
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &r21xr11[0], &p2, &r_sigma12[0], &two, &dmone, &r_sigma2112[0], &p2 );
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &r21xr11[0], &p2, &i_sigma12[0], &two, &beta, &i_sigma2112[0], &p2 );
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &i21xr11[0], &p2, &r_sigma12[0], &two, &alpha, &i_sigma2112[0], &p2 );


					// solve( K[-e, -e] ) = sigma22 - sigma2112
					for( k = 0; k < p2xp2 ; k++ ) 
					{
						r_K22_inv[k] = r_sigma22[k] - r_sigma2112[k];
						i_K22_inv[k] = i_sigma22[k] - i_sigma2112[k];
					}

					// K12 %*% K22_inv
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &i_K12[0], &two, &i_K22_inv[0], &p2, &beta, &i12xi22[0], &two );
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &r_K12[0], &two, &r_K22_inv[0], &p2, &dmone, &i12xi22[0], &two );
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &r_K12[0], &two, &i_K22_inv[0], &p2, &beta, &r12xi22[0], &two );
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &i_K12[0], &two, &r_K22_inv[0], &p2, &alpha, &r12xi22[0], &two );
					// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 		
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &r12xi22[0], &two, &i_K12[0], &two, &beta, &r_K121[0], &two );
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &i12xi22[0], &two, &r_K12[0], &two, &alpha, &r_K121[0], &two );
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &i12xi22[0], &two, &i_K12[0], &two, &beta, &i_K121[0], &two );
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &r12xi22[0], &two, &r_K12[0], &two, &dmone, &i_K121[0], &two );
																			
					// Finished (i,j) = 1------------------------------------------|

					r_a11      = r_K[i * dim + i + txpxp] - r_K121[0]; //k_ii - k_ii^1
					i_a11      = i_K[i * dim + i + txpxp] - i_K121[0]; //k_ii - k_ii^1
					r_sum_diag = r_Dsjj*(r_K022 - r_K121[3]) - (r_Dsij*r_K121[1] - i_Dsij*i_K121[1]) - (r_Dsij*r_K121[2] + i_Dsij*i_K121[2]); //tr(D*(K0-K1))

					mod_Dsjj = sqrt( r_Dsjj * r_Dsjj + i_Dsjj * i_Dsjj );
					mod_a11  = sqrt( r_a11 * r_a11 + i_a11 * i_a11 );
					coef     = ( r_Dsij * r_Dsij + i_Dsij * i_Dsij ) / ( r_Dsjj * r_Dsjj + i_Dsjj * i_Dsjj );
					r_temp   = coef * ( r_a11 * r_Dsjj + i_a11 * i_Dsjj ) + r_sum_diag;
					rate     = ( G[ij - txpxp] ) ? mod_Dsjj / mod_a11 * exp( -r_temp ) : mod_a11 / mod_Dsjj * exp( r_temp );
					
					if( R_FINITE( rate ) )   // Computer the rate in log space
						log_rates[counter++] += log( static_cast<double>( rate ) );
					else
					{
						rates[counter] = max_numeric_limits_ld;
						is_infinite[counter++] = true;
					}
				}
			}
		}// end of frequency loop	
		
		// Selecting an edge based on birth and death rates
		for( t = 0; t < qp; t++ ) // Exponentialize the log rate
		{
			rate = exp( log_rates[t] );
			if( !is_infinite[t] && R_FINITE( rate ) )
				rates[t] = rate;
			else
				rates[t] = max_numeric_limits_ld;
		}

		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];
		
		// If sum_rates = 0, we consider the graph in this state is the true graph
		if( sum_rates == 0 )
		{
			for( t = 0; t < pxpxT; t++ )
			{
				r_K_hat[t] = r_K[t];
				i_K_hat[t] = i_K[t];
			}
			
			for( i = 0; i < pxp; i++ )
				if( G[i] ) p_links[i] = 1.0;
			delete[] csigma, Ind, K;
			return;
		}		
		
		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			// if( sum_rates < constant )
			// sum_rates = constant;
			if( !useful && sum_rates < 1 ) // Set the first huge sum_rates as the common_factor and multiply it when calculating
			{							   // K_hat_Cpp, p_links_Cpp and sum_weights to avoid overflow
				common_factor = sum_rates;
				useful        = TRUE;
			}
			
			for( t = 0; t < pxpxT; t++ )
			{
				r_K_hat_Cpp[t] += r_K[t] * ( common_factor / sum_rates );
				i_K_hat_Cpp[t] += i_K[t] * ( common_factor / sum_rates );
			}
			
			for ( i = 0; i < pxp ; i++ )
				if( G[i] ) p_links_Cpp[i] += 1.0 * ( common_factor / sum_rates );

			sum_weights += 1.0 * ( common_factor / sum_rates );
		} // End of saving result ---------------------------------------------|	

		// Updating G (graph) based on selected edge
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

		// STEP 2: Sampling from G-Wishart for new graph ----------------------|
		for( t = 0; t < T; t++ )
		{
			txpxp = t*pxp;
			rgcwish_sigma(G, &size_node[0], &Ls[4*txpxp], K, &r_sigma[txpxp], 
			&i_sigma[txpxp], csigma, Ind, &b_star[t], &dim, &threshold_C, 
			&r_sigma_start[0], &i_sigma_start[0], &X[0], &Y[0],
			&r_beta_star[0], &i_beta_star[0], &joint[0], &r_sigma_i[0],
			&i_sigma_i[0], r_sigma_start_N_i, r_sigma_start_N_i_2,
			i_sigma_start_N_i, r_sigma_N_i, r_sigma_N_i_2, i_sigma_N_i, N_i, IR, Inv_R);
			
			for( k = 0; k < pxp; k++ )
			{
				r_K[k + txpxp] = K[k].r;
				i_K[k + txpxp] = K[k].i;
			}		
		}
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < pxpxT; i++ )
	{
		r_K_hat[i] = r_K_hat_Cpp[i] / sum_weights;
		i_K_hat[i] = i_K_hat_Cpp[i] / sum_weights;
	}
	
	for( i = 0; i < pxp; i++ )
		p_links[i] = p_links_Cpp[i] / sum_weights;
		
	delete[] csigma, Ind, K;
}
 
/*
 * birth-death MCMC for time series with Gaussian Graphical models  
 * for case D = I_p 
 * it is for maximum a posterior probability estimation (MAP)
*/
void bdmcmc_map_for_multi_dim( int *iter, int *burnin, int G[], double Ls[], double r_K[], double i_K[], int *p, 
			   int *freq, double r_sigma[], double i_sigma[], int all_graphs[], double all_weights[], double r_K_hat[], 
			   double i_K_hat[], char *sample_graphs[], double graph_weights[], int *size_sample_g, int *exit,
			   int *b, int *b_star, double r_Ds[], double i_Ds[], double *threshold )
{
	int iteration = *iter, burn_in = *burnin, b1 = 0, T = *freq;
	int counterallG = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );

	bool this_one;

	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int row, col, rowCol, i, j, t, k, ij, jj, counter, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2;
	int pxpxT = pxp * T, txpxp, n = dim, max_b_star = 0;
	bool useful;

	double r_Dsjj, i_Dsjj, r_Dsij, i_Dsij, sum_weights = 0.0, r_sum_diag, r_K022, i_K022, r_a11, i_a11, r_sigmaj11, i_sigmaj11, threshold_C = *threshold;
	long double mod_Dsjj, mod_a11, coef, r_temp, nu_star, rate, sum_rates, G_prior, common_factor = 1.0;
	
	for ( i = 0; i < T; i++ )
	{
		b1 = b1 + b[i];
		if (b_star[i] > max_b_star)
			max_b_star = b_star[i];
	}
	n = n + max_b_star;
	b1 = b1 / T;			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	double alpha = 1.0, beta = 0.0, dmone = -1.0;
	char transT = 'T', transN = 'N';																	

	// Counting size of nodes
	int ip;
	vector<int> size_node( dim, 0 );      // degrees of vertex
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<long double> log_rates( qp );
	vector<long double> rates( qp );      // the rates of all edges
	vector<int> index_rates_row( qp );    // 0,0,1,0,1,2,...
	vector<int> index_rates_col( qp );    // 1,2,2,3,3,3,...
	vector<bool> is_infinite(qp, 0);
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	vector<double> r_K121( 4 ); 
	vector<double> i_K121( 4 ); 
	vector<double> r_Kj12( p1 );              // K[j, -j]
	vector<double> i_Kj12( p1 );              
	vector<double> r_sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> i_sigmaj12( p1 );           
	vector<double> r_sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> i_sigmaj22( p1xp1 );       
	vector<double> r_Kj22_inv( p1xp1 ); 
	vector<double> i_Kj22_inv( p1xp1 ); 
	vector<double> i12xi22_j( p1 ); 
	vector<double> r12xi22_j( p1 );
	vector<double> i12xi22( p2x2 ); 
	vector<double> r12xi22( p2x2 );	 
	vector<double> r21xr11( p2x2 ); 
	vector<double> i21xr11( p2x2 ); 
	vector<double> r_K12( p2x2 );             // K[e, -e]
	vector<double> i_K12( p2x2 );             
	vector<double> r_sigma11( 4 );            // sigma[e, e]
	vector<double> i_sigma11( 4 );            
	vector<double> r_sigma12( p2x2 );         // sigma[e, -e]
	vector<double> i_sigma12( p2x2 );         
	vector<double> r_sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> i_sigma22( p2xp2 );        
	vector<double> r_sigma11_inv( 4 ); 
	vector<double> i_sigma11_inv( 4 ); 
	vector<double> r_sigma2112( p2xp2 ); 
	vector<double> i_sigma2112( p2xp2 ); 
	vector<double> r_K22_inv( p2xp2 ); 
	vector<double> i_K22_inv( p2xp2 ); 
	vector<double> r_K12xK22_inv( p2x2 );   
	vector<double> i_K12xK22_inv( p2x2 );   
	// ---- for rgcwish_sigma 
	Rcomplex *csigma = new Rcomplex[pxp];
	Rcomplex *Ind = new Rcomplex[pxp];
	Rcomplex *K = new Rcomplex[pxp];
	vector<double> X(dim*n);
	vector<double> Y(dim*n);
	vector<double> joint(2*dim*n);
	vector<double> IR(pxp);
	vector<double> Inv_R(pxp);
	vector<double> r_sigma_start( pxp ); 
	vector<double> i_sigma_start( pxp ); 
	vector<double> r_beta_star( dim ); 
	vector<double> i_beta_star( dim ); 
	vector<double> r_sigma_i( dim ); 
	vector<double> i_sigma_i( dim ); 
	vector<double> r_sigma_start_N_i( dim );     // For dynamic memory used
	vector<double> r_sigma_start_N_i_2( dim );   // For dynamic memory used
	vector<double> i_sigma_start_N_i( dim );     // For dynamic memory used
	vector<double> r_sigma_N_i( pxp );           // For dynamic memory used
	vector<double> r_sigma_N_i_2( pxp );         // For dynamic memory used
	vector<double> i_sigma_N_i( pxp );           // For dynamic memory used
	vector<int> N_i( dim );                      // For dynamic memory used
	// ----------------------------

	long double max_numeric_limits_ld = std::numeric_limits<long double>::max() / 10000;

	GetRNGstate();
	// main loop for birth-death MCMC sampling algorithm for time series ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % 1000 == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 		

		counter = 0;
		useful = FALSE;
		for( j = 1; j < dim; j++ )  // The terms related to G in rate
			for( i = 0; i < j; i++ )
			{
				nu_star = b1 + 0.0;
				for( k = 0; k < dim; k++ ) // nu_star = b + sum( Gf[,i] * Gf[,j] )
					nu_star += G[i * dim + k] * G[j * dim + k];   

				G_prior = lgamma( 0.5 * ( nu_star + 1 ) ) - lgamma( 0.5 * nu_star );
				log_rates[counter++] = ( G[j * dim + i] ) ? G_prior : - G_prior;
			}

		for( j = 0; j < qp; j++ ) 
			is_infinite[j] = 0;   // Record if the rate of edge j is infinite

		for( t = 0; t < T; t++ ) // The loop for the frequencies
		{
			txpxp   = t*pxp;
			counter = 0;		
			// STEP 1: calculating birth and death rates --------------------------|
			for( j = 1; j < dim; j++ )
			{		
				jj     = j * dim + j + txpxp;
				r_Dsjj = r_Ds[jj];
				i_Dsjj = i_Ds[jj];
				
				r_sigmaj11 = r_sigma[jj];        // sigma[j, j]  
				i_sigmaj11 = i_sigma[jj]; 
				sub_matrices1( &r_sigma[txpxp], &r_sigmaj12[0], &r_sigmaj22[0], &j, &dim );
				Hsub_matrices1( &i_sigma[txpxp], &i_sigmaj12[0], &i_sigmaj22[0], &j, &dim );

				// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
				for( row = 0; row < p1; row++ )       
					for( col = 0; col < p1; col++ )
					{
						rowCol             = col * p1 + row;
						r_Kj22_inv[rowCol] = r_sigmaj22[rowCol] - (r_sigmaj11*(r_sigmaj12[row]*r_sigmaj12[col] + i_sigmaj12[row]*i_sigmaj12[col]) + i_sigmaj11*(r_sigmaj12[row]*i_sigmaj12[col] - i_sigmaj12[row]*r_sigmaj12[col]))/(r_sigmaj11*r_sigmaj11 + i_sigmaj11*i_sigmaj11);
						i_Kj22_inv[rowCol] = i_sigmaj22[rowCol] - (r_sigmaj11*(r_sigmaj12[row]*i_sigmaj12[col] - i_sigmaj12[row]*r_sigmaj12[col]) - i_sigmaj11*(r_sigmaj12[row]*r_sigmaj12[col] + i_sigmaj12[row]*i_sigmaj12[col]))/(r_sigmaj11*r_sigmaj11 + i_sigmaj11*i_sigmaj11);
					}		

				for( i = 0; i < j; i++ )
				{
					ij = j * dim + i;
					if( is_infinite[counter] ) // If the rate of edge counter has already been infinite, then we consider the
					{						   // final rate of it as infinite.
						counter++;
						continue;
					}
					
					ij    += txpxp;
					r_Dsij = r_Ds[ij];
					i_Dsij = i_Ds[ij];
					
					// For (i,j) = 0 ----------------------------------------------|	
					sub_row_mins( &r_K[txpxp], &r_Kj12[0], &j, &dim );    // K12 = K[j, -j]  
					Hsub_row_mins( &i_K[txpxp], &i_Kj12[0], &j, &dim );   // K12 = K[j, -j]  
					r_Kj12[ i ] = 0.0;                         // K12[i] = 0
					i_Kj12[ i ] = 0.0;                         // K12[i] = 0

					// Kj12xK22_inv = Kj12 %*% Kj22_inv
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &i_Kj12[0], &one, &i_Kj22_inv[0], &p1, &beta, &i12xi22_j[0], &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &r_Kj12[0], &one, &r_Kj22_inv[0], &p1, &dmone, &i12xi22_j[0], &one );				
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &r_Kj12[0], &one, &i_Kj22_inv[0], &p1, &beta, &r12xi22_j[0], &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &i_Kj12[0], &one, &r_Kj22_inv[0], &p1, &alpha, &r12xi22_j[0], &one );				
					// K022  <- Kj12 %*% solve( K0[-j, -j] ) %*% t(Kj12) = c
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &r12xi22_j[0], &one, &i_Kj12[0], &p1, &beta, &r_K022, &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &i12xi22_j[0], &one, &r_Kj12[0], &p1, &alpha, &r_K022, &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &i12xi22_j[0], &one, &i_Kj12[0], &p1, &beta, &i_K022, &one );
					F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &r12xi22_j[0], &one, &r_Kj12[0], &p1, &dmone, &i_K022, &one );
					// For (i,j) = 1 ----------------------------------------------|
					sub_rows_mins( &r_K[txpxp], &r_K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
					Hsub_rows_mins( &i_K[txpxp], &i_K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
					
					sub_matrices( &r_sigma[txpxp], &r_sigma11[0], &r_sigma12[0], &r_sigma22[0], &i, &j, &dim ); //r_sigma[e,e], r_sigma[e,-e], r_sigma[-e,-e]
					Hsub_matrices( &i_sigma[txpxp], &i_sigma11[0], &i_sigma12[0], &i_sigma22[0], &i, &j, &dim ); //i_sigma[e,e], i_sigma[e,-e], i_sigma[-e,-e]

					// solve( sigma[e, e] )
					cinverse_2x2( &r_sigma11[0], &i_sigma11[0], &r_sigma11_inv[0], &i_sigma11_inv[0] );

					// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &r_sigma12[0], &two, &r_sigma11_inv[0], &two, &beta, &r21xr11[0], &p2 );
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &i_sigma12[0], &two, &i_sigma11_inv[0], &two, &alpha, &r21xr11[0], &p2 );
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &i_sigma12[0], &two, &r_sigma11_inv[0], &two, &beta, &i21xr11[0], &p2 );
					F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &r_sigma12[0], &two, &i_sigma11_inv[0], &two, &dmone, &i21xr11[0], &p2 );				
					// sigma2112 = sigma21xsigma11_inv %*% sigma12 = sigma[-e,e] %*% solve(sigma[e,e]) %*% sigma[e,-e]
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &i21xr11[0], &p2, &i_sigma12[0], &two, &beta, &r_sigma2112[0], &p2 );
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &r21xr11[0], &p2, &r_sigma12[0], &two, &dmone, &r_sigma2112[0], &p2 );
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &r21xr11[0], &p2, &i_sigma12[0], &two, &beta, &i_sigma2112[0], &p2 );
					F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &i21xr11[0], &p2, &r_sigma12[0], &two, &alpha, &i_sigma2112[0], &p2 );


					// solve( K[-e, -e] ) = sigma22 - sigma2112
					for( k = 0; k < p2xp2 ; k++ ) 
					{
						r_K22_inv[k] = r_sigma22[k] - r_sigma2112[k];
						i_K22_inv[k] = i_sigma22[k] - i_sigma2112[k];
					}

					// K12 %*% K22_inv
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &i_K12[0], &two, &i_K22_inv[0], &p2, &beta, &i12xi22[0], &two );
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &r_K12[0], &two, &r_K22_inv[0], &p2, &dmone, &i12xi22[0], &two );
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &r_K12[0], &two, &i_K22_inv[0], &p2, &beta, &r12xi22[0], &two );
					F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &i_K12[0], &two, &r_K22_inv[0], &p2, &alpha, &r12xi22[0], &two );
					// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 		
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &r12xi22[0], &two, &i_K12[0], &two, &beta, &r_K121[0], &two );
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &i12xi22[0], &two, &r_K12[0], &two, &alpha, &r_K121[0], &two );
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &i12xi22[0], &two, &i_K12[0], &two, &beta, &i_K121[0], &two );
					F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &r12xi22[0], &two, &r_K12[0], &two, &dmone, &i_K121[0], &two );
																			
					// Finished (i,j) = 1------------------------------------------|

					r_a11      = r_K[i * dim + i + txpxp] - r_K121[0]; //k_ii - k_ii^1
					r_sum_diag = r_Dsjj * ( r_K022 - r_K121[3] ) - ( r_Dsij * r_K121[1] - i_Dsij * i_K121[1] ) - ( r_Dsij * r_K121[2] + i_Dsij * i_K121[2] ); //tr(D*(K0-K1))

					mod_Dsjj = fabs( r_Dsjj );
					mod_a11  = fabs( r_a11 );
					coef     = ( r_Dsij * r_Dsij + i_Dsij * i_Dsij ) / ( r_Dsjj * r_Dsjj );
					r_temp   = coef * ( r_a11 * r_Dsjj ) + r_sum_diag;
					rate     = ( G[ij - txpxp] ) ? mod_Dsjj / mod_a11 * exp( -r_temp ) : mod_a11 / mod_Dsjj * exp( r_temp );
					
					if( R_FINITE( rate ) )   // Computer the rate in log space
						log_rates[counter++] += log( static_cast<double>( rate ) );
					else
					{
						rates[counter] = max_numeric_limits_ld;
						is_infinite[counter++] = true;
					} 
				}
			}
		}// end of frequency loop	
		
		// Selecting an edge based on birth and death rates
		counter = 0;
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
				char_g[counter++] = G[j * dim + i] + '0'; 
		for( t = 0; t < qp; t++ ) // Exponentialize the log rate
		{
			rate = exp( log_rates[t] );
			if( !is_infinite[t] && R_FINITE( rate ) )
				rates[t] = rate;
			else
				rates[t] = max_numeric_limits_ld;
		}

		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];	
		
		// saving result ------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			if( sum_rates == 0 )
			{
				*exit = i_mcmc + 1;
				for( t = 0; t < pxpxT; t++ )
				{
					r_K_hat[t] = r_K[t];
					i_K_hat[t] = i_K[t];
				}
				
				string_g = string( char_g.begin(), char_g.end() ); // The adjacent matrix of G
				all_weights[counterallG] = max_numeric_limits_ld;

				this_one = false;
				for ( i = 0; i < size_sample_graph ; i++ ) // Check whether G appeared before
					if( sample_graphs_C[i] == string_g )
					{
						graph_weights[i] = all_weights[counterallG];
						all_graphs[counterallG] = i;
						this_one = true;
						break;
					} 			
					
				if( !this_one || size_sample_graph == 0 )  // If not, record a new graph
				{
					sample_graphs_C[size_sample_graph] = string_g;
					graph_weights[size_sample_graph]   = all_weights[counterallG];
					all_graphs[counterallG]            = size_sample_graph; 
					size_sample_graph++;	
					*size_sample_g = size_sample_graph;			
				}

				for( i = 0; i < ( i_mcmc - burn_in ); i++ )
				{
					sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
					sample_graphs[i][qp] = '\0';
				}	
				delete[] csigma, Ind, K;					
				return;				
			}
			else
			{
				if( !useful && sum_rates < 1 ) // Set the first huge sum_rates as the common_factor and multiply it when calculating
				{							   // K_hat_Cpp, p_links_Cpp and sum_weights to avoid overflow
					common_factor = sum_rates;
					useful        = TRUE;
				}
				
				for( t = 0; t < pxpxT; t++ )
				{
					r_K_hat[t] += r_K[t] * ( common_factor / sum_rates );
					i_K_hat[t] += i_K[t] * ( common_factor / sum_rates );
				}

				string_g = string( char_g.begin(), char_g.end() ); // The adjacent matrix of G
				all_weights[counterallG] = 1.0 / sum_rates;

				this_one = false;
				for ( i = 0; i < size_sample_graph ; i++ ) // Check whether G appeared before
					if( sample_graphs_C[i] == string_g )
					{
						graph_weights[i] += all_weights[counterallG];
						all_graphs[counterallG] = i;
						this_one = true;
						break;
					} 			
					
				if( !this_one || size_sample_graph == 0 )  // If not, record a new graph
				{
					sample_graphs_C[size_sample_graph] = string_g;
					graph_weights[size_sample_graph]   = all_weights[counterallG];
					all_graphs[counterallG]            = size_sample_graph; 
					size_sample_graph++;				
				}

				counterallG++;
				sum_weights += 1.0 * ( common_factor / sum_rates );
			}
		} // End of saving result ---------------------------------------------|	

		// If sum_rates = 0, we consider the graph in this state is the true graph
		if( sum_rates == 0 )
		{
			*exit = i_mcmc + 1;
			for( t = 0; t < pxpxT; t++ )
			{
				r_K_hat[t] = r_K[t];
				i_K_hat[t] = i_K[t];
			}
			delete[] csigma, Ind, K;
			return;
		}	

		// Updating G (graph) based on selected edge
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

		// STEP 2: Sampling from G-Wishart for new graph ----------------------|
		for( t = 0; t < T; t++ )
		{
			txpxp = t*pxp;
			rgcwish_sigma(G, &size_node[0], &Ls[4*txpxp], K, &r_sigma[txpxp], 
			&i_sigma[txpxp], csigma, Ind, &b_star[t], &dim, &threshold_C, 
			&r_sigma_start[0], &i_sigma_start[0], &X[0], &Y[0],
			&r_beta_star[0], &i_beta_star[0], &joint[0], &r_sigma_i[0],
			&i_sigma_i[0], r_sigma_start_N_i, r_sigma_start_N_i_2,
			i_sigma_start_N_i, r_sigma_N_i, r_sigma_N_i_2, i_sigma_N_i, N_i, IR, Inv_R);
			for (k = 0; k < pxp; k++)
			{
				r_K[k + txpxp] = K[k].r;
				i_K[k + txpxp] = K[k].i;
			}		
		}
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < size_sample_graph; i++ )
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	
	for( i = 0; i < pxpxT; i++ )
	{
		r_K_hat[i] /= sum_weights;
		i_K_hat[i] /= sum_weights;
	}
	delete[] csigma, Ind, K;
}

} // End of exturn "C"
