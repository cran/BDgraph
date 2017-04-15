// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2016 Mohammadi A. and Wit C. E.
//
//     This file is part of BDgraph package.
//
//     BDgraph is free software: you can redistribute it and/or modify it under 
//     the terms of the GNU General Public License as published by the Free 
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.
//
//     Maintainer:
//     Abdolreza Mohammadi: a.mohammadi@rug.nl or a.mohammadi@uvt.nl
// ----------------------------------------------------------------------------|
#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <math.h>        // isinf, sqrt
#include <R.h>
#include <Rmath.h>
#include <algorithm>     // for transform function
#include <functional>    // for transform function

#include "matrix.h"
#include "rgwish.h"

using namespace std;

extern "C" {
// ----------------------------------------------------------------------------|
// birth-death MCMC for Gaussian Graphical models  
// for case D = I_p 
// it is for Bayesian model averaging (MA)
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_ma( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double K_hat[], double p_links[],
			 int *b, int *b_star, double Ds[], int *print )
{
	int print_c = *print;
	int iteration = *iter, burn_in = *burnin, b1 = *b;

	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij;
	int i, j, k, ij, jj, counter, nu_star, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2, dim1 = dim + 1;

	double Dsjj, Dsij, sum_weights = 0.0, sum_diag, K022, a11, sigmajj_inv, weight_C;
	double log_rate, sum_rates;
	
	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 
	
	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;

	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideL = 'L';																	

	// Counting size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );             // K[-e, e]
	vector<double> sigma21( p2x2 );         // sigma[-e, e]
	vector<double> sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );         // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
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

	GetRNGstate();
	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
		counter = 0;
		// STEP 1: calculating birth and death rates --------------------------|		
		for( j = 1; j < dim; j++ )
		{		
			jj   = j * dim1;
			Dsjj = Ds[jj];
			
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
			sigmajj_inv = - 1.0 / sigma[jj];
			F77_NAME(dsyr)( &sideL, &p1, &sigmajj_inv, &sigmaj12[0], &one, &sigmaj22[0], &p1 );

			for( i = 0; i < j; i++ )
			{
				ij = j * dim + i;
				
				// For (i,j) = 0 ----------------------------------------------|	
				sub_row_mins( K, &Kj12[0], &j, &dim );   // Kj12 = K[j, -j]  
				Kj12[ i ] = 0.0;                         // Kj12[1,i] = 0

				// Kj12xK22_inv = Kj12 %*% Kj22_inv here sigmaj22 instead of Kj22_inv
				F77_NAME(dsymv)( &sideL, &p1, &alpha, &sigmaj22[0], &p1, &Kj12[0], &one, &beta, &Kj12xK22_inv[0], &one );
				
				// K022 = Kj12xK22_inv %*% t(Kj12)
				K022 = F77_NAME(ddot)( &p1, &Kj12xK22_inv[0], &one, &Kj12[0], &one );			

				// For (i,j) = 1 ----------------------------------------------|
				sub_cols_mins( K, &K21[0], &i, &j, &dim );  // K21 = K[-e, e]  
				
				sub_matrices_inv( &sigma[0], &sigma11_inv[0], &sigma21[0], &sigma22[0], &i, &j, &dim );

				// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
				F77_NAME(dgemm)( &transN, &transN, &p2, &two, &two, &alpha, &sigma21[0], &p2, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
				F77_NAME(dgemm)( &transN, &transT, &p2, &p2, &two, &alpha1, &sigma21xsigma11_inv[0], &p2, &sigma21[0], &p2, &beta1, &sigma22[0], &p2 );

				// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
				F77_NAME(dgemm)( &transT, &transN, &two, &p2, &p2, &alpha, &K21[0], &p2, &sigma22[0], &p2, &beta, &K12xK22_inv[0], &two );  
				
				// K121 = K12xK22_inv %*% K21													
				F77_NAME(dgemm)( &transN, &transN, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K21[0], &p2, &beta, &K121[0], &two );		
				// Finished (i,j) = 1------------------------------------------|

				a11      = K[i * dim1] - K121[0];	
				sum_diag = Dsjj * ( K022 - K121[3] ) - Ds[ij] * ( K121[1] + K121[2] );

				// nu_star = b + sum( Gf[,i] * Gf[,j] )
				nu_star = b1;
				//for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k];   
				for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k] * G[j * dim + k];   
				//nu_star = F77_NAME(ddot)( &dim, &G[0] + ixdim, &one, &G[0] + jxdim, &one );
				nu_star = 0.5 * nu_star;

				log_rate = ( G[ij] )   
					? 0.5 * log( 2.0 * Dsjj / a11 ) + lgammafn( nu_star + 0.5 ) - lgammafn( nu_star ) - 0.5 * ( Dsijj[ij] * a11 + sum_diag )
					: 0.5 * log( 0.5 * a11 / Dsjj ) - lgammafn( nu_star + 0.5 ) + lgammafn( nu_star ) + 0.5 * ( Dsijj[ij] * a11 + sum_diag );

				rates[counter++] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;
			}
		}
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

//----- saving result ---------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			for( i = 0; i < pxp ; i++ )
				if( G[i] ) p_links_Cpp[i] += weight_C;
			
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
			
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
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < pxp; i++ )
	{	
		p_links[i] = p_links_Cpp[i] / sum_weights;
		K_hat[i]   = K_hat_Cpp[i] / sum_weights;
	}
}
       
// ----------------------------------------------------------------------------|
// birth-death MCMC for Gaussian Graphical models  
// for case D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_map( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g,
			 int *b, int *b_star, double Ds[], int *print )
{
	int print_c = *print;
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	int count_all_g = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	bool this_one;

	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int i, j, k, ij, jj, counter, nu_star, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2, dim1 = dim + 1;

	double Dsjj, Dsij, sum_weights = 0.0, sum_diag, K022, a11, sigmajj_inv, weight_C;
	double log_rate, sum_rates;
	
	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideL = 'L';																	

	// Counting size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );             // K[e, -e]
	vector<double> sigma21( p2x2 );         // sigma[-e, e]
	vector<double> sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );        // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
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

	GetRNGstate();
	// main loop for birth-death MCMC sampling algorithm ----------------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
		counter = 0;	
		// STEP 1: calculating birth and death rates --------------------------|
		for( j = 1; j < dim; j++ )
		{
			jj   = j * dim1;
			Dsjj = Ds[jj];
			
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
			sigmajj_inv = - 1.0 / sigma[jj];
			F77_NAME(dsyr)( &sideL, &p1, &sigmajj_inv, &sigmaj12[0], &one, &sigmaj22[0], &p1 );

			for( i = 0; i < j; i++ )
			{
				ij = j * dim + i;

				// For (i,j) = 0 ----------------------------------------------|	
				sub_row_mins( K, &Kj12[0], &j, &dim );   // Kj12 = K[j, -j]  
				Kj12[ i ] = 0.0;                         // Kj12[1,i] = 0

				// Kj12xK22_inv = Kj12 %*% Kj22_inv
				F77_NAME(dsymv)( &sideL, &p1, &alpha, &sigmaj22[0], &p1, &Kj12[0], &one, &beta, &Kj12xK22_inv[0], &one );
				
				// K022 = Kj12xK22_inv %*% t(Kj12)
				K022 = F77_NAME(ddot)( &p1, &Kj12xK22_inv[0], &one, &Kj12[0], &one );			
				
				// For (i,j) = 1 ----------------------------------------------|
				sub_cols_mins( K, &K21[0], &i, &j, &dim );  // K21 = K[-e, e]  
				
				sub_matrices_inv( &sigma[0], &sigma11_inv[0], &sigma21[0], &sigma22[0], &i, &j, &dim );

				// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
				F77_NAME(dgemm)( &transN, &transN, &p2, &two, &two, &alpha, &sigma21[0], &p2, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
				F77_NAME(dgemm)( &transN, &transT, &p2, &p2, &two, &alpha1, &sigma21xsigma11_inv[0], &p2, &sigma21[0], &p2, &beta1, &sigma22[0], &p2 );

				// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
				F77_NAME(dgemm)( &transT, &transN, &two, &p2, &p2, &alpha, &K21[0], &p2, &sigma22[0], &p2, &beta, &K12xK22_inv[0], &two );  
				
				// K121 = K12xK22_inv %*% K21													
				F77_NAME(dgemm)( &transN, &transN, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K21[0], &p2, &beta, &K121[0], &two );		
				// Finished (i,j) = 1------------------------------------------|

				a11      = K[i * dim1] - K121[0];	
				sum_diag = Dsjj * ( K022 - K121[3] ) - Ds[ij] * ( K121[1] + K121[2] );

				// nu_star = b + sum( Gf[,i] * Gf[,j] )
				nu_star = b1;
				//for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k];   
				for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k] * G[j * dim + k];   
				nu_star = 0.5 * nu_star;   

				log_rate = ( G[ij] )   
					? 0.5 * log( 2.0 * Dsjj / a11 ) + lgammafn( nu_star + 0.5 ) - lgammafn( nu_star ) - 0.5 * ( Dsijj[ij] * a11 + sum_diag )
					: 0.5 * log( 0.5 * a11 / Dsjj ) - lgammafn( nu_star + 0.5 ) + lgammafn( nu_star ) + 0.5 * ( Dsijj[ij] * a11 + sum_diag );

				rates[counter] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;
				
				char_g[counter] = G[ij] + '0'; 
				counter++; 
			}
		}	
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

//----- saving result ---------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			//for( i = 0; i < pxp; i++ ) K_hat[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat[0], &one );			

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[count_all_g] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[count_all_g];
					all_graphs[count_all_g] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[count_all_g];
				all_graphs[count_all_g]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
    			
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
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < size_sample_graph; i++ ) 
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	
	for( i = 0; i < pxp; i++ )
		K_hat[i] /= sum_weights;
}
        
// ----------------------------------------------------------------------------|
// Multiple birth-death MCMC for Gaussian Graphical models  
// for D = I_p 
// it is for Bayesian model averaging
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_ma_multi_update( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 double K_hat[], double p_links[],
			 int *b, int *b_star, double Ds[], int *multi_update, int *print )
{
	int print_c = *print;
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	int multi_update_C = *multi_update;

	int selected_edge_i, selected_edge_j, selected_edge_ij;
	int i, j, k, ij, jj, counter, nu_star, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2, dim1 = dim + 1;

	double Dsjj, Dsij, sum_weights = 0.0, sum_diag, K022, a11, sigmajj_inv, weight_C;
	double log_rate, sum_rates;

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;

	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideL = 'L';																	

	// Count size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}

	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );             // K[e, -e]
	vector<double> sigma21( p2x2 );         // sigma[-e, e]
	vector<double> sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );        // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
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

	GetRNGstate();
// --- main loop for multiple birth-death MCMC sampling algorithm -------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % print_c < multi_update_C ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
		counter = 0;
		// STEP 1: calculating birth and death rates --------------------------|			
		for( j = 1; j < dim; j++ )
		{
			jj   = j * dim1;
			Dsjj = Ds[jj];
			
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
			sigmajj_inv = - 1.0 / sigma[jj];
			F77_NAME(dsyr)( &sideL, &p1, &sigmajj_inv, &sigmaj12[0], &one, &sigmaj22[0], &p1 );
									
			for( i = 0; i < j; i++ )
			{
				ij = j * dim + i;

				// For (i,j) = 0 ----------------------------------------------|	
				sub_row_mins( K, &Kj12[0], &j, &dim );   // Kj12 = K[j, -j]  
				Kj12[ i ] = 0.0;                         // Kj12[1,i] = 0

				// Kj12xK22_inv = Kj12 %*% Kj22_inv
				F77_NAME(dsymv)( &sideL, &p1, &alpha, &sigmaj22[0], &p1, &Kj12[0], &one, &beta, &Kj12xK22_inv[0], &one );
				
				// K022 = Kj12xK22_inv %*% t(Kj12)
				K022 = F77_NAME(ddot)( &p1, &Kj12xK22_inv[0], &one, &Kj12[0], &one );			
				
				// For (i,j) = 1 ----------------------------------------------|
				sub_cols_mins( K, &K21[0], &i, &j, &dim );  // K21 = K[-e, e]  
				
				sub_matrices_inv( &sigma[0], &sigma11_inv[0], &sigma21[0], &sigma22[0], &i, &j, &dim );

				// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
				F77_NAME(dgemm)( &transN, &transN, &p2, &two, &two, &alpha, &sigma21[0], &p2, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
				F77_NAME(dgemm)( &transN, &transT, &p2, &p2, &two, &alpha1, &sigma21xsigma11_inv[0], &p2, &sigma21[0], &p2, &beta1, &sigma22[0], &p2 );

				// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
				F77_NAME(dgemm)( &transT, &transN, &two, &p2, &p2, &alpha, &K21[0], &p2, &sigma22[0], &p2, &beta, &K12xK22_inv[0], &two );  
				
				// K121 = K12xK22_inv %*% K21													
				F77_NAME(dgemm)( &transN, &transN, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K21[0], &p2, &beta, &K121[0], &two );		
				// Finished (i,j) = 1------------------------------------------|

				a11      = K[i * dim1] - K121[0];	
				sum_diag = Dsjj * ( K022 - K121[3] ) - Ds[ij] * ( K121[1] + K121[2] );

				// nu_star = b + sum( Gf[,i] * Gf[,j] )
				nu_star = b1;
				//for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k];   
				for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k] * G[j * dim + k];   
				nu_star = 0.5 * nu_star;

				log_rate = ( G[ij] )   
					? 0.5 * log( 2.0 * Dsjj / a11 ) + lgammafn( nu_star + 0.5 ) - lgammafn( nu_star ) - 0.5 * ( Dsijj[ij] * a11 + sum_diag )
					: 0.5 * log( 0.5 * a11 / Dsjj ) - lgammafn( nu_star + 0.5 ) + lgammafn( nu_star ) + 0.5 * ( Dsijj[ij] * a11 + sum_diag );

				rates[counter++] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;
			}
		}	
		
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

//----- Saving result ---------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			for( i = 0; i < pxp ; i++ )
				if( G[i] ) p_links_Cpp[i] += weight_C;
			
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	

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
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	} // End of MCMC sampling algorithm ---------------------------------------|
	PutRNGstate();

	for( i = 0; i < pxp; i++ )
	{	
		p_links[i] = p_links_Cpp[i] / sum_weights;
		K_hat[i]   = K_hat_Cpp[i] / sum_weights;
	}
}
    
// ----------------------------------------------------------------------------|
// Multiple birth-death MCMC for Gaussian Graphical models  
// for D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// ----------------------------------------------------------------------------|
void ggm_bdmcmc_map_multi_update( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
			 int all_graphs[], double all_weights[], double K_hat[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g, int *counter_all_g,
			 int *b, int *b_star, double Ds[], int *multi_update, int *print )
{
	int print_c = *print;
	int multi_update_C = *multi_update;
	int iteration = *iter, burn_in = *burnin, b1 = *b;
	int count_all_g = *counter_all_g;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	bool this_one;

	int selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int i, j, k, ij, jj, counter, nu_star, one = 1, two = 2;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p1xp1 = p1 * p1, p2 = dim - 2, p2xp2 = p2 * p2, p2x2 = p2 * 2, dim1 = dim + 1;

	double Dsjj, Dsij, sum_weights = 0.0, sum_diag, K022, a11, sigmajj_inv, weight_C;
	double log_rate, sum_rates;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideL = 'L';																	

	// Counting size of notes
	int ip;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ip + j];
	}
	
	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	vector<double> Dsijj( pxp ); 
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
			
			// for calculating the birth/death rates
			ij        = j * dim + i;
			Dsij      = Ds[ij];
			Dsijj[ij] = Dsij * Dsij / Ds[j * dim + j]; 
		}

	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1xp1 );       // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );             // K[e, -e]
	vector<double> sigma21( p2x2 );         // sigma[-e, e]
	vector<double> sigma22( p2xp2 );        // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );        // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
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

	GetRNGstate();
	// main loop for multiple birth-death MCMC sampling algorithm -------------| 
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % print_c < multi_update_C ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
		counter = 0;
		// STEP 1: calculating birth and death rates
		for( j = 1; j < dim; j++ )
		{
			jj   = j * dim1;
			Dsjj = Ds[jj];
			
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
			sigmajj_inv = - 1.0 / sigma[jj];
			F77_NAME(dsyr)( &sideL, &p1, &sigmajj_inv, &sigmaj12[0], &one, &sigmaj22[0], &p1 );
									
			for( i = 0; i < j; i++ )
			{
				ij = j * dim + i;

				// For (i,j) = 0 ----------------------------------------------|	
				sub_row_mins( K, &Kj12[0], &j, &dim );   // Kj12 = K[j, -j]  
				Kj12[ i ] = 0.0;                         // Kj12[1,i] = 0

				// Kj12xK22_inv = Kj12 %*% Kj22_inv
				F77_NAME(dsymv)( &sideL, &p1, &alpha, &sigmaj22[0], &p1, &Kj12[0], &one, &beta, &Kj12xK22_inv[0], &one );
				
				// K022 = Kj12xK22_inv %*% t(Kj12)
				K022 = F77_NAME(ddot)( &p1, &Kj12xK22_inv[0], &one, &Kj12[0], &one );			
				
				// For (i,j) = 1 ----------------------------------------------|
				sub_cols_mins( K, &K21[0], &i, &j, &dim );  // K21 = K[-e, e]  
				
				sub_matrices_inv( &sigma[0], &sigma11_inv[0], &sigma21[0], &sigma22[0], &i, &j, &dim );

				// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
				F77_NAME(dgemm)( &transN, &transN, &p2, &two, &two, &alpha, &sigma21[0], &p2, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

				// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
				F77_NAME(dgemm)( &transN, &transT, &p2, &p2, &two, &alpha1, &sigma21xsigma11_inv[0], &p2, &sigma21[0], &p2, &beta1, &sigma22[0], &p2 );

				// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
				F77_NAME(dgemm)( &transT, &transN, &two, &p2, &p2, &alpha, &K21[0], &p2, &sigma22[0], &p2, &beta, &K12xK22_inv[0], &two );  
				
				// K121 = K12xK22_inv %*% K21													
				F77_NAME(dgemm)( &transN, &transN, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K21[0], &p2, &beta, &K121[0], &two );		
				// Finished (i,j) = 1------------------------------------------|

				a11      = K[i * dim1] - K121[0];	
				sum_diag = Dsjj * ( K022 - K121[3] ) - Ds[ij] * ( K121[1] + K121[2] );

				// nu_star = b + sum( Gf[,i] * Gf[,j] )
				nu_star = b1;
				//for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k];   
				for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k] * G[j * dim + k];   
				nu_star = 0.5 * nu_star;

				log_rate = ( G[ij] )   
					? 0.5 * log( 2.0 * Dsjj / a11 ) + lgammafn( nu_star + 0.5 ) - lgammafn( nu_star ) - 0.5 * ( Dsijj[ij] * a11 + sum_diag )
					: 0.5 * log( 0.5 * a11 / Dsjj ) - lgammafn( nu_star + 0.5 ) + lgammafn( nu_star ) + 0.5 * ( Dsijj[ij] * a11 + sum_diag );

				rates[counter] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;
				
				char_g[counter] = G[ij] + '0'; 
				counter++; 
			}
		}	
		
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

//----- Saving result ---------------------------------------------------------|	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			//for( i = 0; i < pxp; i++ ) K_hat[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat[0], &one );			

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[count_all_g] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[count_all_g];
					all_graphs[count_all_g] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[count_all_g];
				all_graphs[count_all_g]          = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
			
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
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	} // End of MCMC sampling algorithm ---------------------------------------| 
	PutRNGstate();

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	*counter_all_g = count_all_g;
	
	for( i = 0; i < pxp; i++ ) K_hat[i] /= sum_weights;		
}
              
} // End of exturn "C"
