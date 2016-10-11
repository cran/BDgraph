#include "matrix.h"

// Takes square matrix A (p x p) and 
// retrieves square sub_matrix B (p_sub x p_sub), dictated by vector sub
void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int ixp, subixp, psub = *p_sub, pdim = *p;
	
	for( int i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[i] * pdim;
		
		for( int j = 0; j < psub; j++ )
			sub_A[ixp + j] = A[subixp + sub[j]]; 
	}
}

// Takes square matrix A (p x p) and 
// retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
void sub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int i, subj = *sub, pdim = *p, subxp = subj * pdim;

	memcpy( sub_A,        A + subxp,            sizeof( double ) * subj );		
	memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R ONLY FOR SYMMETRIC MATRICES
void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col;
	int sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}
}

// Takes symmatric matrix A (p x p) and 
// retrieves A_jj, A12(1x(p-1)), A21((p-1)x1), and A22((p-1)x(p-1))
// Like A11=A[j, j], A12=A[j, -j], and A22=A[-j, -j] in R
void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, j, ixpdim, ij, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim;

	memcpy( A12,        A + subxp,            sizeof( double ) * psub );	
	memcpy( A12 + psub, A + subxp + psub + 1, sizeof( double ) * ( pdim - psub - 1 ) );	

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		
		for( j = 0; j < psub; j++ )
			A22[j * p1 + i] = A[ixpdim + j];

		for( j = psub + 1; j < pdim; j++ )
		{
			ij = ixpdim + j;
			A22[(j - 1) * p1 + i] = A[ij];
			A22[i * p1 + j - 1]   = A[ij];
		}
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
			
		for( j = psub + 1; j < pdim; j++ )
			A22[(j - 1) * p1 + i - 1] = A[ixpdim + j];
	}
}

// Takes square matrix A (p x p) and 
// retrieves A11(2x2), A12(2x(p-2)), and A22((p-2)x(p-2))
// Like A11=A[e, e], A12=A[e, -e], and A22=A[-e, -e] in R
void sub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, j, ixp, ij, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[0] = A[sub0 * pdim + sub0];
	A11[1] = A[sub0 * pdim + sub1];
	A11[2] = A11[1];                   // for symmetric matrices
	A11[3] = A[sub1 * pdim + sub1];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[i + i]     = A[ixp + sub0];
		A12[i + i + 1] = A[ixp + sub1];
	
		for( j = 0; j < sub0; j++ )
			A22[j * p2 + i] = A[ixp + j];

		for( j = sub0 + 1; j < sub1; j++ )
		{
			ij = ixp + j;
			A22[(j - 1) * p2 + i] = A[ij];
			A22[i * p2 + j - 1]   = A[ij];
		}
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[(j - 2) * p2 + i] = A[ij];
			A22[i * p2 + j - 2]   = A[ij];
		}
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		
		A12[i + i - 2] = A[ixp + sub0];
		A12[i + i - 1] = A[ixp + sub1];
	
		for( j = sub0 + 1; j < sub1; j++ )
			A22[(j - 1) * p2 + i - 1] = A[ixp + j];
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[(j - 2) * p2 + i - 1] = A[ij];
			A22[(i - 1) * p2 + j - 2] = A[ij];
		}
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
				
		A12[i + i - 4]     = A[ixp + sub0];
		A12[i + i - 3]       = A[ixp + sub1];
		
		for( j = sub1 + 1; j < pdim; j++ )
			A22[(j - 2) * p2 + i - 2] = A[ixp + j];
	}
}
   
// Multiplies (p_i x p_k) matrix by (p_k x p_j) matrix to give (p_i x p_j) matrix
// C := A %*% B
void multiply_matrix( double A[], double B[], double C[], int *p_i, int *p_j, int *p_k )
{
	double alpha = 1.0, beta  = 0.0;
	char trans   = 'N';																	
	F77_NAME(dgemm)( &trans, &trans, p_i, p_j, p_k, &alpha, A, p_i, B, p_k, &beta, C, p_i );
}

// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
void inverse( double A[], double A_inv[], int *p )
{
	int info, dim = *p;
	char uplo = 'U';

	// creating an identity matrix
	for( int i = 0; i < dim; i++ )
		for( int j = 0; j < dim; j++ )
			A_inv[j * dim + i] = (i == j);
	
	// LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
	F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info );
}

// inverse function for symmetric (2 x 2)
void inverse_2x2( double B[], double B_inv[] )
{
	double detB = B[0] * B[3] - B[1] * B[1];
	B_inv[0]    = B[3] / detB;
	B_inv[1]    = - B[1] / detB;
	B_inv[2]    = B_inv[1];
	B_inv[3]    = B[0] / detB;
}
    
// Cholesky decomposition of symmetric positive-definite matrix
// WARNING: Matrix you pass is overwritten with the result
// A = U' %*% U
void cholesky( double A[], double U[], int *p )
{
	char uplo = 'U';
	int j, info, dim = *p, i, pxp = dim * dim;
	
	memcpy( U, A, sizeof( double ) * pxp );	
	
	F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info );	
	
	for( i = 0; i < dim; i++ )
		for( j = 0; j < i; j++ )
			U[j * dim + i] = 0.0;
}
  
// To select an edge for BDMCMC algorithm  
void select_edge( long double rates[], int *index_selected_edge, long double *sum_rates, int *qp )
{
	long double random_value;
	int qp_star = *qp;

	vector<long double> sum_sort_rates( qp_star );
	sum_sort_rates[0] = rates[0];
	for ( int i = 1; i < qp_star; i++ )
		sum_sort_rates[i] = sum_sort_rates[ i - 1 ] + rates[i];
	
	*sum_rates = sum_sort_rates[qp_star - 1];
	random_value = *sum_rates * runif( 0, 1 );

	// To start, find the subscript of the middle position.
	int position;
	int lower_bound = 0;
	int upper_bound = qp_star - 1;
	position = upper_bound / 2;      // ( lower_bound + upper_bound ) / 2;

	while( upper_bound - lower_bound > 1 )
	{
		if ( sum_sort_rates[position] > random_value )    
			upper_bound = position;    
		else                                                
			lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	*index_selected_edge = ( sum_sort_rates[position] < random_value ) ? ++position : position;
} 
    
// To simultaneously select multiple edges for BDMCMC algorithm  
void select_multi_edges( long double rates[], int index_selected_edges[], int *size_index, long double *sum_rates, int *multi_update, int *qp )
{
	int qp_star = *qp, lower_bound, upper_bound, position;
	long double max_bound, random_value;

	vector<long double> sum_sort_rates( qp_star );
	sum_sort_rates[0] = rates[0];
	for ( int i = 1; i < qp_star; i++ )
		sum_sort_rates[i] = sum_sort_rates[ i - 1 ] + rates[i];
	
	max_bound = sum_sort_rates[qp_star - 1];
	
	// ---------- for first edge ----------------------------------------------|
	// To start, find the subscript of the middle position.
	lower_bound = 0;
	upper_bound = qp_star - 1;
	position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

	random_value = max_bound * runif( 0, 1 );

	while( upper_bound - lower_bound > 1 )
	{
		if ( sum_sort_rates[position] > random_value )    
			upper_bound = position;    
		else                                                
			lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	if ( sum_sort_rates[position] < random_value ) ++position;
	index_selected_edges[0] = position;
	// ------------------------------------------------------------------------|

	int counter = 1;
	int same;
	for ( int it = 0; it < 200 * *multi_update; it++ )
	{
		if ( counter == *multi_update ) break;
		
		random_value = max_bound * runif( 0, 1 );
	
		// To start, find the subscript of the middle position.
		lower_bound = 0;
		upper_bound = qp_star - 1;
		position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

		while( upper_bound - lower_bound > 1 )
		{
			if ( sum_sort_rates[position] > random_value )    
				upper_bound = position;    
			else                                                
				lower_bound = position;     
			
			position = ( lower_bound + upper_bound ) / 2;
		}
		
		if ( sum_sort_rates[position] < random_value ) ++position;
		
		same = 0;
		for ( int i = 0; i < counter; i++ )
			if( index_selected_edges[i] == position )
				++same;

		if ( same == 0 ) index_selected_edges[counter++] = position;
	}
	
	*size_index = counter;
	*sum_rates  = max_bound;
} 
         
// computing birth/death rate or alpha for element (i,j)
// it is for double Metropolis-Hasting algorihtms
void log_H_ij( double K[], double sigma[], double *log_Hij, int *selected_edge_i, int *selected_edge_j,
               double Kj22_inv[], double Kj12[], double Kj12xK22_inv[], double *K022, double K12[], double K22_inv[], double K12xK22_inv[], double K121[], 
               double sigmaj12[], double sigmaj22[], double sigma11[], double sigma12[], double sigma22[], double sigma11_inv[], double sigma21xsigma11_inv[], double sigma2112[],
               int *dim, int *p1, int *p2, int *p2xp2, int *jj,
               double *Dsijj, double *Dsij, double *Dsjj )
{
	int p1_here = *p1, p2_here = *p2, dim_here = *dim, one = 1, two = 2;
	double alpha = 1.0, beta = 0.0;
	char transT = 'T', transN = 'N';																	
	
	double sigmaj11 = sigma[*jj];        // sigma[j, j]  
	sub_matrices1( sigma, sigmaj12, sigmaj22, selected_edge_j, &dim_here );

	// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
	for( int row = 0; row < p1_here; row++ )
		for( int col = 0; col < p1_here; col++ )
		{
			int rowCol = col * p1_here + row;
			Kj22_inv[rowCol] = sigmaj22[rowCol] - sigmaj12[row] * sigmaj12[col] / sigmaj11;
		}

	// For (i,j) = 0 ------------------------------------------------------|
	sub_row_mins( K, Kj12, selected_edge_j, &dim_here );  // K12 = K[j, -j]  
	Kj12[ *selected_edge_i ] = 0.0;                        // K12[1,i] = 0

	// K12 %*% K22_inv
	F77_NAME(dgemm)( &transN, &transN, &one, &p1_here, &p1_here, &alpha, Kj12, &one, Kj22_inv, &p1_here, &beta, Kj12xK22_inv, &one );

	// K022  <- K_12 %*% solve( K0[-j, -j] ) %*% t(K_12) 
	F77_NAME(dgemm)( &transN, &transT, &one, &one, &p1_here, &alpha, Kj12xK22_inv, &one, Kj12, &one, &beta, K022, &one );			

	// For (i,j) = 1 ------------------------------------------------------|
	// K12 = K[e, -e] 
	sub_rows_mins( K, K12, selected_edge_i, selected_edge_j, &dim_here );  
	
	sub_matrices( sigma, sigma11, sigma12, sigma22, selected_edge_i, selected_edge_j, &dim_here );

	// solve( sigma[e, e] )
	inverse_2x2( sigma11, sigma11_inv );

	// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
	F77_NAME(dgemm)( &transT, &transN, &p2_here, &two, &two, &alpha, sigma12, &two, sigma11_inv, &two, &beta, sigma21xsigma11_inv, &p2_here );

	// sigma21xsigma11_inv %*% sigma12
	F77_NAME(dgemm)( &transN, &transN, &p2_here, &p2_here, &two, &alpha, sigma21xsigma11_inv, &p2_here, sigma12, &two, &beta, sigma2112, &p2_here );

	// solve( K[-e, -e] ) = sigma22 - sigma2112
	for( int k = 0; k < *p2xp2 ; k++ ) 
		K22_inv[k] = sigma22[k] - sigma2112[k];	
	
	// K12 %*% K22_inv
	F77_NAME(dgemm)( &transN, &transN, &two, &p2_here, &p2_here, &alpha, K12, &two, K22_inv, &p2_here, &beta, K12xK22_inv, &two );

	// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e])															
	F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2_here, &alpha, K12xK22_inv, &two, K12, &two, &beta, K121, &two );		
	// Finished (i,j) = 1--------------------------------------------------|

	double a11      = K[*selected_edge_i * dim_here + *selected_edge_i] - K121[0];	
	double sum_diag = *Dsjj * ( *K022 - K121[3] ) - *Dsij * ( K121[1] + K121[2] );

	// Dsijj = Dsii - Dsij * Dsij / Dsjj;
	*log_Hij = ( log( static_cast<double>(*Dsjj) ) - log( static_cast<double>(a11) ) + *Dsijj * a11 - sum_diag ) / 2;
}    

// -------------- NEW for Lang codes -------------------------------------------
// For Hermitian matrix
void Hsub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int i, l = 0, subj = *sub, pdim = *p, subxp = subj * pdim;

	for( i = 0; i < subj; i++ )
		sub_A[l++] = -A[subxp + i];
	
	for( i = subj + 1; i < pdim; i++ )
		sub_A[l++] = -A[subxp + i];
}
      
// For Hermitian matrix
void Hsub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col;
	int sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[l++] = -A[sub0p + i]; 
		sub_A[l++] = -A[sub1p + i]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[l++] = -A[sub0p + i]; 
		sub_A[l++] = -A[sub1p + i]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[l++] = -A[sub0p + i]; 
		sub_A[l++] = -A[sub1p + i]; 
	}
}
       
// sub_matrices1 for Hermitian matrix
void Hsub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, j, ixpdim, ij, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;

	for( i = 0; i < psub; i++ )
		A12[i] = -A[subxp + i];
	for( i = psub; i < pdim - 1; i++ )
		A12[i] = -A[subxp + i + 1];

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		memcpy( A22 + i * p1, A + ixpdim, sizeof( double ) * psub );
		memcpy( A22 + i * p1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
		memcpy( A22 + ( i - 1 ) * p1, A + ixpdim, sizeof( double ) * psub);
                memcpy( A22 + ( i - 1 ) * p1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}
}
        
// sub_matrices for Hermitian matrix
void Hsub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, i1, i2, j, ixp, ij, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[0] = A[sub0 * pdim + sub0];
	A11[1] = A[sub0 * pdim + sub1];
	A11[2] = -A11[1];                   // for symmetric matrices
	A11[3] = A[sub1 * pdim + sub1];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[i + i]     = A[ixp + sub0];
		A12[i + i + 1] = A[ixp + sub1];

		memcpy( A22 + i * p2, A + ixp, sizeof( double ) * sub0 );
		memcpy( A22 + i * p2 + sub0, A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		i1 = i - 1;

		A12[i + i - 2] = A[ixp + sub0];
		A12[i + i - 1] = A[ixp + sub1];

		memcpy( A22 + i1 * p2, A + ixp, sizeof( double ) * sub0 );
		memcpy( A22 + i1 * p2 + sub0, A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i1 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
		i2 = i - 2;
				
		A12[i + i - 4]     = A[ixp + sub0];
		A12[i + i - 3]     = A[ixp + sub1];

		memcpy( A22 + i2 * p2, A + ixp, sizeof( double ) * sub0 );
		memcpy( A22 + i2 * p2 + sub0, A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i2 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );		
	}
}
   
// inverse function for Hermitian (2 x 2)
void cinverse_2x2( double r_B[], double i_B[], double r_B_inv[], double i_B_inv[] )
{
	double r_det = r_B[0] * r_B[3] - i_B[0] * i_B[3] - ( r_B[1] * r_B[1] + i_B[1] * i_B[1] );
	double i_det = r_B[0] * i_B[3] + i_B[0] * r_B[3];
	double mod   = r_det * r_det + i_det * i_det;
	
	r_B_inv[0] =  ( r_B[3] * r_det + i_B[3] * i_det ) / mod;
	i_B_inv[0] =  ( r_det * i_B[3] - r_B[3] * i_det ) / mod;
	r_B_inv[1] = -( r_B[1] * r_det + i_B[1] * i_det ) / mod;
	i_B_inv[1] = -( r_det * i_B[1] - r_B[1] * i_det ) / mod;
	r_B_inv[2] = -( r_B[1] * r_det - i_B[1] * i_det ) / mod;
	i_B_inv[2] =  ( r_det * i_B[1] + r_B[1] * i_det ) / mod;
	r_B_inv[3] =  ( r_B[0] * r_det + i_B[0] * i_det ) / mod;
	i_B_inv[3] =  ( r_det * i_B[0] - r_B[0] * i_det ) / mod;
}

// For generating scale-free graphs: matrix G (p x p) is an adjacency matrix
void scale_free( int *G, int *p )
{
	int i, j, tmp, total, dim = *p, p0 = 2;
	double random_value;
	std::vector<int> size_a( dim ); 

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
	   random_value = (double) total * runif( 0, 1 );
	   
		tmp = 0;
		j   = 0;
		
		while( tmp < random_value && j < i ) 
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
