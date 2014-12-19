# R code for sampling from G-Wishart AND Wishart distribution
################################################################################
# sampling from G-Wishart distribution
rgwish = function( n = 1, G = NULL, b = 3, D = NULL )
{
	if ( is.null(G) ) stop( "Adjacency matrix G should be determined" )
	G <- as.matrix(G)
	if ( sum( (G == 1) * (G == 0) ) != 0 ) stop( "Elements of matrix G should be zero or one" )	
	if ( !isSymmetric(G) )
	{
		G[ lower.tri( G, diag(TRUE) ) ] <- 0
		G  = G + t(G)
	}
	
	p <- nrow(G)  
	
	if ( is.null(D) ) D <- diag(p)	
	Ti = chol( solve(D) ) 
	
	samples <- array( 0, c( p, p, n ) )
	K  = matrix( 0, p, p )
	
	for ( i in 1 : n )
	{
		# rgwish ( int G[], double T[], double K[], int *b, int *p )
		result = .C( "rgwish", as.integer(G), as.double(Ti), K = as.double(K), 
					 as.integer(b), as.integer(p)
					 , PACKAGE = "BDgraph" )
		samples[,,i] = matrix ( result $ K, p, p ) 		
	}	

	return( samples )   
}
