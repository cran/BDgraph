# R code for sampling from COMPLEX G-Wishart
################################################################################
# sampling from COMPLEX G-Wishart distribution
rgcwish = function( n = 1, adj.g = NULL, b = 3, D = NULL )
{
	if( b <= 2 ) stop( "In complex G-Wishart distribution parameter 'b' has to be more than 2" )
	if( is.null( adj.g ) ) stop( "Adjacency matrix should be determined" )

	G <- as.matrix( adj.g )
	if( sum( ( G == 1 ) * ( G == 0 ) ) != 0 ) stop( "Elements of matrix G should be zero or one" )	

	if( !isSymmetric(G) )
	{
		G[ lower.tri( G, diag( TRUE ) ) ] <- 0
		G  = G + t(G)
	}
	
	p <- nrow(G)  
	
	if( is.null(D) ) D <- diag(p)
	
	if( dim( D )[1] != p ) stop( "Dimension of matrix G and D must to be the same." )
		
	inv_D     = solve( D )
	row1      = cbind( Re(inv_D), -Im(inv_D) )
	row2      = cbind( Im(inv_D), Re(inv_D) )
	sig       = rbind( row1, row2 ) / 2
	Ls        = t( chol(sig) )
	samples   = array( 0, c( p, p, n ) )
	K         = matrix( 0, p, p )
	threshold = 1e-8
	
	for( i in 1 : n )
	{
		result       = .C( "rgcwish", as.integer(G), as.double(Ls), K = as.complex(K), as.integer(b), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
		samples[,,i] = matrix( result $ K, p, p ) 		
	}	

	return( samples )   
}
  
