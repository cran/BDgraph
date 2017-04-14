# sampling from COMPLEX Wishart distribution
rcwish = function( n = 1, p = 2, b = 3, D = diag(p) )
{
	inv_D   = solve( D )
	row1    = cbind( Re(inv_D), -Im(inv_D) )
	row2    = cbind( Im(inv_D), Re(inv_D) )
	sig     = rbind( row1, row2 ) / 2
	Ls      = t( chol(sig) )
	samples = array( 0, c( p, p, n ) )
	K       = matrix( 0, p, p )
	
	for ( i in 1 : n )
	{
		result       = .C( "rcwish_c", as.double(Ls), K = as.complex(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
		samples[,,i] = matrix( result $ K, p, p ) 		
	}	

	return( samples )   
}
    
