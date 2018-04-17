## ----------------------------------------------------------------------------|
# sampling from Wishart distribution
## ----------------------------------------------------------------------------|
rwish = function( n = 1, p = 2, b = 3, D = diag(p) )
{
	Ti      = chol( solve( D ) ) 
	K       = matrix( 0, p, p )

	if( n > 1 )
	{
		samples = array( 0, c( p, p, n ) )
		
		for ( i in 1 : n )
		{
			result       = .C( "rwish_c", as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
			samples[,,i] = matrix( result $ K, p, p ) 		
		}
	}else{
		result  = .C( "rwish_c", as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
		samples = matrix( result $ K, p, p ) 				
	}	

	return( samples )   
}
    
