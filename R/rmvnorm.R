## ------------------------------------------------------------------------------------------------|
# Data generator from multivarate normal distribution X N_p( mu, sig )
## ------------------------------------------------------------------------------------------------|
rmvnorm = function( n = 10, mean = rep( 0, length = ncol( sigma ) ), sigma = diag( length( mean ) ) )
{
    if( !isSymmetric( sigma, tol = sqrt( .Machine$double.eps ), check.attributes = FALSE ) ) 
        stop( "sigma must be a symmetric matrix" )
    
    sigma <- as.matrix( sigma )
    p     <- nrow( sigma )
    if( typeof( mean ) == "double" ) mean <- rep( mean, p )
    if( length( mean ) != nrow( sigma ) ) stop( "mean and sigma have non-conforming size" )
    
    #--- generate multivariate normal data ----------------------------------------------------|
    chol_sig <- chol( sigma )
    z        <- matrix( rnorm( p * n ), p, n )
    data     <- t( chol_sig ) %*% z + mean
    data     <- t( data )

    return( data )
}
   