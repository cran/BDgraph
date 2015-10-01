# plot of graph size to check the convergency of BDMCMC algorithm
traceplot = function( x, acf = FALSE, pacf = FALSE, main = NULL, ... )
{
	if( !is.null( x $ phat ) ) stop( "Function needs output of 'bdgraph' with option save.all = TRUE" )  

	sampleGraphs     = x $ sampleGraphs
    allGraphs        = x $ allGraphs
	graphWeights     = x $ graphWeights
	sizeSampleGraphs = sapply( sampleGraphs, function(x) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )  
	sizeAllGraphs    = sizeSampleGraphs[ allGraphs ]
	which_G_max      = which( max( graphWeights ) == graphWeights )
	size_selected_G  = sizeAllGraphs[ which_G_max ] 

	if( is.null( main ) ) main = "Trace of graph size"
	
	x_vec = 1 : length( allGraphs )
	
	if ( acf == FALSE & pacf == FALSE )
	{
		plot( x = x_vec, sizeAllGraphs, type = "l", main = main, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_G, col = "red" )	   
	}
	
	if ( acf == TRUE & pacf == TRUE )
	{
		op = par( mfrow = c( 2, 2 ), pty = "s" )  
		plot( x = x_vec, sizeAllGraphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_G, col = "red" )	  
		acf( sizeAllGraphs,  main = "ACF for graph size" )
		pacf( sizeAllGraphs, main = "PACF for graph size" )
		par( op )
	}
	
	if ( acf == TRUE & pacf == FALSE )
	{
		op <- par( mfrow = c( 1, 2 ), pty = "s" ) 
		plot( x = x_vec, sizeAllGraphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_G, col = "red" )	  
		acf( sizeAllGraphs, main = "ACF for graph size" )
		par( op )
	}
	
	if ( acf == FALSE & pacf == TRUE )
	{
		op <- par( mfrow = c( 1, 2 ), pty = "s" ) 
		plot( x = x_vec, sizeAllGraphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_G, col = "red" )	  
		pacf( sizeAllGraphs, main = "PAIC for graph size" )
		par(op)
	}		
}  
      
