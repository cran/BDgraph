# function for ROC plot
outRoc = function( G, prob, cut )
{
	G[ lower.tri( G, diag = TRUE ) ]       <- 0
	prob[ lower.tri( prob, diag = TRUE ) ] <- 0
	p          = nrow(prob)
	sumEdges   = sum(G)
	sumNoEdges = p * ( p - 1 ) / 2 - sumEdges
	
	tp = c( 1, rep( 0, cut ) )
	fp = tp

	cutPoint = ( 0:cut ) / cut
	
	for ( i in 2:cut )
	{
		# checking for cut pints
		estG = 0 * G
		estG[prob > cutPoint[i]] = 1

		tp.all <- ( G != 0 ) * ( estG != 0 ) 
		fp.all <- ( G == 0 ) * ( estG != 0 ) 	
		tp[i]  <- sum( tp.all ) / sumEdges
		fp[i]  <- sum( fp.all ) / sumNoEdges
	}
	
	return( list( tp = tp, fp = fp ) )
}
    
# To plot ROC curve
plotroc = function( G, prob, prob2 = NULL, cut = 20, smooth = FALSE )
{
    if ( class(G)     == "sim" ) G <- as.matrix( G $ G )
    if ( class(prob)  == "bdgraph" )
    {
		phat = prob $ phat
		if( is.null( phat ) ) phat = phat( prob, round = 10 )
		prob = as.matrix( phat )
	} 
    
    output = outRoc( G = G, prob = prob, cut = cut )
    x      = output $ fp
    y      = output $ tp
 	
	if ( smooth == TRUE )
	{
		fit = smooth.spline( x = x, y = y )
		x   = c( 0, fit $ x )
		y   = c( 0, fit $ y )
	}
	
	# par( mar = c( 3.8, 4.2, 1.8, 1 ) )
    plot( x = x, y = y, type = "l", col = "black", lty = 1, cex.lab = 1.3, cex.main = 2, cex.axis = 1.2,
         main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate", ylim = c(0,1) )
  
    if( !is.null( prob2 ) )
    {
		if ( class(prob2)  == "bdgraph" )
		{
			phat2 = prob2 $ phat
			if( is.null( phat2 ) ) phat2 = phat( prob2, round = 10 )
			prob2 = as.matrix( phat2 )
		} 

        output2 = outRoc( G = G, prob = prob2, cut = cut )
		x2      = output2 $ fp
		y2      = output2 $ tp

		if ( smooth == TRUE )
		{
			fit2 = smooth.spline( x = x2, y = y2 )
			x2   = c( 0, fit2 $ x )
			y2   = c( 0, fit2 $ y )
		}
		
        points( x = x2, y = y2, type = "l", col = "blue", lty = 2, lw = 2 )
    }
}
       
