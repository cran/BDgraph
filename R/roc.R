# To compare the result
roc = function ( G, est ) 
{
	G   <- as.matrix(G)     # G is the adjacency matrix of true graph 
	est <- as.matrix(est)   # est is the adjacency matrix of estimated graph 
	G[ lower.tri( G, diag = TRUE ) ]     <- 0
	est[ lower.tri( est, diag = TRUE ) ] <- 0
	p      <- nrow(G)
	
	tp.all <- ( G != 0 ) * ( est != 0 ) 
	fp.all <- ( G == 0 ) * ( est != 0 ) 
	fn.all <- ( G != 0 ) * ( est == 0 ) 
	tn.all <- ( G == 0 ) * ( est == 0 )
	
	tp     <- sum( tp.all )
	fp     <- sum( fp.all )
	fn     <- sum( fn.all )
	tn     <- sum( tn.all[ upper.tri( tn.all == 1 ) ] )
	
	# positive predictive value  
	Precision <- tp / ( tp + fp ) 
	
	if ( is.na(Precision) ) Precision <- 0
	# Precision is the probability that a randomly selected link is relevant
	Recall <- tp / ( tp + fn ) # also called TPR
	
	if ( is.na(Recall) ) Recall <- 0
	# Recall is the probability that a randomly selected relevant link 
	# is retrieved in a search 
	FPR <- fp / ( fp + tn ) # False positive rate
	
	if ( is.na(FPR) ) FPR <- 0
	Accuracy <- ( tp + tn ) / (tp + tn + fp + fn)
	
	if ( is.na(Accuracy) ) Accuracy <- 0
	# Specificity <- tn / (tn + fp) # or 1 - false positive rate 
	# Sensitivity <- tp / (tp + fn) # or true positive rate
	# F1score <- 2 * (Precision * Recall) / (Precision + Recall)
	F1score <- ( 2 * tp ) / ( 2 * tp + fp + fn )
	if ( is.na(F1score) ) F1score <- 0
	# harmonic mean of precision and recall, called F-measure or balanced F-score:
	roc.matrix <- matrix( c(tp, tn, fp, fn, Recall, FPR, Accuracy, F1score, Precision), 9, 1 )
	
	return( round( roc.matrix, 3 ) )
}
# function for ROC plot
outRoc = function( G, prob, cut )
{
	G[ lower.tri( G, diag = TRUE ) ]     <- 0
	prob[ lower.tri( prob, diag = TRUE ) ] <- 0
	p = nrow(prob)
	pp = p * ( p - 1 ) / 2
	sumEdges = sum(G)
	sumNoEdges = pp - sumEdges
	
	tp = c( rep( 0, cut + 1 ) )
	tp[1] = 1
	fp = c( rep( 0, cut + 1 ) )
	fp[1] = 1

	cutPoint = (0:cut) / cut
	
	for ( i in 2:cut )
	{
		# checking for cut pints
		estG = matrix( 0, p, p )
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
    if ( class(G)     == "simulate" ) G     <- as.matrix( G $ G )
    if ( class(prob)  == "bdgraph" )  prob  <- as.matrix( phat( prob, round = 10 ) ) 
    
    output = outRoc( G = G, prob = prob, cut = cut )
    x      = output $ fp
    y      = output $ tp
 	
	if ( smooth == TRUE )
	{
		fit = smooth.spline( x = x, y = y )
		x   = c( 0, fit $ x )
		y   = c( 0, fit $ y )
	}
	
	par( mar = c( 3.8, 4.2, 1.8, 1 ) )
    plot( x = x, y = y, type = "l", col = "black", lty = 1, cex.lab = 1.6, cex.main = 2.5, cex.axis = 1.7,
         main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate", ylim = c(0,1) )
  
    if( !is.null(prob2) )
    {
        if ( class(prob2)  == "bdgraph" ) prob2 <- as.matrix( phat( prob2, round = 10 ) ) 
        output2 = outRoc( G = G, prob = prob2, cut = cut )
		x2      = output2 $ fp
		y2      = output2 $ tp

		if ( smooth == TRUE )
		{
			fit2 = smooth.spline( x = x2, y = y2 )
			x2   = c( 0, fit2 $ x )
			y2   = c( 0, fit2 $ y )
		}
		
        points( x = x2, y = y2, type = "l", col = "blue", lty = 2, lw = 2,
                main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate", ylim = c(0,1) )
    }
}
