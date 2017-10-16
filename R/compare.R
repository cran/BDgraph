# To compare the result according to the true graph
compare = function( sim.obj, bdgraph.obj, bdgraph.obj2 = NULL, bdgraph.obj3 = NULL, colnames = NULL, vis = FALSE ) 
{
	if( is.matrix( sim.obj ) )      G    = sim.obj
	if( is.matrix( bdgraph.obj ) )  est  = bdgraph.obj
	
	if( !is.matrix( sim.obj )      && ( class( sim.obj )      == "sim"     ) )  G    <- sim.obj $ G 
	if( !is.matrix( bdgraph.obj )  && ( class( bdgraph.obj )  == "bdgraph" ) )  est  <- select( bdgraph.obj ) 
	if( class( bdgraph.obj )  == "select"  )  est  <- bdgraph.obj $ refit
   
	G   = as.matrix( G )        # G is the adjacency matrix of true graph 
	est = as.matrix( est )      # est is the adjacency matrix of estimated graph 
	p   = nrow( G )
	G[ lower.tri( G, diag = TRUE ) ]     = 0
	est[ lower.tri( est, diag = TRUE ) ] = 0
   	   
	result = matrix( 1, 8, 2 )
	result[ , 2 ] = compute_measures( G = G, est_G = est )
   
    if( !is.null( bdgraph.obj2 ) )
    {
		if( is.matrix( bdgraph.obj2 ) ) est2 = bdgraph.obj2
		if( !is.matrix( bdgraph.obj2 ) && ( class( bdgraph.obj2 ) == "bdgraph" ) )  est2 <- select( bdgraph.obj2 ) 
		if( class( bdgraph.obj2 ) == "select"  )  est2 <- bdgraph.obj2 $ refit

		est2 = as.matrix( est2 )       
		est2[ lower.tri( est2, diag = TRUE ) ] = 0

		result = cbind( result, compute_measures( G = G, est_G = est2 ) )
	} 
	
    if( !is.null( bdgraph.obj3 ) )
    { 
		if( is.matrix( bdgraph.obj3 ) ) est3 = bdgraph.obj3
		if( !is.matrix( bdgraph.obj3 ) && ( class( bdgraph.obj3 ) == "bdgraph" ) )  est3 <- select( bdgraph.obj3 ) 
		if( class( bdgraph.obj3 ) == "select"  )  est3 <- bdgraph.obj3 $ refit
		  
		est3 = as.matrix( est3 )       
		est3[ lower.tri( est3, diag = TRUE ) ] = 0

		result = cbind( result, compute_measures( G = G, est_G = est3 ) )
	} 
	
	result[ c( 3, 4 ), 1 ] = 0
	result[ 1, 1 ] = sum( G )
	result[ 2, 1 ] = p * ( p - 1 ) / 2 - result[ 1, 1 ]  

	result[ is.na( result ) ] = 0

	if( is.null( colnames ) ) 
	{
		colnames = c( "True", "estimate" )
		if( !is.null( bdgraph.obj2 ) ) colnames = c( colnames, "estimate2" )
		if( !is.null( bdgraph.obj3 ) ) colnames = c( colnames, "estimate3" )
	}
		
	colnames( result ) <- colnames

	rownames( result ) <- c( "true positive", "true negative", "false positive", "false negative", 
                             "F1-score", "specificity", "sensitivity", "MCC" )
				
   if( vis )
   {
		G_igraph   <- graph.adjacency( G,   mode = "undirected", diag = FALSE )
		est_igraph <- graph.adjacency( est, mode = "undirected", diag = FALSE )
		if ( p < 20 ) sizev = 15 else sizev = 2

		row_plot = ifelse( is.null( bdgraph.obj2 ), 1, 2 )
		op       = par( mfrow = c( row_plot, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )

		plot.igraph( G_igraph,   layout = layout.circle, main = colnames[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		plot.igraph( est_igraph, layout = layout.circle, main = colnames[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		 
		if( !is.null( bdgraph.obj2 ) )
		{
			est2_igraph <- graph.adjacency( as.matrix(est2), mode = "undirected", diag = FALSE )
			plot.igraph( est2_igraph, layout = layout.circle, main = colnames[3], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		if( !is.null( bdgraph.obj3 ) )
		{ 
			est3_igraph <- graph.adjacency( as.matrix(est3), mode = "undirected", diag = FALSE ) 
			plot.igraph( est3_igraph, layout = layout.circle, main = colnames[4], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		par(op)
   }
   
	return( round( result, 3 ) )
}
    
# To compare the result
compute_measures = function( G, est_G ) 
{
	upper_G     = G[     upper.tri( G     ) ]
	upper_est_G = est_G[ upper.tri( est_G ) ]
		
	tp = sum( ( upper_G != 0 ) * ( upper_est_G != 0 ) ) 
	tn = sum( ( upper_G == 0 ) * ( upper_est_G == 0 ) )
	fp = sum( ( upper_G == 0 ) * ( upper_est_G != 0 ) ) 
	fn = sum( ( upper_G != 0 ) * ( upper_est_G == 0 ) )
		
	# harmonic mean of precision and recall, called F-measure or balanced F-score
	F1score = ( 2 * tp ) / ( 2 * tp + fp + fn )

	specificity  = tn / ( tn + fp )
	sensitivity  = tp / ( tp + fn )
	# Matthews Correlation Coefficients (MCC)
	mcc          = ( ( tp * tn ) - ( fp * fn ) ) / ( sqrt( ( tp + fp ) * ( tp + fn ) ) * sqrt( ( tn + fp ) * ( tn + fn ) ) )
	
	return( c( tp, tn, fp, fn, F1score, specificity, sensitivity, mcc ) )
}
   
