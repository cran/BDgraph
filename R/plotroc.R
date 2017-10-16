# To plot ROC curve
plotroc = function( sim.obj, bdgraph.obj, bdgraph.obj2 = NULL, bdgraph.obj3 = NULL, 
                    bdgraph.obj4 = NULL, cut = 20, smooth = FALSE, label = TRUE, main = "ROC Curve" )
{
    if ( class( sim.obj ) == "sim" ) G = as.matrix( sim.obj $ G ) else G = as.matrix( sim.obj )
	G[ lower.tri( G, diag = TRUE ) ] = 0
	
    output_tp_fp = compute_tp_fp( G = G, bdgraph.obj = bdgraph.obj, cut = cut, smooth = smooth )
    fp           = output_tp_fp $ fp
    tp           = output_tp_fp $ tp
 	
	# par( mar = c( 3.8, 4.2, 1.8, 1 ) )
    plot( NA, type = "l", col = "black", cex.lab = 1.3, cex.main = 2, cex.axis = 1.2,
          main = main, xlab = "False Postive Rate", ylab = "True Postive Rate", 
          ylim = c( 0, 1 ), xlim = c( 0, 1 ) )
    points( x = fp, y = tp, type = "l", col = "black", lty = 1, lw = 2 )
  
    if( !is.null( bdgraph.obj2 ) )
    {
        output_tp_fp = compute_tp_fp( G = G, bdgraph.obj = bdgraph.obj2, cut = cut, smooth = smooth )
		fp_2         = output_tp_fp $ fp
		tp_2         = output_tp_fp $ tp
	
        points( x = fp_2, y = tp_2, type = "l", col = "blue", lty = 2, lw = 2 )
    }
    
    if( !is.null( bdgraph.obj3 ) )
    {   
        output_tp_fp = compute_tp_fp( G = G, bdgraph.obj = bdgraph.obj3, cut = cut, smooth = smooth )
		fp_3         = output_tp_fp $ fp
		tp_3         = output_tp_fp $ tp
   		
        points( x = fp_3, y = tp_3, type = "l", col = "green", lty = 3, lw = 2 )
    }
    
    if( !is.null( bdgraph.obj4 ) )
    {   
        output_tp_fp = compute_tp_fp( G = G, bdgraph.obj = bdgraph.obj4, cut = cut, smooth = smooth )
		fp_4         = output_tp_fp $ fp
		tp_4         = output_tp_fp $ tp
   		
        points( x = fp_4, y = tp_4, type = "l", col = "red", lty = 4, lw = 2 )
    }
    
	if ( label )
	{ 
		if( !is.null( bdgraph.obj2 ) && is.null( bdgraph.obj3 ) )  legend( "bottomright", c( "bdgraph.obj", "bdgraph.obj2" ), lty = 1:2, col = c( "black", "blue" ), lwd = c( 2, 2 ), cex = 1.5 )
		if( !is.null( bdgraph.obj3 ) && is.null( bdgraph.obj4 ) )  legend( "bottomright", c( "bdgraph.obj", "bdgraph.obj2", "bdgraph.obj3" ), lty = 1:3, col = c( "black", "blue", "green" ), lwd = c( 2, 2 ), cex = 1.5 )
		if( !is.null( bdgraph.obj4 ) ) legend( "bottomright", c( "bdgraph.obj", "bdgraph.obj2", "bdgraph.obj3", "bdgraph.obj4" ), lty = 1:4, col = c( "black", "blue", "green", "red" ), lwd = c( 2, 2 ), cex = 1.5 )
	}   
}
       
# function for ROC plot
compute_tp_fp = function( G, bdgraph.obj, cut, smooth )
{
	p           = nrow( G )
	upper_G     = G[ upper.tri( G ) ]
	sum_edges   = sum( upper_G )
	sum_no_dges = p * ( p - 1 ) / 2 - sum_edges

    if( class( bdgraph.obj ) != "huge" )
    {
		if( class( bdgraph.obj ) == "bdgraph" )
		{
			p_links = bdgraph.obj $ p_links
			if( is.null( p_links ) ) p_links = plinks( bdgraph.obj, round = 10 )
			p_links = as.matrix( p_links )
		}else{
			p_links = as.matrix( bdgraph.obj )
		}
		
		tp = c( 1, rep( 0, cut ) )
		fp = tp

		cut_points = ( 0 : cut ) / cut
		
		for( i in 2 : cut )
		{
			# checking for cut pints
			est_G = matrix( 0, p, p )
			est_G[ p_links > cut_points[i] ] = 1
			upper_est_G = est_G[ upper.tri( est_G ) ]

			tp[i] = sum( ( upper_G != 0 ) * ( upper_est_G != 0 ) ) / sum_edges
			fp[i] = sum( ( upper_G == 0 ) * ( upper_est_G != 0 ) ) / sum_no_dges
		}
	}else{
		path = bdgraph.obj $ path
		tp   = numeric( length( path ) )
		fp   = tp

		for( i in 1 : length( path ) )
		{
			est_G       = as.matrix( path[[ i ]] )
			upper_est_G = est_G[ upper.tri( est_G ) ]
			
			tp[i] = sum( ( upper_G != 0 ) * ( upper_est_G != 0 ) ) / sum_edges
			fp[i] = sum( ( upper_G == 0 ) * ( upper_est_G != 0 ) ) / sum_no_dges
		}
		
		tp = c( tp, 1 )
		fp = c( fp, 1 )
	}
	
	if ( smooth == TRUE )
	{
		fit = smooth.spline( x = fp, y = tp )
		fp  = c( 0, fit $ x )
		tp  = c( 0, fit $ y )
	}	
	
	return( list( tp = tp, fp = fp ) )
}
    
