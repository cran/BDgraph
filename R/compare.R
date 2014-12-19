# To compare the result according to the true graph
compare = function ( G, est, est2 = NULL, colnames = NULL, vis = FALSE ) 
{
	if ( class(G)    == "simulate" ) G    <- G $ G 
	if ( class(G)    == "bdgraph" ){ es <- select(G) ; G <- est $ G ; est <- es }
	if ( class(est)  == "bdgraph" )  est  <- select( est ) 
	if ( class(est2) == "bdgraph" )  est2 <- select( est2 ) 
	if ( class(est)  == "select" )   est  <- est $ refit
	if ( class(est2) == "select" )   est2 <- est2 $ refit

	compare.true <- roc( G = G, est = G )
	compare.g    <- roc( G = G, est = est )

	if ( !is.null(est2) )
	{
		compare.g2  <- roc( G = G, est = est2 )
		compare.all <- cbind( compare.true, compare.g, compare.g2 )
		if ( is.null(colnames) ) colnames <- c("True graph", "estimate", "estimate2")
	} else {
		compare.all <- cbind( compare.true, compare.g )
		if ( is.null(colnames) ) colnames <- c( "True graph", "estimate" )
	}

   if ( vis == TRUE )
   {
		G   <- graph.adjacency( G,   mode = "undirected", diag = FALSE )
		est <- graph.adjacency( est, mode = "undirected", diag = FALSE )

		if ( is.null(est2) )
		{
			op <- par( mfrow = c(1, 2), pty = "s", omi = c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3) )
			plot.igraph( as.matrix(G), layout = layout.circle, main = colnames[1] )
			plot.igraph( as.matrix(est), layout = layout.circle, main = colnames[2] )
		} else {
			op   <- par( mfrow = c(2, 2), pty = "s", omi = c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3) )
			est2 <- graph.adjacency(as.matrix(est2), mode = "undirected", diag = FALSE)
			plot.igraph( G,    layout = layout.circle, main = colnames[1] )
			plot.igraph( est,  layout = layout.circle, main = colnames[2] )
			plot.igraph( est2, layout = layout.circle, main = colnames[3] )			
		}
		par(op)
   }
   
   colnames( compare.all ) <- colnames
   rownames( compare.all ) <- c("true positive", "true negative", "false positive", "false negative", 
                "true positive rate", "false positive rate", "accuracy", "balanced F-score", "positive predictive")
				
   return( compare.all )
}
   
