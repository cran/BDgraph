# To compare the result according to the true graph
compare = function ( G, est, est2 = NULL, est3 = NULL, colnames = NULL, vis = FALSE ) 
{
	if ( class(G)    == "sim" ) G    <- G $ G 
	if ( class(G)    == "bdgraph" ){ es <- select(G) ; G <- est $ G ; est <- es }
	if ( class(est)  == "bdgraph" )  est  <- select( est ) 
	if ( class(est2) == "bdgraph" )  est2 <- select( est2 ) 
	if ( class(est3) == "bdgraph" )  est3 <- select( est3 ) 
	if ( class(est)  == "select" )   est  <- est $ refit
	if ( class(est2) == "select" )   est2 <- est2 $ refit
	if ( class(est3) == "select" )   est3 <- est3 $ refit

	compare.true <- roc( G = G, est = G )
	compare.g    <- roc( G = G, est = est )

	if ( is.null(est2) & is.null(est3) )
	{
		compare.all <- cbind( compare.true, compare.g )
		if ( is.null(colnames) ) colnames <- c( "True graph", "estimate" )
	}

	if ( !is.null(est2) & is.null(est3) )
	{
		compare.g2  <- roc( G = G, est = est2 )
		compare.all <- cbind( compare.true, compare.g, compare.g2 )
		if ( is.null(colnames) ) colnames <- c( "True graph", "estimate", "estimate2" )
	} 
	
	if ( !is.null(est3) & !is.null(est3) )
	{
		compare.g2  <- roc( G = G, est = est2 )
		compare.g3  <- roc( G = G, est = est3 )
		compare.all <- cbind( compare.true, compare.g, compare.g2, compare.g3 )
		if ( is.null(colnames) ) colnames <- c( "True graph", "estimate", "estimate2", "estimate3" )
	} 
	
   if ( vis )
   {
		p = dim(G)[1]
		G   <- graph.adjacency( G,   mode = "undirected", diag = FALSE )
		est <- graph.adjacency( est, mode = "undirected", diag = FALSE )
		if ( p < 20 ) sizev = 15 else sizev = 2

		if ( is.null(est2) )
		{
			op <- par( mfrow = c(1, 2), pty = "s", omi = c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3) )
			plot.igraph( as.matrix(G), layout = layout.circle, main = colnames[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( as.matrix(est), layout = layout.circle, main = colnames[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		}
		 
		if ( !is.null(est2) & is.null(est3) )
		{
			op   <- par( mfrow = c(2, 2), pty = "s", omi = c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3) )
			est2 <- graph.adjacency( as.matrix(est2), mode = "undirected", diag = FALSE )
			plot.igraph( G,    layout = layout.circle, main = colnames[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est,  layout = layout.circle, main = colnames[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est2, layout = layout.circle, main = colnames[3], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		if ( !is.null(est2) & !is.null(est3) )
		{
			op   <- par( mfrow = c(2, 2), pty = "s", omi = c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3) )
			est2 <- graph.adjacency( as.matrix(est2), mode = "undirected", diag = FALSE ) 
			est3 <- graph.adjacency( as.matrix(est3), mode = "undirected", diag = FALSE ) 
			plot.igraph( G,    layout = layout.circle, main = colnames[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est,  layout = layout.circle, main = colnames[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est2, layout = layout.circle, main = colnames[3], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
			plot.igraph( est3, layout = layout.circle, main = colnames[4], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		par(op)
   }
   
   colnames( compare.all ) <- colnames
   rownames( compare.all ) <- c("true positive", "true negative", "false positive", "false negative", 
                "true positive rate", "false positive rate", "accuracy", "balanced F-score", "positive predictive")
				
   return( compare.all )
}
         
