## ----opts, echo = FALSE, message = FALSE, warning = FALSE---------------------
knitr::opts_chunk $ set( collapse = TRUE, comment = " ", fig.width = 7, fig.height = 7, fig.align = "center" )

## ----eval = FALSE-------------------------------------------------------------
#  install.packages( "BDgraph" )

## ----loadpkg, message = FALSE, warning = FALSE--------------------------------
library( BDgraph )

## ----fig.align = 'center'-----------------------------------------------------
set.seed( 20 )

data.sim = bdgraph.sim( n = 200, p = 15, graph = "scale-free", vis = TRUE )

## -----------------------------------------------------------------------------
bdgraph.obj = bdgraph( data = data.sim, method = "ggm", iter = 5000, verbose = FALSE )

## ----fig.align = 'center', fig.width = 3, fig.height = 3----------------------
conf.mat( actual = data.sim, pred = bdgraph.obj, cutoff = 0.5 )

conf.mat.plot( actual = data.sim, pred = bdgraph.obj, cutoff = 0.5 )

## ----fig.align = 'center'-----------------------------------------------------
compare( data.sim, bdgraph.obj, main = c( "Target", "BDgraph" ), vis = TRUE )

## ----fig.align = 'center', fig.width = 3, fig.height = 3----------------------
bdgraph.mpl.obj = bdgraph.mpl( data = data.sim, method = "ggm", iter = 5000, verbose = FALSE )

conf.mat( actual = data.sim, pred = bdgraph.mpl.obj )
conf.mat.plot( actual = data.sim, pred = bdgraph.mpl.obj )

## ----fig.align = 'center'-----------------------------------------------------
compare( list( bdgraph.obj, bdgraph.mpl.obj ), data.sim, 
         main = c( "Target", "BDgraph", "BDgraph.mpl" ), vis = TRUE )

## ----fig.align = 'center'-----------------------------------------------------
plotroc( list( bdgraph.obj, bdgraph.mpl.obj ), data.sim, cut = 200,
         labels = c( "BDgraph", "BDgraph.mpl" ), color = c( "blue", "red" ) )

## ----fig.align = 'center'-----------------------------------------------------
set.seed( 2 )

data.sim = bdgraph.sim( n = 300, p = 10, type = "mixed", graph = "random", vis = TRUE )

## -----------------------------------------------------------------------------
bdgraph.obj = bdgraph( data = data.sim, method = "gcgm", iter = 5000, verbose = FALSE )

## ----fig.align = 'center'-----------------------------------------------------
compare( bdgraph.obj, data.sim, main = c( "Target", "BDgraph" ), vis = TRUE )

## ----fig.align = 'center'-----------------------------------------------------
plotroc( bdgraph.obj, data.sim, labels = "BDgraph", color = "blue" )

