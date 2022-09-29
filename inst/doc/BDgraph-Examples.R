## ----opts, echo = FALSE, message = FALSE, warning = FALSE---------------------
knitr::opts_chunk $ set( collapse = TRUE, comment = " ", fig.width = 7, fig.height = 7, fig.align = "center" )

## ----eval = FALSE-------------------------------------------------------------
#  install.packages( "BDgraph" )

## ----loadpkg, message = FALSE, warning = FALSE--------------------------------
library( BDgraph )

library( pROC )    
library( ggplot2 )  

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
compare( data.sim, list( bdgraph.obj, bdgraph.mpl.obj ), 
         main = c( "Target", "BDgraph", "BDgraph.mpl" ), vis = TRUE )

## ----fig.align = 'center'-----------------------------------------------------
roc.bdgraph     = BDgraph::roc( pred = bdgraph.obj,     actual = data.sim )
roc.bdgraph.mpl = BDgraph::roc( pred = bdgraph.mpl.obj, actual = data.sim )

pROC::ggroc( list( BDgraph = roc.bdgraph, BDgraph.mpl = roc.bdgraph.mpl ), size = 0.8 ) + 
    theme_minimal() + ggtitle( "ROC plots with AUC" ) +
  scale_color_manual( values = c( "red", "blue" ), 
    labels = c( paste( "AUC=", round( auc( roc.bdgraph ), 3 ), "; BDgraph; " ),
                paste( "AUC=", round( auc( roc.bdgraph.mpl ), 3 ), "; BDgraph.mpl" ) ) ) +
  theme( legend.title = element_blank() ) +
  theme( legend.position = c( .7, .3 ), text = element_text( size = 17 ) ) + 
    geom_segment( aes( x = 1, xend = 0, y = 0, yend = 1 ), color = "grey", linetype = "dashed" )

## ----fig.align = 'center'-----------------------------------------------------
set.seed( 2 )

data.sim = bdgraph.sim( n = 300, p = 10, type = "mixed", graph = "random", vis = TRUE )

## -----------------------------------------------------------------------------
bdgraph.obj = bdgraph( data = data.sim, method = "gcgm", iter = 5000, verbose = FALSE )

## ----fig.align = 'center'-----------------------------------------------------
compare( data.sim, bdgraph.obj, main = c( "Target", "BDgraph" ), vis = TRUE )

