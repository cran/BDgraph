## ----opts, echo = FALSE, message = FALSE, warning = FALSE---------------------
knitr::opts_chunk $ set( collapse = TRUE, comment = " ", fig.width = 7, fig.height = 7, fig.align = "center" )

## ----eval = FALSE-------------------------------------------------------------
#  install.packages( "BDgraph" )
#  
#  library( BDgraph )

## ----pressure, echo = FALSE, out.width = '85%'--------------------------------
knitr::include_graphics( "Figure_1.png" )

## ----eval = FALSE-------------------------------------------------------------
#  bdgraph( data, n = NULL, method = "ggm", algorithm = "bdmcmc", iter = 5000,
#           burnin = iter / 2, not.cont = NULL, g.prior = 0.5, df.prior = 3,
#           g.start = "empty", jump = NULL, save = FALSE,
#           cores = NULL, threshold = 1e-8, verbose = TRUE )

## ----eval = FALSE-------------------------------------------------------------
#  plinks( bdgraph.obj, round = 2, burnin = NULL )

## ----eval = FALSE-------------------------------------------------------------
#  select( bdgraph.obj, cut = NULL, vis = FALSE )

## ----eval = FALSE-------------------------------------------------------------
#  plotcoda( bdgraph.obj, thin = NULL, control = TRUE, main = NULL,
#            verbose = TRUE, ... )

## ----eval = FALSE-------------------------------------------------------------
#  traceplot( bdgraph.obj, acf = FALSE, pacf = FALSE, main = NULL, ... )

## ----eval = FALSE-------------------------------------------------------------
#  compare( pred, actual, main = NULL, vis = FALSE )

## ----eval = FALSE-------------------------------------------------------------
#  plotroc = function( pred, actual, cut = 20, smooth = FALSE, ... )

## ----eval = FALSE-------------------------------------------------------------
#  bdgraph.sim( p = 10, graph = "random", n = 0, type = "Gaussian", prob = 0.2,
#               size = NULL, mean = 0, class = NULL, cut = 4, b = 3,
#               D = diag( p ), K = NULL, sigma = NULL,
#               q = exp(-1), beta = 1, vis = FALSE, rewire = 0.05,
#               range.mu = c( 3, 5 ), range.dispersion = c( 0.01, 0.1 ) )

## ----eval = FALSE-------------------------------------------------------------
#  graph.sim( p = 10, graph = "random", prob = 0.2, size = NULL, class = NULL,
#             vis = FALSE, rewire = 0.05 )

## -----------------------------------------------------------------------------
library( BDgraph )

set.seed( 5 )

data.sim <- bdgraph.sim( n = 60, p = 8, graph = "scale-free", type = "Gaussian" )
round( head( data.sim $ data, 4 ), 2 ) 

## ----eval = TRUE--------------------------------------------------------------
sample.bdmcmc <- bdgraph( data = data.sim, method = "ggm", algorithm = "bdmcmc", 
                          iter = 5000, save = TRUE, verbose = FALSE )

## -----------------------------------------------------------------------------
summary( sample.bdmcmc )

## -----------------------------------------------------------------------------
sample.rjmcmc <- bdgraph( data = data.sim, method = "ggm", algorithm = "rjmcmc", 
                          iter = 5000, save = TRUE, verbose = FALSE )

## ----eval = FALSE-------------------------------------------------------------
#  plotroc( list( sample.bdmcmc, sample.rjmcmc ), data.sim, smooth = TRUE,
#           labels = c( "BDMCMC", "RJMCMC" ), color = c( "blue", "red" ) )

## -----------------------------------------------------------------------------
compare( list( sample.bdmcmc, sample.rjmcmc ), data.sim, 
          main = c( "True graph", "BDMCMC", "RJMCMC" ), vis = TRUE )

## -----------------------------------------------------------------------------
plotcoda( sample.bdmcmc, verbose = FALSE )
plotcoda( sample.rjmcmc, verbose = FALSE )

