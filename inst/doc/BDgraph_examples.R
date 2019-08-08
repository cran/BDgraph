## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#   Some examples for the BDgraph package
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

if( !require( "BDgraph" ) )
    install.packages( "BDgraph" )

## Loading the "BDgraph" package
library( "BDgraph" )
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  Example 1: Sampling from G-Wishart distribution
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

adj <- matrix( c( 0, 0, 1, 0, 0, 0, 1, 0, 0 ), 3, 3 )
adj   # adjacency matrix 
   
# Sampling from G-Wishart distribution with parameters b and D
sample <- rgwish( n = 1, adj = adj, b = 3, D = diag( 3 ) )
round( sample, 2 )
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  Example 2: An example on simulated data
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

# Generating multivariate Gaussian data based on a 'scale-free' graph structure
data.sim <- bdgraph.sim( n = 60, p = 8, graph = "scale-free", type = "Gaussian" )  
round( head( data.sim $ data, 4 ), 2 )
    
# Running the BDMCMC sampling algorithm for Gaussian graphical models (ggm) 
sample.bdmcmc <- bdgraph( data = data.sim, method = "ggm", algorithm = "bdmcmc", 
                          iter = 5000, save = TRUE )
   
summary( sample.bdmcmc ) 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
   
# Running the RJMCMC sampling algorithm for Gaussian graphical models (ggm) 
sample.rjmcmc <- bdgraph( data = data.sim, method = "ggm", algorithm = "rjmcmc", 
                          iter = 5000, save = TRUE )
   
plotroc( data.sim, sample.bdmcmc, sample.rjmcmc, smooth = TRUE, label = FALSE )
legend( "bottomright", c( "BDMCMC", "RJMCMC" ), lty = 1:2, col = c( 1, 4 ), lwd = c( 2, 2 ), cex = 1.5 )
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
   
# Comparing the performance of the BDMCMC and RJMCMC algorithms
compare( data.sim, sample.bdmcmc, sample.rjmcmc, main = c( "Target", "BDMCMC", "RJMCMC" ) )
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  Example 3: Application to labor force survey data
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

data( "surveyData", package = "BDgraph" )     # Loading the labor force survey dataset
head( surveyData, 5 )

## The data is stored in an integer matrix with the ordinal variables
## already having suitable scores assigned. Alternatively also a
## data.frame where ordinal variables are coded as factors could be
## bused and might be preferable.

# Running the BDMCMC sampling algorithm for Gaussian copula graphical models (gcgm) 
sample.bdmcmc <- bdgraph( data = surveyData, method = "gcgm", iter = 10000, burnin = 7000 )      
summary(sample.bdmcmc)
      
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  Example 4: Application to human gene expression
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

data( "geneExpression", package = "BDgraph" )  # Loading the human gene expression dataset
dim( geneExpression )
    
# Running the BDMCMC sampling algorithm for Gaussian copula graphical models (gcgm) 
sample.bdmcmc <- bdgraph( data = geneExpression, method = "gcgm", g.prior = 0.1, 
                          iter = 10000, burnin = 7000 )      
   
## - - - - Figure 8 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
if( !require( "Matrix" ) ) install.packages( "Matrix" )
library( "Matrix" ) 
p_links <- Matrix( plinks( sample.bdmcmc ) ) 

image( p_links + t( p_links ), useAbs = F, col.regions = grey( seq( 1, 0, length = 256 ) ), 
       xlab = NULL, ylab = NULL, sub = NULL, lwd = 0.01, 
       main = "Posterior Probabilities of all Links" )
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
    
selected_graph <- select( sample.bdmcmc, cut = 0.5 )
      
colnames( selected_graph ) <- substr( colnames( geneExpression ), 1, 7 ) 
    
# for removing nodes without links
nodes_size      <- apply( selected_graph + t( selected_graph ), 2, sum )
nodes_size_zero <- which( nodes_size == 0 )
selected_graph  <- selected_graph[ -nodes_size_zero, -nodes_size_zero ] 

plot.graph( selected_graph, layout = igraph::layout.fruchterman.reingold, edge.color = "gray20", 
            vertex.color = "red", vertex.size = 2, vertex.label.font = 1, 
            vertex.label.color = "gray50", vertex.label.dist = 0.15, 
            margin = -0.05 , edge.arrow.width = 5 )

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
