\name{plot.bdgraph}
\alias{plot.bdgraph}

\title{ Plot function for \code{S3} class \code{"bdgraph"} }

\description{ Visualizes structure of the selected graphs which could be a graph with links for which their estimated posterior probabilities are more that 0.5
or graph with the highest posterior probability. }

\usage{ \method{plot}{bdgraph}( x, cut = 0.5, number.g = 1, layout = layout.circle, ... ) }

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  
  \item{cut}{ Threshold for including the links in the selected graph based on the estimated posterior probabilities of the links; See the examples. } 

  \item{number.g}{ The number of graphs with the highest probabilities (default is 1).}
  
  \item{layout}{The vertex placement algorithm which is according to \pkg{igraph} package. The default is \code{"layout.circle"}.}
  
  \item{\dots}{System reserved (no specific usage).}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
   
output   <- bdgraph( data = data.sim )
   
plot( output )
   
output   <- bdgraph( data = data.sim, save.all = TRUE )
   
plot( output, number.g = 4 )
  
plot( output, cut = 0.4 )
}
}
