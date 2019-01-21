\name{plotcoda}
\alias{plotcoda}

\title{ Convergence plot }

\description{
Visualizes the cumulative occupancy fractions of all possible links in the graph.
It can be used for monitoring the convergence of the sampling algorithms, BDMCMC and RJMCMC.
}

\usage{ plotcoda( bdgraph.obj, thin = NULL, control = TRUE, main = NULL, ... ) }

\arguments{
  \item{bdgraph.obj}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph]{ssgraph}} of \code{R} package \code{\link[ssgraph]{ssgraph}}.
  }
  
  \item{thin}{ An option for getting fast result for a cumulative plot according to part of the iteration.}
    
  \item{control}{ Logical: if TRUE (default) and the number of nodes is greater than 15, then 100 links randomly is selected for visualization. } 
  
  \item{main}{ Graphical parameter (see plot). }
  
  \item{\dots}{ System reserved (no specific usage). }  
}

\details{
Note that a spending time for this function depends on the number of nodes. 
%It should be slow for high-dimensional graphs. 
For fast result, you can choose bigger value for the \code{'thin'} option.
}

\references{
Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and Wit, E. C. (2017). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108v5} 

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{traceplot}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )
  
bdgraph.obj <- bdgraph( data = data.sim, iter = 10000, burnin = 0 , save = TRUE )
   
plotcoda( bdgraph.obj )  
}
}
