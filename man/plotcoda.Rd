\name{plotcoda}
\alias{plotcoda}

\title{ Convergence plot }

\description{
Displays a plot which provides the cumulative occupancy fractions of all possible links in the graph.
The plot can be used for monitoring the convergence of the BD-MCMC sampling algorithm.
}

\usage{ plotcoda( output, thin = NULL, main = NULL, links = TRUE, ... ) }

\arguments{
  \item{output}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  
  \item{thin}{ An option for getting fast result for a cumulative plot according to part of the iteration.}
    
  \item{main}{ Graphical parameter (see plot). }
 
  \item{links}{ Logical: if TRUE and number of variables (p) is more that 15, you will see a plot result for 100 selected links. The default is TRUE. } 
  
  \item{\dots}{ System reserved (no specific usage). }  
}

\details{
Note that a spending time for this function depends on the number of nodes. 
It should be slow for the high-dimensional graphs. For fast faster, you can choose bigger value than other alternatives for the 'thin' option.
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108v2} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Gaussian Copula Graphical Models, \emph{Arxiv preprint arXiv:1501.04849v2} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )

output <- bdgraph( data = data.sim, iter = 10000, burnin = 0  )

plotcoda( output )
}
}
