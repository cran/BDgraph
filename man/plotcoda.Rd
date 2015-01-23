\name{plotcoda}
\alias{plotcoda}

\title{ Convergence plot }

\description{
Displays a plot which provides the cumulative occupancy fractions of all possible links in the graph.
The plot can be used for monitoring the convergence of the BD-MCMC sampling algorithm.
}

\usage{ plotcoda( output, thin = NULL, main = NULL, ... ) }

\arguments{
  \item{output}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  
  \item{thin}{ An option for getting fast result for a cumulative plot according to part of the iteration.}
    
  \item{main}{ Graphical parameter (see plot). }
  
  \item{\dots}{ System reserved (no specific usage). }  
}
\details{
Note that a spending time for this function depends on the graph. It should be slow for the 
high-dimensional graphs. To make it faster, you can choose bigger value for 'thin' option.
}
\references{
Mohammadi, A. and E. Wit (2014). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis} (Forthcoming). \url{http://arxiv.org/abs/1210.5371v6}

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108}

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling 1}, 225-230
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )

output <- bdgraph( data = data.sim )

plotcoda( output )
}
}
