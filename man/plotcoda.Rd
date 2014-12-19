\name{plotcoda}
\alias{plotcoda}

\title{ Convergence plot }

\description{
This function according to "bdgraph" output gives us a plot that can be used to check the convergence of the BD-MCMC algorithm.
}

\usage{ plotcoda( output, thin = NULL, main = NULL, ... ) }

\arguments{
  \item{output}{ An object with S3 class "bdgraph" }
  
  \item{thin}{ An option for getting fast result for a cumulative plot according to part of the iteration.}
    
  \item{main}{ Graphical parameter (see plot). }
  
  \item{\dots}{ System reserved (no specific usage). }  
}
\details{
Note that a spending time for this function depends on the graph. It should be slow for the 
high-dimensional graphs. To make it faster you can choose bigger value for 'thin' option.
}
\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating synthetic multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )

output <- bdgraph( data = data.sim )

plotcoda( output )
}
}
