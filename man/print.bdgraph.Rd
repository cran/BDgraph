\name{print.bdgraph}
\alias{print.bdgraph}

\title{Print function for "bdgraph" output}

\description{According to output of the BD-MCMC algorithm, this function gives us the print of the result.  }

\usage{\method{print}{bdgraph}( x, round = 3, Khat = FALSE, phat = FALSE, ... )}

\arguments{
  \item{x}{ An object with S3 class "bdgraph". }
  \item{round}{A number for rounding all probabilities to the specified number of decimal places (default 3).}
  
  \item{Khat}{Logical: if TRUE you will see the estimation of precision matrix. The default is FALSE.}
  
  \item{phat}{Logical: if TRUE you will see the posterior edge inclusion probabilities. The default is FALSE.}
 
  \item{\dots}{System reserved (no specific usage). }
  
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
  # # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 50, p = 6, vis = TRUE )
  
  output   <- bdgraph( data = data.sim, iter = 2000 )
  
  print(output)
  
  print( output, Khat = TRUE, phat = TRUE )
  }
}
