\name{print.bdgraph}
\alias{print.bdgraph}

\title{Print function for \code{S3} class \code{"bdgraph"} }

\description{
Print the information about the best graph which is the graph with the highest posterior probability.
It provides adjacency matrix, size and posterior probability of the best graph.
}

\usage{\method{print}{bdgraph}( x, round = 3, Khat = FALSE, phat = FALSE, ... )}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{round}{A value to round the probabilities to the specified number of decimal places (default is 3).}
  
  \item{Khat}{Logical: if TRUE you will see the estimation of the precision matrix. The default is FALSE.}
  
  \item{phat}{Logical: if TRUE you will see the posterior link probabilities. The default is FALSE.}
 
  \item{\dots}{System reserved (no specific usage). }
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
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, vis = TRUE )

output   <- bdgraph( data = data.sim )

print( output )

print( output, Khat = TRUE, phat = TRUE )
}
}
