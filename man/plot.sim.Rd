\name{plot.sim}
\alias{plot.sim}

\title{ Plot function for \code{S3} class \code{"sim"} }

\description{ Visualizes structure of the true graph. }
\usage{
\method{plot}{sim}( x, main = NULL, layout = layout.circle, ... )
}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}.  }
  
  \item{main}{ Graphical parameter (see plot).}
  
  \item{layout}{The vertex placement algorithm which is according to \pkg{igraph} package. }
  
  \item{\dots}{System reserved (no specific usage).}
}
\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108v2} 

Dobra, A. and A. Mohammadi (2017). Loglinear Model Selection and Human Mobility, \emph{arXiv preprint arXiv:1711.02623}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 

Mohammadi, A., Massam H., and G. Letac (2017). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416} 

Mohammadi, A. and A. Dobra (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# Generating synthetic multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 10, p = 15 )
  
plot( data.sim )
}
}
