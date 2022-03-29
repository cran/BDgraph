\name{bdgraph.npn}
\alias{bdgraph.npn}

\title{Nonparametric transfer}

\description{ Transfers non-Gaussian data to Gaussian.  }

\usage{ bdgraph.npn( data, npn = "shrinkage", npn.thresh = NULL ) }

\arguments{
  \item{data}{ (\eqn{n \times p}{n x p}) \code{matrix} or a \code{data.frame} corresponding to the data (\eqn{n} is the sample size and \eqn{p} is the number of variables). }
  \item{npn}{ 
  character with three options \code{"shrinkage"} (default), \code{"truncation"}, and \code{"skeptic"}. Option \code{"shrinkage"} is for the 
  shrunken transformation, option \code{"truncation"} is for the truncated transformation and option \code{"skeptic"} is for the 
  non-paranormal skeptic transformation. For more details see references.
  }
  \item{npn.thresh}{ truncation threshold; it is only for the truncated transformation (\code{npn=} \code{"truncation"}). 
  The default value is \eqn{1/(4n^{1/4} \sqrt{\pi \log(n)})}. }
}

\value{ (\eqn{n \times p}) matrix of transferred data, if \code{npn = "shrinkage" or "truncation"}, and 
  a non-paranormal correlation (\eqn{p \times p}) matrix, if \code{npn = "skeptic"}.  
  }

\references{
Liu, H., et al (2012). High Dimensional Semiparametric Gaussian Copula Graphical Models, \emph{Annals of Statistics}, 40(4):2293-2326

Zhao, T. and Liu, H. (2012). The \pkg{huge} Package for High-dimensional Undirected Graph Estimation in \code{R}, \emph{Journal of Machine Learning Research}, 13:1059-1062
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{bdgraph.sim}}, \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 6, p = 4, size = 4 )

data     <- ( data.sim $ data - 3 ) ^ 4
data
   
# Transfer the data by truncation 
bdgraph.npn( data, npn = "truncation" )
  
# Transfer the data by shrunken 
bdgraph.npn( data, npn = "shrunken" )
  
# Transfer the data by skeptic 
bdgraph.npn( data, npn = "skeptic" )
}
}

\keyword{data preprocessing}
