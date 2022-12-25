\name{summary.bdgraph}

\alias{summary.bdgraph}

\title{ Summary function for \code{S3} class "\code{bdgraph}" }

\description{
Provides a summary of the results for function \code{\link{bdgraph}}.  
}

\usage{
\method{summary}{bdgraph}( object, round = 2, vis = TRUE, ... )
}

\arguments{
  \item{object}{ object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. }
  \item{round}{ value for rounding all probabilities to the specified number of decimal places. }
  \item{vis}{ visualize the results. }

  \item{\dots}{ additional plotting parameters for the case \code{vis = TRUE}. See \code{\link[BDgraph]{plot.graph}}. }
}

\value{
	\item{selected_g}{adjacency matrix corresponding to the selected graph which has the highest posterior probability.}
	\item{p_links}{upper triangular matrix corresponding to the posterior probabilities of all possible links.}
	\item{K_hat}{estimated precision matrix.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim )
   
summary( bdgraph.obj )
   
bdgraph.obj <- bdgraph( data = data.sim, save = TRUE )
   
summary( bdgraph.obj )
   
summary( bdgraph.obj, vis = FALSE )
}
}
