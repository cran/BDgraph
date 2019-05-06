\name{summary.bdgraph}
\alias{summary.bdgraph}

\title{ Summary function for \code{S3} class \code{"bdgraph"} }

\description{
Provides a summary of the results for function \code{\link{bdgraph}}.  
}

\usage{
\method{summary}{bdgraph}( object, round = 2, vis = TRUE, ... )
}

\arguments{
  \item{object}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{round}{ A value for rounding all probabilities to the specified number of decimal places. }
  \item{vis}{ Visualize the results. }
  \item{\dots}{System reserved (no specific usage).}
}

\value{
	\item{selected_g}{The adjacency matrix corresponding to the selected graph which has the highest posterior probability.}
	\item{p_links}{An upper triangular matrix corresponding to the posterior probabilities of all possible links.}
	\item{K_hat}{The estimated precision matrix.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 
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
