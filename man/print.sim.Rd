\name{print.sim}
\alias{print.sim}

\title{ Print function for \code{S3} class "\code{sim}" }

\description{ Prints the information about the type of data, the sample size, the graph type, the number of nodes, number of links and sparsity of the true graph. }

\usage{ \method{print}{sim}( x, ... ) }

\arguments{
  \item{x}{ object of \code{S3} class "\code{sim}", from function \code{\link{bdgraph.sim}}.  }
  \item{\dots}{ system reserved (no specific usage). }
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{graph.sim}}, \code{\link{bdgraph.sim}} }

\examples{
# Generating multivariate normal data from a 'random' graph 
data.sim <- bdgraph.sim( n = 20, p = 10, vis = TRUE )
   
print( data.sim )
}

\keyword{print}
