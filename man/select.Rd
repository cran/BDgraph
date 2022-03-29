\name{select}
\alias{select}

\title{ Graph selection }

\description{
Provides the selected graph which, based on input, could be a graph with links for which their estimated posterior probabilities are greater than 0.5 (default)
or a graph with the highest posterior probability; see examples.
%For more specific selection of graphs consult the 'prob' function.  
}

\usage{ select( bdgraph.obj, cut = NULL, vis = FALSE ) }

\arguments{
  \item{bdgraph.obj}{ 
  matrix in which each element response to the weight of the links.
  It can be an object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}.
  It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
					}
					
  \item{cut}{ threshold for including the links in the selected graph based on the estimated posterior probabilities of the links; see the examples. }  
  \item{vis}{ visualize the selected graph structure. }
}

\value{	An adjacency matrix corresponding to the selected graph.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
  
Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
  
bdgraph.obj <- bdgraph( data = data.sim )
   
select( bdgraph.obj )
  
bdgraph.obj <- bdgraph( data = data.sim, save = TRUE )
  
select( bdgraph.obj )
  
select( bdgraph.obj, cut = 0.5, vis = TRUE )
}
}

\keyword{structure learning}
