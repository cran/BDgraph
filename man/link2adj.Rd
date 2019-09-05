\name{link2adj}
\alias{link2adj}

\title{ Extract links from an adjacency matrix }

\description{ Extract links from an adjacency matrix or an object of calsses \code{"sim"} from function \code{\link{bdgraph.sim}} and
 \code{"graph"} from function \code{\link{graph.sim}}.}

\usage{ link2adj( link, p = NULL ) }

\arguments{
  \item{link}{ An (\eqn{2 \times p}) \code{matrix} or a \code{data.frame} corresponding to the links from the graph structure. }
  \item{p}{The number of nodes of the graph.}
}

\value{ An adjacency matrix corresponding to the graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Pensar, J. et al (2017) Marginal pseudo-likelihood learning of discrete Markov network structures, \emph{Bayesian Analysis}, 12(4):1195-215
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{adj2link}}, \code{\link{graph.sim}} }

\examples{
# Generating a 'random' graph 
adj <- graph.sim( p = 6, vis = TRUE )

link <- adj2link( adj )

link2adj( link, p = 6 )
   
}

\keyword{graphs}
