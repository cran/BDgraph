\name{adj2link}
\alias{adj2link}

\title{ Extract links from an adjacency matrix }

\description{ Extract links from an adjacency matrix or an object of calsses "\code{sim}" from function \code{\link{bdgraph.sim}} and
 \code{"graph"} from function \code{\link{graph.sim}}.}

\usage{ adj2link( adj ) }

\arguments{
  \item{adj}{     
  adjacency matrix corresponding to a graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
    }
}

\value{ matrix corresponding to the extracted links from graph structure.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03}

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
  
Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Pensar, J. et al (2017) Marginal pseudo-likelihood learning of discrete Markov network structures, \emph{Bayesian Analysis}, 12(4):1195-215, \doi{10.1214/16-BA1032}
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{link2adj}}, \code{\link{graph.sim}} }

\examples{
# Generating a 'random' graph 
adj <- graph.sim( p = 6, vis = TRUE )

adj2link( adj )
}

\keyword{graphs}
