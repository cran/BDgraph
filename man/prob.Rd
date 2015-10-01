\name{prob}
\alias{prob}

\title{ Posterior probabilities of the graphs }

\description{
Provides the estimated posterior probabilities for the most likely graphs or a specific graph, \code{'G'}.
}

\usage{prob( x, number.g = 4, G = NULL )}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{number.g}{ The number of graphs with the highest posterior probabilities to be shown (default is 4).
            This option is ignored if \code{'G'} is specified. }
  \item{G}{ Adjacency matrix corresponding to a graph structure. It is an upper triangular matrix in which 
            \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}.
            It also can be an object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}. }
}

\value{
\item{ selected_G }{ The graphs with the highest posterior probabilities.}

\item{ prob_G }{ A vector of the posterior probabilities of the graphs corresponding to \code{'selected\_G'}.}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{ \code{\link{bdgraph}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 6, vis = TRUE )
   
output   <- bdgraph( data = data.sim, save.all = TRUE )
   
# Estimated posterior probability of the true graph
prob( output, G = data.sim )
   
# Estimated posterior probability of the first and second graphs with highest probabilities
prob( output, number.g = 2 )
}
}
