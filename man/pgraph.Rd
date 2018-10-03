\name{pgraph}
\alias{pgraph}

\title{ Posterior probabilities of the graphs }

\description{
Provides the estimated posterior probabilities for the most likely graphs or a specific graph.
}

\usage{ pgraph( bdgraph.obj, number.g = 4, adj_g = NULL ) }

\arguments{
  \item{bdgraph.obj}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{number.g}{ The number of graphs with the highest posterior probabilities to be shown.
            This option is ignored if \code{'adj_g'} is specified. }
  \item{adj_g}{ An adjacency matrix corresponding to a graph structure. It is an upper triangular matrix in which 
            \eqn{a_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{a_{ij}=0}.
            It also can be an object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}. }
}

\value{
\item{ selected_g }{ The graphs with the highest posterior probabilities.}

\item{ prob_g }{ A vector of the posterior probabilities of the graphs corresponding to \code{'selected\_g'}.}
}

\references{
Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and Wit, E. C. (2017). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108v5} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Mohammadi, A., Massam H., and Letac, G. (2017). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph}} and \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 6, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, save.all = TRUE )
   
# Estimated posterior probability of the true graph
pgraph( bdgraph.obj, adj_g = data.sim )
   
# Estimated posterior probability of first and second graphs with highest probabilities
pgraph( bdgraph.obj, number.g = 2 )
}
}
