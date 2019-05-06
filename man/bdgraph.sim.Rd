\name{bdgraph.sim}
\alias{bdgraph.sim}

\title{ Graph data simulation }

\description{
Simulating multivariate distributions with different types of underlying graph structures, including 
\code{"random"}, \code{"cluster"}, \code{"scale-free"}, \code{"lattice"}, \code{"hub"}, \code{"star"}, \code{"circle"}, \code{"AR(1)"}, and \code{"AR(2)"}.
Based on the underling graph structure, it generates four different types of datasets, including \emph{multivariate Gaussian}, \emph{non-Gaussian}, \emph{discrete}, or \emph{mixed} data. 
This function can be used also for only simulating graphs by option \code{n=0}, as a default. 
}

\usage{
bdgraph.sim( p = 10, graph = "random", n = 0, type = "Gaussian", prob = 0.2, 
             size = NULL, mean = 0, class = NULL, cut = 4, b = 3,
             D = diag( p ), K = NULL, sigma = NULL, vis = FALSE )
}

\arguments{
  \item{p}{The number of variables (nodes).}
  \item{graph}{The graph structure with options 
                \code{"random"}, \code{"cluster"}, \code{"scale-free"}, \code{"lattice"}, \code{"hub"}, \code{"star"}, \code{"circle"}, 
                 \code{"AR(1)"}, and \code{"AR(2)"}.
     It also could be an adjacency matrix corresponding to a graph structure (an upper triangular matrix in which 
     \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}). 
    }
  \item{n}{The number of samples required. Note that for the case \code{n = 0}, only graph is generated. }
  \item{type}{Type of data with four options \code{"Gaussian"} (default), \code{"non-Gaussian"}, \code{"discrete"}, \code{"mixed"}, and \code{"binary"}.
	  For option \code{"Gaussian"}, data are generated from multivariate normal distribution.
	  For option \code{"non-Gaussian"}, data are transfered multivariate normal distribution to continuous multivariate non-Gaussian distribution.
	  For option \code{"discrete"}, data are transfered from multivariate normal distribution to discrete multivariat distribution.
	  For option \code{"mixed"}, data are transfered from multivariate normal distribution to mixture of 'count', 'ordinal', 'non-Gaussian', 'binary' and 'Gaussian', respectively. 
	  For option \code{"binary"}, data are generated directly from the joint distribution, in this case \eqn{p} must be less than \eqn{17}. 
	}

  \item{prob}{ If \code{graph="random"}, it is the probability that a pair of nodes has a link.}
  \item{size}{The number of links in the true graph (graph size).}
  \item{mean}{A vector specifies the mean of the variables.}
  \item{class}{ If \code{graph="cluster"}, it is the number of classes. }
  \item{cut}{ If \code{type="discrete"}, it is the number of categories for simulating discrete data.}
 
  \item{b}{The degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}.}
  \item{D}{The positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. The default is an identity matrix.}

  \item{K}{ If \code{graph="fixed"}, it is a positive-definite symmetric matrix specifies as a true precision matrix. }
  \item{sigma}{ If \code{graph="fixed"}, it is a positive-definite symmetric matrix specifies as a true covariance matrix.}
 
  \item{vis}{Visualize the true graph structure.}
}

\value{
	An object with \code{S3} class \code{"sim"} is returned:
	\item{data}{Generated data as an (\eqn{n \times p}{n x p}) matrix.}
	\item{sigma}{The covariance matrix of the generated data.}
	\item{K}{The precision matrix of the generated data.}
	\item{G}{The adjacency matrix corresponding to the true graph structure.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Pensar, J. et al (2017) Marginal pseudo-likelihood learning of discrete Markov network structures, \emph{Bayesian Analysis}, 12(4):1195-215
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{graph.sim}}, \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( p = 10, n = 50, prob = 0.3, vis = TRUE )
print( data.sim )
     
# Generating multivariate normal data from a 'hub' graph
data.sim <- bdgraph.sim( p = 6, n = 3, graph = "hub", vis = FALSE )
round( data.sim $ data, 2 )
     
# Generating mixed data from a 'hub' graph 
data.sim <- bdgraph.sim( p = 8, n = 10, graph = "hub", type = "mixed" )
round( data.sim $ data, 2 )

# Generating only a 'scale-free' graph (with no data) 
graph.sim <- bdgraph.sim( p = 8, graph = "scale-free" )
plot( graph.sim )
graph.sim $ G
}
}
