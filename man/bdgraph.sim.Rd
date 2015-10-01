\name{bdgraph.sim}
\alias{bdgraph.sim}

\title{ Synthetic graph data generator }

\description{
Implements a synthetic graph data generation for multivariate distributions with different types of underlying graph structures, including 
\code{"random"}, \code{"cluster"}, \code{"scale-free"}, \code{"hub"}, \code{"fixed"}, and \code{"circle"}. 
Based on the underling graph structure, it generates four different types of datasets, including \emph{multivariate Gaussian}, \emph{non-Gaussian}, \emph{discrete}, or \emph{mixed} data.  
}

\usage{
bdgraph.sim( n = 2, p = 10, type = "Gaussian", graph = "random", prob = 0.2, 
			 size = NULL, mean = 0, class = NULL, cut = 4, b = 3, 
			 D = diag(p), K = NULL, sigma = NULL, vis = FALSE )
}

\arguments{
  \item{n}{The number of samples required. The default value is 2.}
  \item{p}{The number of variables (nodes). The default value is 10.}
  \item{type}{Type of data with four options \code{"Gaussian"} (as a default), \code{"non-Gaussian"}, \code{"discrete"}, and \code{"mixed"}.
	  For option \code{"Gaussian"}, data are generated from multivariate normal distribution.
	  For option \code{"non-Gaussian"}, data are transfered multivariate normal distribution to continuous multivariate non-Gaussian distribution.
	  For option \code{"discrete"}, data are transfered from multivariate normal distribution to discrete multivariat distribution.
	  For option \code{"mixed"}, data are transfered from multivariate normal distribution to mixture of 'count', 'ordinal', 'non-Gaussian', 'binary' and 'Gaussian', respectively. 
	}

  \item{graph}{The graph structure with option \code{"random"} (default), \code{"cluster"}, \code{"scale-free"}, \code{"hub"}, \code{"fixed"}, and \code{"circle"}. 
     It also could be an adjacency matrix corresponding to a graph structure (an upper triangular matrix in which 
     \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}). 
    }
  \item{prob}{ If \code{graph="random"}, it is the probability that a pair of nodes has a link. The default value is 0.2.}
  \item{size}{The number of links in the true graph (graph size).}
  \item{mean}{A vector specifies the mean of the variables. The default value is a zero vector.}
  \item{class}{ If \code{graph="cluster"}, it is the number of classes. }
  \item{cut}{ If \code{type="discrete"}, it is the number of categories for simulating discrete data. The default value is 4.}
 
  \item{b}{The degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. The default is 3.}
  \item{D}{The positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. The default is an identity matrix.}

  \item{K}{ If \code{graph="fixed"}, it is a positive-definite symmetric matrix specifies as a true precision matrix. }
  \item{sigma}{ If \code{graph="fixed"}, it is a positive-definite symmetric matrix specifies as a true covariance matrix.}
 
  \item{vis}{Visualize the true graph structure. The default value is FALSE.}
}

\value{
	An object with \code{S3} class \code{"sim"} is returned:
	\item{data}{Generated data as a (\eqn{n \times p}{n x p}) matrix.}
	\item{sigma}{The covariance matrix of the generated data.}
	\item{K}{The precision matrix of the generated data.}
	\item{G}{The adjacency matrix corresponding to the true graph structure.}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 10, prob = 0.3, vis = TRUE )
print( data.sim )
   
# Generating multivariate normal data from a 'hub' graph
data.sim <- bdgraph.sim( n = 3, p = 6, graph = "hub", vis = FALSE )
round( data.sim $ data, 2 )
   
# Generating mixed data from a 'hub' graph 
data.sim <- bdgraph.sim( n = 10, p = 8, graph = "hub", type = "mixed" )
round( data.sim $ data, 2 )
}
}
