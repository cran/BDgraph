\name{bdgraph.sim}
\alias{bdgraph.sim}

\title{ Synthetic graph data generator }

\description{
This function implements a synthetic graph data generation for multivariate distributions with different types of underling graph structures, including 
\code{"random"}, \code{"cluster"}, \code{"scale-free"}, \code{"hub"}, \code{"fixed"}, and \code{"circle"}. 
Based on the underling graph structure, it generates four different types of datasets, including \emph{multivariate Gaussian}, \emph{non-Gaussian}, \emph{discrete}, or \emph{mixed} data.  
}

\usage{
bdgraph.sim( n = 2, p = 10, graph = "random", size = NULL, prob = 0.2, 
             class = NULL, type = "Gaussian", cut = 4, b = 3, D = diag(p), 
             K = NULL, sigma = NULL, mean = 0, vis = FALSE )
}

\arguments{
  \item{n}{The number of samples required. The default value is 2.}
  \item{p}{The number of variables (nodes). The default value is 10.}
  \item{graph}{The graph structure with option \code{"random"} (as a default), \code{"cluster"}, \code{"scale-free"}, \code{"hub"}, \code{"fixed"}, and \code{"circle"}. 
     It also could be adjacency matrix corresponding to a graph structure (an upper triangular matrix in which 
     \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}). }
  \item{size}{The number of links in the true graph (graph size).}
  \item{prob}{For "random" graph, it is the probability that a pair of nodes has a link. The default value is 0.2.}
  \item{class}{Number of classes for option \code{"cluster"}.}
  
  \item{type}{Type of data with four options \code{"Gaussian"} (as a default), \code{"non-Gaussian"}, \code{"discrete"}, and \code{"mixed"}.
  For option \code{"Gaussian"}, data are generated from multivariate normal distribution.
  For option \code{"non-Gaussian"}, data are transfered multivariate normal distribution to continuous multivariate non-Gaussian distribution.
  For option \code{"discrete"}, data are transfered from multivariate normal distribution to discrete multivariat distribution.
  For option \code{"mixed"}, data are transfered from multivariate normal distribution to mixture of 'count', 'ordinal', 'non-Gaussian', 'binary' and 'Gaussian', respectively. 
  }
  
  \item{cut}{Number of categories for simulating discrete data (\code{type = "discrete"}). The default value is 4.}

  \item{b}{The degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. The default is 3.}
  \item{D}{The positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. The default is an identity matrix.}

  \item{K}{A positive-definite symmetric matrix specifying the precision matrix. It is for option \code{graph = "fixed"}. }
  \item{sigma}{A positive-definite symmetric matrix specifying the covariance matrix. It is for option \code{graph = "fixed"}.}
  
  \item{mean}{A vector specifying the mean of the variables. The default value is a zero vector.}
  \item{vis}{Visualize the true graph pattern. The default value is FALSE.}
}

\value{
	An object with \code{S3} class \code{"simulate"} is returned:
	\item{data}{Generated data as a \eqn{n \times p}{n x p} matrix.}
	\item{sigma}{The covariance matrix of the generated data.}
	\item{K}{The precision matrix of the generated data.}
	\item{G}{The adjacency matrix corresponding to the true graph structure.}
}

\references{
Mohammadi, A. and E. Wit (2014). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis} (Forthcoming). \url{http://arxiv.org/abs/1210.5371v6}

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108}

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling 1}, 225-230
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 10, prob = 0.3, vis = TRUE )
print(data.sim)

# generating multivariate normal data from a 'hub' graph
data.sim <- bdgraph.sim( n = 3, p = 6, graph = "hub", vis = FALSE )
data.sim $ data

# generating mixed data from a 'hub' graph 
data.sim <- bdgraph.sim( n = 10, p = 8, graph = "hub", type = "mixed" )
data.sim $ data
}
}
