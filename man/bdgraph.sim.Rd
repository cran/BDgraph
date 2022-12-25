\name{bdgraph.sim}
\alias{bdgraph.sim}

\title{ Graph data simulation }

\description{
Simulating multivariate distributions with different types of underlying graph structures, including 
"\code{random}", "\code{cluster}", "\code{smallworld}", "\code{scale-free}", "\code{lattice}", "\code{hub}", "\code{star}", "\code{circle}", "\code{AR(1)}", and "\code{AR(2)}".
Based on the underlying graph structure, the function generates different types of \emph{multivariate} data, including "\emph{Gaussian}", "\emph{non-Gaussian}", "\emph{categorical}", "\emph{pois}" (Poisson), "\emph{nbinom}" (negative binomial), "\emph{dweibull}" (discrete Weibull), "\emph{binary}", "\emph{t}" (t-distribution), "\emph{alternative-t}", or "\emph{mixed}" data. 
This function can be used also for simulating only graphs by setting the option \code{n=0} (default). 
}

\usage{
bdgraph.sim( p = 10, graph = "random", n = 0, type = "Gaussian", prob = 0.2, 
             size = NULL, mean = 0, class = NULL, cut = 4, b = 3,
             D = diag( p ), K = NULL, sigma = NULL, 
             q = exp(-1), beta = 1, vis = FALSE, rewire = 0.05,
             range.mu = c( 3, 5 ), range.dispersion = c( 0.01, 0.1 ), nu = 1 )
}

\arguments{
  \item{p}{number of variables (nodes).}
  \item{graph}{graph structure with options 
                "\code{random}", "\code{cluster}", "\code{smallworld}", "\code{scale-free}", "\code{lattice}", "\code{hub}", "\code{star}", "\code{circle}", 
                 "\code{AR(1)}", and "\code{AR(2)}".
     It could also be an adjacency matrix corresponding to a graph structure (an upper triangular matrix in which 
     \eqn{g_{ij}=1} if there is a link between nodes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}). 
    }
  \item{n}{number of samples required. Note that for the case \code{n = 0}, only the graph is generated. }
  \item{type}{type of data with options "\code{Gaussian}" (default), "\code{non-Gaussian}", "\code{categorical}", "\code{pois}", "\code{nbinom}", "\code{dweibull}", "\code{binary}", "\code{mixed}", "\code{t}", and "\code{alternative-t}".    
	  For the option "\code{Gaussian}", data are generated from a multivariate normal distribution.
	  For the option "\code{non-Gaussian}", data are transfered from a multivariate normal distribution to a continuous multivariate non-Gaussian distribution via Exponential marginals.
	  For the option "\code{categorical}", data are transfered from a multivariate normal distribution to multivariate 'categorical' data.
      For the option "\code{pois}", data are transfered from a multivariate normal distribution to a multivariate Poisson distribution.	  
      For the option "\code{nbinom}", data are transfered from a multivariate normal distribution to a multivariate Negative Binomial distribution.	  
      For the option "\code{dweibull}", data are transfered from a multivariate normal distribution to a multivariate discrete Weibull distribution with parameters \code{q} and \code{beta}.	  
	  For the option "\code{binary}", data are generated directly from the joint distribution, in this case \eqn{p} must be less than \eqn{17}. 
	  For the option "\code{mixed}", data are transfered from a multivariate normal distribution to a mixture of 'categorical', 'non-Gaussian', 'binary' and 'Gaussian', respectively. 
	}

  \item{prob}{ if \code{graph} = "\code{random}", it is the probability that a pair of nodes has a link.}
  \item{size}{number of links in the true graph (graph size).}
  \item{mean}{vector specifying the mean of the variables.}
  \item{class}{ if \code{graph} = "\code{cluster}", it is the number of classes. }
  \item{cut}{ if \code{type} = "\code{categorical}", it is the number of categories for simulating 'categorical' data.}
 
  \item{b}{degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}.}
  \item{D}{positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. The default is an identity matrix.}

  \item{K}{     if \code{graph} = "\code{fixed}", it is a positive-definite symmetric matrix, corresponding to the true precision matrix. }
  \item{sigma}{ if \code{graph} = "\code{fixed}", it is a positive-definite symmetric matrix corresponding to the true covariance matrix.}

  \item{q, beta}{ if \code{type} = "\code{dweibull}", they are the parameters of the discrete Weibull distribution with density 
        \deqn{p( x, q, \beta ) = q^{x^{\beta}}-q^{(x+1)^{\beta}}, \quad \forall x = \{ 0, 1, 2, \ldots \}.} They can be given either as a vector of length p or as an  (\eqn{n \times p}{n x p}) matrix, e.g. if covariates are available and a regression model is used.}

  \item{vis}{visualize the true graph structure.}
  \item{rewire}{rewiring probability for smallworld network. Must be between 0 and 1.}
  \item{range.mu, range.dispersion}{if \code{type} = "\code{nbinom}", vector with two elements specifying the range of parameters for the Negative Binomial distribution.}
    \item{nu}{ if \code{type} = "\code{t}" or "\code{alternative-t}", it is the parameter of the t distribution with density.}

}

\value{
	An object with \code{S3} class "\code{sim}" is returned:
	\item{data}{generated data as an (\eqn{n \times p}{n x p}) matrix.}
	\item{sigma}{covariance matrix of the generated data.}
	\item{K}{precision matrix of the generated data.}
	\item{G}{adjacency matrix corresponding to the true graph structure.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}, Pariya Behrouzi, Veronica Vinciotti, Ernst Wit, and Alexander Christensen }

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

\keyword{ datagen }
\keyword{ simulation }
