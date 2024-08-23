\name{posterior.predict}

\alias{posterior.predict}

\title{ Posterior Predictive Samples }

\description{
Provides samples from the posterior predictive distribution.  
}

\usage{
posterior.predict( object, iter = 1, ... )
}

\arguments{
  \item{object}{ object of \code{S3} class "\code{bdgraph}", from function \code{\link[BDgraph]{bdgraph}}. }
  \item{iter}{ number of predictions. }

  \item{\dots}{ additional parameters. }
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2023). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
  
Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Vogels, L., Mohammadi, R., Schoonhoven, M., and Birbil, S.I. (2023) Bayesian Structure Learning in Undirected Gaussian Graphical Models: Literature Review with Empirical Comparison, \emph{arXiv preprint}, \doi{10.48550/arXiv.2307.02603}

Mohammadi, R., Schoonhoven, M., Vogels, L., and Birbil, S.I. (2023) Large-scale Bayesian Structure Learning for Gaussian Graphical Models using Marginal Pseudo-likelihood, \emph{arXiv preprint}, \doi{10.48550/arXiv.2307.00127}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171} 
}

\value{ 
a \code{matrix} containing the predicted datasets, corresponding to the samples from the joint posterior disribtuion. 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link[BDgraph]{bdgraph}}, \code{\link[BDgraph]{bdgraph.mpl}}, \code{\link[BDgraph]{bdgraph.dw}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim )
   
posterior.predict( bdgraph.obj, iter = 20 )
}
}
