\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}
\title{
Bayesian model selection in Gaussian graphical models based on BDMCMC algorithm
}
\description{
The R package BDgraph is a statistical tool for Bayesian model determination in 
undirected Gaussian graphical models based on birth-death MCMC methodology.      
}
\details{
  This package provides an implementation of the procedures described in Mohammadi and Wit (2012). 
The main function is 'bdmcmc' which is birth-death MCMC algorithm for Bayesian model selection 
in Gaussian graphical models.

  Functions:
  \preformatted{
bdmcmc           BDMCMC algorithm for undirected Gaussain graphical model determination
bdmcmc.high      BDMCMC algorithm for high-dimensional graphs
bdmcmc.low       BDMCMC algorithm for low-dimensional graphs
I.g              Computing normalizing constant of G-Wishart distribution
phat             Posterior edge inclusion probabilities 
plotConvergency  Cumulative occupancy fractions for checking the convergency
plotLinks        Plot of posterior distribution for graphs according to number of their links 
prob.allg        Posterior probability of all possible graphs
prob.g           Posterior probability for one special graph
sample.gwishart  Sampling from G-Wishart distribution
select.g         Selecting the best graphical models based on BDMCMC algorithm
  }
}
\author{
Abdolreza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, arXiv:1210.5371. \url{http://arxiv.org/abs/1210.5371}

Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the 
marginal likelihood in nondecomposable gaussian graphical models. Biometrika 
92(2), 317-335.
}
\keyword{ package }
