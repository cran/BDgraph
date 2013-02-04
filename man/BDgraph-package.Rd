\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}
\title{
Bayesian undirected graph estimation with BDMCMC algorithm
}
\description{
The R package BDgraph is a statistical tool for Bayesian model selection in 
undirected Gaussian graphical models based on birth-death MCMC methodology.      
}
\details{
  This package provides an implementation of the procedures described in Mohammadi and Wit (2012). 
The main function is 'bdgraph' which is birth-death MCMC algorithm for Bayesian model selection 
in Gaussian graphical models.

  Functions:
  \preformatted{
bdgraph          Graph estimation besed on BDMCMC algorithm 
bdgraph.sim      Data generator according to graph structure  
compare          Comparing the result according to the true graph
I.g              Computing normalizing constant of G-Wishart distribution
phat             Posterior edge inclusion probabilities
plot.bdgraph     Plot function for "bdgraph" output
plotcoda         Cumulative occupancy fractions for checking the convergency
densplot         Plot of posterior distribution for graphs according to their size
prob             Posterior probabilities for most highest graphs or specific graph
print.bdgraph    Print function for "bdgraph" output
summary.bdgraph  Summary function for "bdgraph" output
select           Selecting the best graphs
traceplot        Trace plot of graph size from "bdgraph" output
  }
}
\author{
Abdolreza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}

Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the 
marginal likelihood in nondecomposable gaussian graphical models. Biometrika 92(2), 317-335.

Wang, H. and S. Li (2012). Efficient Gaussian graphical model determination under 
G-Wishart prior distributions. Electronic Journal of Statistics
6, 168-198.
}
\keyword{ package }
