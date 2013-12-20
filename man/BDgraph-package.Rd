\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}

\title{Graph selection based on birth-death MCMC}

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
bdgraph          Graph selection based on birth-death MCMC algorithm
bdgraph.sim      Data generator according to graph structure 
bdgraph.npt      Nonparametric transfer 
compare          Comparing the result 
I.g              Normalizing constant of G-Wishart distribution
phat             Posterior edge inclusion probabilities
plot.bdgraph     Plot function for "bdgraph" output
plotcoda         Convergency plots
print.bdgraph    Print function for "bdgraph" output
prob             Posterior probabilities of the graphs
rGWishart        Sampling from G-Wishart distribution
select           Selecting the best graphs
CellSignal       A flow cytometry dataset 
summary.bdgraph  Summary function for "bdgraph" output
traceplot        Trace plot of graph size from "bdgraph" output
  }
}

\author{Abdolreza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit}

\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}

Atay-Kayis, A. and H. Massam (2005). A Monte Carlo method for computing the 
marginal likelihood in nondecomposable Gaussian graphical models. Biometrika 92(2), 317-335.

Wang, H. and S. Li (2012). Efficient Gaussian graphical model determination under 
G-Wishart prior distributions. Electronic Journal of Statistics
6, 168-198.
}

\keyword{ package }
