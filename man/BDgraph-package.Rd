\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}

\title{ Graph selection based on birth-death MCMC }

\description{
The R package BDgraph provides a statistical tool for Bayesian structure learning in 
undirected graphical models based on birth-death MCMC methodology.      
}
\details{
  This package provides an implementation of the procedures described in Mohammadi and Wit (2013). 
The main function is 'bdgraph' which is a birth-death MCMC algorithm for Bayesian inference 
in graphical models.

  Functions:
  \preformatted{
bdgraph          birth-death MCMC algorithm for graphical models
bdgraph.sim      Synthetic graph data generator 
bdgraph.npt      Nonparametric transfer 
CellSignal       A flow cytometry dataset 
compare          Comparing the result 
geneExpression   Human gene expression dataset
I.g              Normalizing constant of G-Wishart distribution
phat             Posterior edge inclusion probabilities
plot.bdgraph     Plot function for "bdgraph" output
plotcoda         Convergence plot
plotroc          ROC plot
print.bdgraph    Print function for "bdgraph" output
prob             Posterior probabilities of the graphs
rgwish           Sampling from G-Wishart distribution
rwish            Sampling from Wishart distribution
select           Selecting the best graph
surveyData       Labor force survey data
summary.bdgraph  Summary function for "bdgraph" output
traceplot        Trace plot of graph size from "bdgraph" output
  }
}

\author{Abdolreza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit}

\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}


Lenkoski, A. (2013). A direct sampler for G-Wishart variates, Stat 2, 119-128.

Wang, H. and S. Li (2012). Efficient Gaussian graphical model determination under 
G-Wishart prior distributions. Electronic Journal of Statistics 6, 168-198.
}

\keyword{ package }
