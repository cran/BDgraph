\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}

\title{ Bayesian Structure Learning in Graphical Models }

\description{
The \code{R} package \pkg{BDgraph} provides statistical tools for Bayesian structure learning in undirected graphical models. %based on birth-death MCMC method.      
The package is implemented the recent improvements in the Bayesian literature, including Mohammadi and Wit (2015) and Mohammadi et al. (2015).
The computationally intensive tasks of the package is implemented in \code{C++} and interfaced with \code{R}, to speed up the computations. 
Besides, the package contains several functions for simulation and visualization, as well as two multivariate datasets taken from the literature.
}
\details{
The package includes 10 main functions: 

\preformatted{
bdgraph        Search algorithm in graphical models
bdgraph.sim    Synthetic graph data generator 
bdgraph.npn    Nonparametric transfer 
compare        Comparing the results 
plinks         Estimated posterior link probabilities
plotcoda       Convergence plot
plotroc        ROC plot
rgwish         Sampling from G-Wishart distribution
select         Graph selection
traceplot      Trace plot of graph size
}
}

\author{Abdolreza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 

Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2:119-128
}

\keyword{ package }
