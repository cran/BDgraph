\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}

\title{ Graph selection based on birth-death MCMC }

\description{
The \code{R} package \pkg{BDgraph} provides statistical tools for Bayesian structure learning in undirected graphical models based on birth-death MCMC method.      
It implements the recent improvements in the Bayesian literature, including Mohammadi and Wit (2015) and Mohammadi et al. (2015).
}
\details{
The package includes 10 main functions: 

\preformatted{
bdgraph        birth-death MCMC sampling algorithm for graphical models
bdgraph.sim    Synthetic graph data generator 
bdgraph.npn    Nonparametric transfer 
compare        Comparing the result 
phat           Posterior link probabilities
plotcoda       Convergence plot
plotroc        ROC plot
rgwish         Sampling from G-Wishart distribution
select         Selecting the best graph
traceplot      Trace plot of graph size
}
}

\author{Abdolreza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108v2} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Gaussian Copula Graphical Models, \emph{Arxiv preprint arXiv:1501.04849v2} 

Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2:119-128
}

\keyword{ package }
