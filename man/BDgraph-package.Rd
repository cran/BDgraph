\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}

\title{ Graph selection based on birth-death MCMC }

\description{
The \code{R} package \pkg{BDgraph} provides statistical tools for Bayesian structure learning in undirected graphical models based on birth-death MCMC methodology.      
It implements the recent improvements in the Bayesian literature, including Mohammadi and Wit (2014) and Mohammadi et al. (2015).
}
\details{
The package includes 10 main functions: 

Functions:
\preformatted{
bdgraph        birth-death MCMC sampling algorithm for graphical models
bdgraph.sim    Synthetic graph data generator 
bdgraph.npt    Nonparametric transfer 
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
Mohammadi, A. and E. Wit (2014). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis} (Forthcoming). \url{http://arxiv.org/abs/1210.5371v6}

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108}

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling 1}, 225-230

Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat 2}, 119-128
}

\keyword{ package }
