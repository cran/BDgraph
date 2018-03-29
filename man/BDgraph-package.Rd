\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}

\alias{ compute_tp_fp }
\alias{ compute_measures }

\alias{ hill_climb_mpl }
\alias{ local_mb_hc }
\alias{ global_hc }
\alias{ log_mpl_disrete }

\alias{ hill_climb_mpl_binary }
\alias{ local_mb_hc_binary }
\alias{ global_hc_binary }
\alias{ log_mpl_binary }

\alias{ detect_cores }
\alias{ check.os }

\alias{ sample_ug }
\alias{ generate_clique_factors }
\alias{ calc_joint_dist }

\docType{package}

\title{ Bayesian Structure Learning in Graphical Models }

\description{
The \code{R} package \pkg{BDgraph} provides statistical tools for Bayesian structure learning in undirected graphical models. %based on birth-death MCMC method.      
The package is implemented the recent improvements in the Bayesian graphical models literature, including Mohammadi and Wit (2015), Mohammadi et al. (2017), Dobra and Mohammadi (2017).
The computationally intensive tasks of the package are implemented in parallel using \pkg{OpenMP} in \code{C++} and interfaced with \code{R}, to speed up the computations. 
Besides, the package contains several functions for simulation and visualization, as well as three multivariate datasets taken from the literature.
}
\details{
The package includes 10 main functions: 

\preformatted{
bdgraph        Search algorithm in graphical models
bdgraph.mpl    Search algorithm in graphical models using marginal pseudo-likehlihood
bdgraph.sim    Graph data generator 
bdgraph.npn    Nonparametric transfer 
compare        Graph structure comparison 
plinks         Estimated posterior link probabilities
plotcoda       Convergence plot
plotroc        ROC plot
rgwish         Sampling from G-Wishart distribution
rwish          Sampling from Wishart distribution
select         Graph selection
traceplot      Trace plot of graph size
}
}

\section{How to cite this package}{
Whenever using this package, please cite as
\preformatted{
Mohammadi A. and E. C. Wit (2017). BDgraph: Bayesian Structure Learning in 
Graphical Models using Birth-Death MCMC, R package version 3.45, 
https://CRAN.R-project.org/package=BDgraph
}
}

\author{Reza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 

Dobra, A. and A. Mohammadi (2017). Loglinear Model Selection and Human Mobility, \emph{arXiv preprint arXiv:1711.02623}

Mohammadi, A., Massam H., and G. Letac (2017). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416} 

Mohammadi, A. and A. Dobra (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16

Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2:119-128
}

\keyword{ Package, BDgraph, Bayesian Structure Learning, Gaussian Graphical Models, Gaussian Copula Graphical Models, BDMCMC, G-Wishart, Marginal Pseudo-likelihood }
