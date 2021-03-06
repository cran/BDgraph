\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}

\alias{ compute_tp_fp }
\alias{ compute_measures }

\alias{ get_graph }
\alias{ get_g_prior }
\alias{ get_g_start }
\alias{ get_K_start }
\alias{ get_S_n_p }
\alias{ get_cores }

\alias{ hill_climb_mpl }
\alias{ local_mb_hc }
\alias{ global_hc }
\alias{ log_mpl_disrete }

\alias{ hill_climb_mpl_binary }
\alias{ local_mb_hc_binary }
\alias{ global_hc_binary }
\alias{ log_mpl_binary }

\alias{ detect_cores }

\alias{ sample_ug }
\alias{ generate_clique_factors }
\alias{ calc_joint_dist }

\docType{package}

\title{ Bayesian Structure Learning in Graphical Models }

\description{
The \code{R} package \pkg{BDgraph} provides statistical tools for Bayesian structure learning in undirected graphical models for continuous, count, binary, and mixed data. The package is implemented the recent improvements in the Bayesian graphical models literature, including Mohammadi and Wit (2015), Mohammadi et al. (2017), Dobra and Mohammadi (2018), and Letac et al. (2018). The computationally intensive tasks of the package are implemented in parallel using \pkg{OpenMP} in \code{C++} and interfaced with \code{R}, to speed up the computations. Besides, the package contains several functions for simulation and visualization, as well as several multivariate datasets taken from the literature.
}

\section{How to cite this package}{
To cite \pkg{BDgraph} in publications use:
\preformatted{
Mohammadi, R. and Wit, E. C. (2019). BDgraph: An R Package for Bayesian 
Structure Learning in Graphical Models, Journal of Statistical Software, 
89(3):1-30. doi:10.18637/jss.v089.i03
}
}

\author{

	Reza Mohammadi \cr
	Amsterdam Business School \cr
	University of Amsterdam

	Maintainer: Reza Mohammadi \email{a.mohammadi@uva.nl}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30. 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A., et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Dobra, A. and Lenkoski, A. (2011). Copula Gaussian graphical models and their application to modeling functional disability data, \emph{The Annals of Applied Statistics}, 5(2A):969-93

Dobra, A., et al. (2011). Bayesian inference for general Gaussian graphical models with application to multivariate lattice data. \emph{Journal of the American Statistical Association}, 106(496):1418-33

Mohammadi, A. and Dobra, A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16

Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2(1):119-28

Pensar, J. et al (2017) Marginal pseudo-likelihood learning of discrete Markov network structures, \emph{Bayesian Analysis}, 12(4):1195-215
}

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{bdgraph.sim}}, \code{\link{compare}}, \code{\link{rgwish}} }

\examples{
\dontrun{
library( BDgraph )

# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 70, p = 6, size = 7, vis = TRUE )

# Running algorithm based on GGMs
bdgraph.obj <- bdgraph( data = data.sim, iter = 5000 )

summary( bdgraph.obj )

# To compare the result with true graph
compare( data.sim, bdgraph.obj, main = c( "Target", "BDgraph" ), vis = TRUE )

# Running algorithm based on GGMs and marginal pseudo-likelihood
bdgraph.obj_mpl <- bdgraph.mpl( data = data.sim, iter = 5000 )

summary( bdgraph.obj_mpl )

# To compare the results of both algorithms with true graph
compare( data.sim, bdgraph.obj, bdgraph.obj_mpl, 
         main = c( "Target", "BDgraph", "BDgraph_mpl" ), vis = TRUE )
}
}

\keyword{package}
