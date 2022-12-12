\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}

\alias{ get_Ds_tgm_R }
\alias{ get_Ts_R }
\alias{ update_tu_R }
\alias{ update_mu_R }

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

\alias{ get_bounds_dw }
\alias{ ddweibull_reg }
\alias{ log_post_cond_dw }
\alias{ near_positive_definite }

\docType{package}

\title{ Bayesian Structure Learning in Graphical Models }

\description{
The \code{R} package \pkg{BDgraph} provides statistical tools for Bayesian structure learning in undirected graphical models for continuous, ordinal/count/dicrete, binary, and mixed data. The package is implemented the recent improvements in the Bayesian graphical models' literature, including Mohammadi and Wit (2015), Mohammadi et al. (2021), Mohammadi et al. (2017), and Dobra and Mohammadi (2018). To speed up the computations, the intensive tasks of the package are implemented in parallel using \pkg{OpenMP} in \code{C++} and interfaced with \code{R}. Besides, the package contains several functions for simulation and visualization, as well as several multivariate datasets taken from the literature.
}

\section{How to cite this package}{
To cite \pkg{BDgraph} in publications use:

Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03}
}

\author{
    Reza Mohammadi [aut, cre] (<https://orcid.org/0000-0001-9538-0648>), \cr
    Ernst Wit [aut] (<https://orcid.org/0000-0002-3671-9610>), \cr
    Adrian Dobra [ctb] (<https://orcid.org/0000-0001-7793-2197>). \cr
	Maintainer: Reza Mohammadi \email{a.mohammadi@uva.nl}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03}

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 

Mohammadi, A., et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Vinciotti, V., Behrouzi, P., and Mohammadi, R. (2022) Bayesian structural learning of microbiota systems from count metagenomic data, \emph{arXiv preprint}, \doi{10.48550/arXiv.2203.10118}

Dobra, A. and Lenkoski, A. (2011). Copula Gaussian graphical models and their application to modeling functional disability data, \emph{The Annals of Applied Statistics}, 5(2A):969-93, \doi{10.1214/10-AOAS397}

Dobra, A., et al. (2011). Bayesian inference for general Gaussian graphical models with application to multivariate lattice data. \emph{Journal of the American Statistical Association}, 106(496):1418-33, \doi{10.1198/jasa.2011.tm10465}

Mohammadi, A. and Dobra, A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16

Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2(1):119-28, \doi{10.1002/sta4.23}

Pensar, J. et al (2017) Marginal pseudo-likelihood learning of discrete Markov network structures, \emph{Bayesian Analysis}, 12(4):1195-215, \doi{10.1214/16-BA1032}
}

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{bdgraph.dw}}, \code{\link{bdgraph.sim}}, \code{\link{compare}}, \code{\link{rgwish}} }

\examples{
\dontrun{
library( BDgraph )

set.seed( 10 )

# Generating multivariate normal data from a 'scale-free' graph
data.sim <- bdgraph.sim( n = 100, p = 10, graph = "scale-free", vis = TRUE )

# Running algorithm based on GGMs
bdgraph.obj <- bdgraph( data = data.sim, iter = 5000 )

summary( bdgraph.obj )

# To compare the result with true graph
compare( bdgraph.obj, data.sim, main = c( "Target", "BDgraph" ), vis = TRUE )

# Confusion Matrix
conf.mat( actual = data.sim, pred = bdgraph.obj )

conf.mat.plot( actual = data.sim, pred = bdgraph.obj )

# Running algorithm based on GGMs and marginal pseudo-likelihood
bdgraph.mpl.obj <- bdgraph.mpl( data = data.sim, iter = 5000 )

summary( bdgraph.mpl.obj )

# Confusion Matrix
conf.mat( actual = data.sim, pred = bdgraph.mpl.obj )

conf.mat.plot( actual = data.sim, pred = bdgraph.mpl.obj )

# To compare the results of both algorithms with true graph
compare( list( bdgraph.obj, bdgraph.mpl.obj ), data.sim, 
         main = c( "Target", "BDgraph", "BDgraph_mpl" ), vis = TRUE )
}
}

\keyword{package}
