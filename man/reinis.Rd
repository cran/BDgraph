\name{reinis}
\alias{reinis}
\docType{data}

\title{ Risk factors of coronary heart disease }

\description{ 
The dataset consist of \eqn{6} discrete variables as the potential risk factors of coronary heart disease.
The data collected from \eqn{1841} men employed of a car factory in Czechoslovakia (Reinis et al. 1981). }

\usage{ data( reinis ) }

\format{ The format is a matrix with \eqn{1841} rows (number of individuals) and \eqn{6} column (number of variables). }

\references{
Edwards and Havranek (1985). A fast procedure for model search in multidimensional contingency tables, \emph{Biometrika}, 72:339-351

Reinis et al (1981). Prognostic significance of the risk profile in the prevention of coronary heart disease, \emph{Bratis. lek. Listy}, 76:137-150

Mohammadi, A. and Dobra, A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16
}

\examples{
data( reinis )

summary( reinis )
}

\keyword{ datasets }
