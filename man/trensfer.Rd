\name{transfer}
\alias{transfer}

\title{ transfer for count data }

\description{ Transfers count data, by counting the duplicated rows.  }

\usage{ transfer( r_data ) }

\arguments{
  \item{r_data}{ (\eqn{n \times p}{n x p}) \code{matrix} or a \code{data.frame} corresponding to the data (\eqn{n} is the sample size and \eqn{p} is the number of variables). }
}

\value{ (\eqn{n \times p+1}) matrix of transferred data, in which the last column is the frequency of duplicated rows.
}

\references{
Mohammadi, A. and Dobra, A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Adrian Dobra }

\seealso{ \code{\link{bdgraph.mpl}}, \code{\link{bdgraph.sim}} }

\examples{
# Generating multivariate binary data from a 'random' graph
data.sim <- bdgraph.sim( n = 12, p = 4, size = 4, type = "binary" )
r_data   <- data.sim $ data
r_data
   
# Transfer the data
transfer( r_data )
}

\keyword{data preprocessing}
