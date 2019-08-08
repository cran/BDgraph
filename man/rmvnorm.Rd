\name{rmvnorm}
\alias{rmvnorm}

\title{ Generate data from the multivariate Normal distribution }

\description{
    Random generation function from the multivariate Normal distribution with mean equal to \eqn{mean} and covariance matrix \eqn{sigma}.
}

\usage{ 
    rmvnorm( n = 10, mean = rep( 0, length = ncol( sigma ) ), 
             sigma = diag( length( mean ) ) ) 
}

\arguments{
    \item{n}{ Number of observations. }
    \item{mean}{ Mean vector, default is \eqn{rep( 0, length = ncol( sigma ) )}. }
    \item{sigma}{ positive definite covariance matrix, default is \eqn{ diag( length( mean ) ) }. }
}

\value{ A numeric matrix with rows equal to \eqn{n} and columns equal to \eqn{length( mean )}. }

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{\code{\link{bdgraph.sim}}, \code{\link{rwish}}, \code{\link{rgwish}} }

\examples{
mean  <- c( 5, 20 )       
sigma <- matrix( c( 4, 2,
                    2, 5 ), 2, 2 )  # covariance matrix

sample <- rmvnorm( n = 500, mean = mean, sigma = sigma )
plot( sample )
}

\keyword{ distribution }
