\name{Discrete Weibull}

\alias{ddweibull}
\alias{pdweibull}
\alias{qdweibull}
\alias{rdweibull}

\title{ The Discrete Weibull Distribution (Type 1) }

\description{
Density, distribution function, quantile function and random generation for the discrete Weibull distribution (type I) with parameters \eqn{q} and \eqn{\beta}.
}

\usage{
ddweibull( x, q = exp( -1 ), beta = 1, zero = TRUE )
pdweibull( x, q = exp( -1 ), beta = 1, zero = TRUE )
qdweibull( p, q = exp( -1 ), beta = 1, zero = TRUE )
rdweibull( n, q = exp( -1 ), beta = 1, zero = TRUE )
}

\arguments{
  \item{x}{ vector of quantiles.}
  \item{p}{ vector of probabilities. }
  \item{q, beta}{ shape and scale parameters, the latter defaulting to 1. }
  
  \item{zero}{ logical; if \code{TRUE} (default), the support contains \eqn{0}; \code{FALSE} otherwise. }
  \item{n}{ number of observations. If \code{length(n) > 1}, the length is taken to be the number required. }
}

\details{
The discrete Weibull distribution has density given by 
\deqn{
f(x) = q^{x^\beta} - q^{(x+1)^\beta}, x = 0, 1, 2, \ldots
}{
f(x) = q^x^\beta - q^(x+1)^\beta, x = 0, 1, 2, \ldots
}
For the case \code{zero} = \code{FALSE}:
\deqn{
f(x) = q^{(x-1)^\beta} - q^{x^\beta}, x = 1, 2, \ldots
}{
f(x) = q^(x-1)^\beta - q^x^\beta, x = 1, 2, \ldots
}

Cumulative distribution function
\deqn{
F(x) = 1-q^{(x+1)^\beta}
}{
F(x) = 1-q^(x+1)^\beta
}
For the case \code{zero} = \code{FALSE}, \code{x+1} should replaced by \code{x}.
}

\value{
\code{ddweibull} gives the density, \code{pdweibull} gives the distribution function, \code{qdweibull} gives the quantile function, and \code{rdweibull} generates random values.
}

\references{
Nakagawa, T. and Osaki, S. (1975). The Discrete Weibull Distribution. \emph{IEEE Transactions on Reliability}, R-24, 300-301, \doi{10.1109/TR.1975.5214915}
} 

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}, Pariya Behrouzi, Veronica Vinciotti }

\seealso{ \code{\link{dweibull}}, \code{\link{bdw.reg}}, \link{bdgraph.dw} }

\examples{
n    = 1000
q    = 0.4
beta = 0.8

set.seed( 7 )

rdw = rdweibull( n = n, q = q, beta = beta )

plot( prop.table( table( rdw ) ), type = "h", col = "gray50" )

x = 0:max( rdw )

lines( x, ddweibull( x = x, q = q, beta = beta ), type = "o", col = "blue", lwd = 2 )

hist( pdweibull( x = rdw, q = q, beta = beta ) )

plot( ecdf( rdw ) )
lines( x, pdweibull( x, q = q, beta = beta ), col = "blue", lwd = 2, type = "s" )
}

\keyword{distribution}

