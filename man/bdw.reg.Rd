\name{bdw.reg}

\alias{bdw.reg}

\title{	Bayesian estimation of (zero-inflated) Discrete Weibull regression }
\description{
  Bayesian estimation of the parameters for Discrete Weibull (DW) regression. The conditional distribution of the response given the predictors is assumed to be DW with parameters q and beta, dependent on the predictors, and,  with an additional parameter pi under zero inflation.
}
\usage{
bdw.reg( data, formula = NA, iter = 5000, burnin = NULL, 
         dist.q = dnorm, dist.beta = dnorm, 
         par.q = c( 0, 1 ), par.beta = c( 0, 1 ), par.pi = c( 1, 1 ), 
         initial.q = NULL, initial.beta = NULL, initial.pi = NULL, 
         ZI = FALSE, scale.proposal = NULL, adapt = TRUE, print = TRUE )
}

\arguments{
	\item{data}{\code{data.frame} or \code{matrix} corresponding to the data, containing the variables in the model.}

	\item{formula}{ object of class \link{formula} as a symbolic description of the model to be fitted. For the case of \code{data.frame}, it is taken as the model frame (see \code{\link{model.frame})}.}

	\item{iter}{ number of iterations for the sampling algorithm.}
	\item{burnin}{number of burn-in iterations for the sampling algorithm.}

  \item{dist.q}{ Prior density for the regression coefficients associated to the parameter \code{q}. The default is a Normal distribution (\code{dnorm}). 
  Any density function which has two parameters and can support the \code{log = TRUE} flag can be used, e.g. \code{\link{dnorm}}, \code{\link{dlnorm}}, \code{\link{dunif}} etc. 
}
  \item{dist.beta}{ Prior density for the regression coefficients associated to the parameter \code{beta}. The default is a Normal distribution (\code{\link{dnorm}}). 
  Any density function which has two parameters and can support the \code{log = TRUE} flag can be used, e.g. \code{\link{dnorm}}, \code{\link{dlnorm}}, \code{\link{dunif}} etc. 
}
  \item{par.q}{ vector of length two corresponding to the parameters of \code{dist.q}. }
  \item{par.beta}{ vector of length two corresponding to the parameters of \code{dist.beta}. }
  \item{par.pi}{ vector of length two corresponding to the parameters of the \code{\link{beta}} prior density on \code{pi}. }

  \item{initial.q, initial.beta, initial.pi}{ vector of initial values for the regression coefficients and for \code{pi} (if \code{ZI} = \code{TRUE}). }

  \item{ZI}{ 
		logical: if FALSE (default), the conditional distribution of the response given the predictors is assumed to be DW with parameters \code{q} and \code{beta}. 
		If TRUE, a zero-inflated DW distribution will be applied.
	}

  \item{scale.proposal }{ scale of the proposal function. Setting to lower values results in an increase in the acceptance rate of the sampler. }
  \item{adapt }{ logical: if TRUE (default), the proposals will be adapted. If FALSE, no adapting will be applied. }
	\item{print }{ logical: if TRUE (default), tracing information is printed.}
}

\details{
The regression model uses a logit link function on \code{q} and a log link function on \code{beta}, the two parameters of a DW distribution, with probability mass function given by 
\deqn{
DW(y) = q^{y^\beta} - q^{(y+1)^\beta}, y = 0, 1, 2, \ldots
}{
DW(y) = q^y^\beta - q^(y+1)^\beta, y = 0, 1, 2, \ldots
}

For the case of zero inflation (\code{ZI} = \code{TRUE}), a zero-inflated DW is considered:
\deqn{
f(y) = ( 1 - pi) I(y = 0 ) + pi DW(y)
}{
f(y) = ( 1 - pi )  I(y = 0 ) + pi DW(y)
}
where \eqn{0 \leq pi \leq 1} and \eqn{I(y = 0 )} is an indicator for the point mass at zero for the response \code{y}. 
}

\value{
	\item{sample}{ MCMC samples }
	\item{q.est}{ posterior estimates of \code{q} }
	\item{beta.est}{ posterior estimates of \code{beta} }
	\item{pi.est}{ posterior estimates of \code{pi} }
	\item{accept.rate}{ acceptance rate of the MCMC algorithm }
}

\references{
Vinciotti, V., Behrouzi, P., and Mohammadi, R. (2022) Bayesian structural learning of microbiota systems from count metagenomic data, \emph{arXiv preprint}, \doi{10.48550/arXiv.2203.10118}

Peluso, A., Vinciotti, V., and Yu, K. (2018) Discrete Weibull generalized additive model: an application to count fertility, \emph{Journal of the Royal Statistical Society: Series C}, 68(3):565-583, \doi{10.1111/rssc.12311}  

Haselimashhadi, H., Vinciotti, V. and Yu, K. (2018) A novel Bayesian regression model for counts with an application to health data, \emph{Journal of Applied Statistics,} 45(6):1085-1105, \doi{10.1080/02664763.2017.1342782}
}

\author{ Veronica Vinciotti, Reza Mohammadi \email{a.mohammadi@uva.nl}, and Pariya Behrouzi }

\seealso{ \code{\link{bdgraph.dw}}, \link{bdgraph}, \link{ddweibull}, \code{\link{bdgraph.sim}} }

\examples{
\dontrun{
# - - Example 1

q    = 0.6
beta = 1.1
n    = 500

y = BDgraph::rdweibull( n = n, q = q, beta = beta )

output = bdw.reg( data = y, y ~ ., iter = 5000 )

output $ q.est
output $ beta.est

traceplot( output $ sample[ , 1 ], acf = T, pacf = T )
traceplot( output $ sample[ , 2 ], acf = T, pacf = T )

# - - Example 2

q    = 0.6
beta = 1.1
pii  = 0.8
n    = 500

y_dw = BDgraph::rdweibull( n = n, q = q, beta = beta )
z = rbinom( n = n, size = 1, prob = pii ) 
y = z * y_dw

output = bdw.reg( data = y, iter = 5000, ZI = TRUE )

output $ q.est
output $ beta.est
output $ pi.est

traceplot( output $ sample[ , 1 ], acf = T, pacf = T )
traceplot( output $ sample[ , 2 ], acf = T, pacf = T )
traceplot( output $ sample[ , 3 ], acf = T, pacf = T )

# - - Example 3

theta.q    = c( 0.1, -0.1, 0.34 )  # true parameter
theta.beta = c( 0.1, -.15, 0.5  )  # true parameter

n  = 500

x1 = runif( n = n, min = 0, max = 1.5 )
x2 = runif( n = n, min = 0, max = 1.5 )

reg_q = theta.q[ 1 ] + x1 * theta.q[ 2 ] + x2 * theta.q[ 3 ]
q     = 1 / ( 1 + exp( - reg_q ) )

reg_beta = theta.beta[ 1 ] + x1 * theta.beta[ 2 ] + x2 * theta.beta[ 3 ]
beta     = exp( reg_beta )

y = BDgraph::rdweibull( n = n, q = q, beta = beta )

data = data.frame( x1, x2, y ) 

output = bdw.reg( data, y ~. , iter = 5000 )

# - - Example 4

theta.q    = c( 1, -1, 0.8 )  # true parameter
theta.beta = c( 1, -1, 0.3 )  # true parameter
pii = 0.8

n  = 500

x1 = runif( n = n, min = 0, max = 1.5 )
x2 = runif( n = n, min = 0, max = 1.5 )

reg_q = theta.q[ 1 ] + x1 * theta.q[ 2 ] + x2 * theta.q[ 3 ]
q     = 1 / ( 1 + exp( - reg_q ) )

reg_beta = theta.beta[ 1 ] + x1 * theta.beta[ 2 ] + x2 * theta.beta[ 3 ]
beta     = exp( reg_beta )

y_dw = BDgraph::rdweibull( n = n, q = q, beta = beta )
z    = rbinom( n = n, size = 1, prob = pii ) 
y    = z * y_dw

data = data.frame( x1, x2, y ) 

output = bdw.reg( data, y ~. , iter = 5000 )
}	  
}

\keyword{Discrete Weibull}
\keyword{sampling algorithms}
\keyword{iteration}
