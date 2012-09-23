\name{plot_links}
\alias{plot_links}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of posterior distribution for graph according to number of their links
}
\description{
This function is for ploting the posterior distribution of all possible graphs according to number
of their links.
}
\usage{
plot_links(output, min = 0, max = (nrow(output$As[[1]])) * (nrow(output$As[[1]]) - 1)/2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
A list which is the result of BDMCMC algorithm
}
  \item{min}{
Minimum value for number of links 
}
  \item{max}{
Maximum value for number of links
}
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, Journal of the Royal Statistical Society: Series B,
submitted.
}
\author{
Abdolreza Mohammadi and Ernst Wit
}

\seealso{
\code{\link{BDMCMC}}
}
\examples{
n=200; p=8
truK=diag(p)
for (i in 1:(p-1)) truK[i,i+1]=truK[i+1,i]=0.5
truK[1,p]=truK[p,1]=0.4
data=mvrnorm(n,c(rep(0,p)),solve(truK))
\dontshow{
  data(output.high)
  output <- output.high}
\dontrun{output = BDMCMC(data) }
plot_links(output)
}
