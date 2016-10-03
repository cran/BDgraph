\name{surveyData}
\alias{surveyData}
\docType{data}
\title{ Labor force survey data }

\description{ 
The survey dataset concerns 1002 males in the U.S labor force, described by Hoff (2007). 
The seven observed variables which have been measured on various scales are as follow: 
the income (income), degree (degree), the number of children (children), parents income (pincome), 
parents degree (pdegree), number of parents children (pchildren), and age (age).  
}

\usage{ data( surveyData ) }

\format{ The format is a matrix with 1002 rows (number of individuals) and 7 column (number of variables). }

\references{
Hoff, P. (2007). Extending the rank likelihood for semiparametric copula estimation, \emph{The Annals of Applied Statistics, 1(1)}, 265-283.

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108} 
}

\examples{
\dontrun{
data( surveyData )
dim( surveyData )
head( surveyData )
boxplot( surveyData )
}
}

\keyword{ datasets }
