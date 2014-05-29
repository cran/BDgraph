\name{surveyData}
\alias{surveyData}
\docType{data}
\title{ Labor force survey data }

\description{ 
This survey data concerns 1002 males in the U.S labor force, described in Hoff (2007). 
The seven observed variables which have been measured on various scales are as follow: 
the income (income), degree (degree), number of children (children), parents income (pincome), 
parents degree (pdegree), number of parents children (pchildren), and age (age).  
}

\usage{ data(surveyData) }

\format{ The format is a matrix with 1002 rows (number of individuals) and 7 column (number of variables). }

\references{
Hoff, P. (2007). Extending the rank likelihood for semiparametric copula estimation,
The Annals of Applied Statistics, 1(1), 265-283.
}


%\author{Abdolreza Mohammadi and Ernst Wit}

\examples{
  data( surveyData )
  dim( surveyData )
  head( surveyData )
  boxplot( surveyData )
}

\keyword{ datasets }
