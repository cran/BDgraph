\name{churn}
\alias{churn}
\docType{data}
\title{ Churn data set }

\description{ 
The data can be downloaded from IBM Sample Data Sets. Customer \emph{churn} occurs when customers stop doing business with a company, also known as customer attrition. The data set contains \eqn{3333} rows (customers) and \eqn{20} columns (features). The "Churn" column is our target which indicate whether customer churned (left the company) or not.   
}

\usage{ data( churn ) }

\format{ The data set, as a data frame, contains \eqn{3333} rows (customers) and \eqn{20} columns (features).  }

\references{
Larose, D. T. and Larose, C. D. (2014). Discovering knowledge in data: an introduction to data mining. \emph{John Wiley & Sons}.
}

\examples{
data( churn )

summary( churn )
}

\keyword{ datasets }
