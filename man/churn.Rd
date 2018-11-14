\name{churn}
\alias{churn}
\docType{data}
\title{ Churn data set }

\description{ 
The data can be downloaded from IBM Sample Data Sets. Customer \emph{churn} occurs when customers stop doing business with a company, also known as customer attrition. The data set contains \eqn{3333} rows (customers) and \eqn{20} columns (features). The "Churn" column is our target which indicate whether customer churned (left the company) or not.   
}

\usage{ data( churn ) }

\format{ 

The \code{churn} dataset, as a data frame, contains \eqn{3333} rows (customers) and \eqn{20} columns (variables/features). The \eqn{20} variables are:
 \itemize{

\item \code{State}: Categorical, for the \eqn{50} states and the District of Columbia.
\item \code{Account.Length}: count, how long account has been active.
\item \code{Area.Code}: Categorical.
\item \code{Int.l.Plan}: Categorical, yes or no, international plan.
\item \code{VMail.Plan}: Categorical, yes or no, voice mail plan.
\item \code{VMail.Message}: Count, number of voice mail messages.
\item \code{Day.Mins}: Continuous, minutes customer used service during the day.
\item \code{Day.Calls}: Count, total number of calls during the day.
\item \code{Day.Charge}: Continuous, total charge during the day.
\item \code{Eve.Mins}: Continuous, minutes customer used service during the evening.
\item \code{Eve.Calls}: Count, total number of calls during the evening.
\item \code{Eve.Charge}: Continuous, total charge during the evening.
\item \code{Night.Mins}: Continuous, minutes customer used service during the night.
\item \code{Night.Calls}: Count, total number of calls during the night.
\item \code{Night.Charge}: Continuous, total charge during the night.
\item \code{Intl.Mins}: Continuous, minutes customer used service to make international calls.
\item \code{Intl.Calls}: Count, total number of international calls.
\item \code{Intl.Charge}: Continuous, total international charge.
\item \code{CustServ.Calls}: Count, number of calls to customer service.
\item \code{Churn}: Categorical, True or False. Indicator of whether the customer has left the company (True or False).    
  }
}

\references{
Larose, D. T. and Larose, C. D. (2014). Discovering knowledge in data: an introduction to data mining. \emph{John Wiley & Sons}.
}

\examples{
data( churn )

summary( churn )
}

\keyword{ datasets }
