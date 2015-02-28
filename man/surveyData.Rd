\name{surveyData}
\alias{surveyData}
\docType{data}
\title{ Labor force survey data }

\description{ 
The survey data concerns 1002 males in the U.S labor force, described in Hoff (2007). 
The seven observed variables which have been measured on various scales are as follow: 
the income (income), degree (degree), number of children (children), parents income (pincome), 
parents degree (pdegree), number of parents children (pchildren), and age (age).  
}

\usage{ data(surveyData) }

\format{ The format is a matrix with 1002 rows (number of individuals) and 7 column (number of variables). }

\references{
Hoff, P. (2007). Extending the rank likelihood for semiparametric copula estimation, \emph{The Annals of Applied Statistics, 1(1)}, 265-283.

%Mohammadi, A. and E. Wit (2015). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108} 

%Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling}, 1:225-230
}

%\author{Abdolreza Mohammadi and Ernst Wit}

\examples{
data( surveyData )
dim( surveyData )
head( surveyData )
boxplot( surveyData )
}

\keyword{ datasets }
