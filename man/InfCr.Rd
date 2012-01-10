\name{InfCr}
\alias{InfCr}

\title{Information Criteria}

\description{ 
\code{InfCr} can be used to calculate the Akaike information criterion (AIC), corrected AIC and Bayesian information criterion.}


\usage{

InfCr(x,cr="AIC")

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{cr}{The criterion to be used. Default is \code{AIC}; alternatives are \code{AICc} and \code{BIC}.} 

}


\details{
When comparing fitted objects, the smaller the chosen criterion, the better the fit.
}


\value{

It returns a numeric value corresponding to the selected criterion.
}



\author{

Maintainer: Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

}




\seealso{
\code{\link{AT}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{information criterion} 








