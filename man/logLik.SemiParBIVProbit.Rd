\name{logLik.SemiParBIVProbit}
\alias{logLik.SemiParBIVProbit}

\title{Extract the log likelihood for a fitted SemiParBIVProbit}


\description{It extracts the log-likelihood for a fitted \code{SemiParBIVProbit} model.

}

\usage{
\method{logLik}{SemiParBIVProbit}(object, ...)
}


\arguments{ 
  
\item{object}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{...}{Un-used for this function.}


}


\details{

Modification of the classic \code{logLik} which accounts for the estimated degrees of freedom used in \code{SemiParBIVProbit} objects.
This function is provided so that information criteria work correctly with \code{SemiParBIVProbit} objects by using the correct degrees 
of freedom. 

}



\value{

Standard \code{\link{logLik}} object.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Radice R., Marra G. and M. Wojtys (submitted), Copula Regression Spline Models for Binary Outcomes.

}


\seealso{
\code{\link{AIC}}, \code{\link{BIC}}
}



\examples{
## see examples for SemiParBIVProbit
}

\keyword{logLik} 
\keyword{AIC} 
\keyword{BIC}






