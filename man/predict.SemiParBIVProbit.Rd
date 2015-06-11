\name{predict.SemiParBIVProbit}
\alias{predict.SemiParBIVProbit}

\title{Prediction from a fitted SemiParBIVProbit model}

\description{

It takes a fitted \code{SemiParBIVProbit} object produced by \code{SemiParBIVProbit()} and, for each equation, produces predictions 
for a new set of values of the model covariates or the original values used for the model fit. 
Standard errors of predictions can be produced and are based on the posterior distribution of the model coefficients. This function is a wrapper for \code{predict.gam()} in \code{mgcv}. Please see the documentation of \code{predict.gam()} for full details. 
}

\usage{

\method{predict}{SemiParBIVProbit}(object, eq, ...)


}


\arguments{ 

\item{object}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{eq}{The equation to be considered for prediction.}

\item{...}{Other arguments as in \code{predict.gam()} in \code{mgcv}.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\seealso{
\code{\link{SemiParBIVProbit}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}
}


\keyword{prediction}







