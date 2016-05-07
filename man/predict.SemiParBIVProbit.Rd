\name{predict.SemiParBIVProbit}
\alias{predict.SemiParBIVProbit}

\title{Prediction from a fitted SemiParBIVProbit model}

\description{

It takes a fitted \code{SemiParBIVProbit},\code{copulaReg},\code{copulaSampleSel},\code{SemiParTRIVProbit} object produced by \code{SemiParBIVProbit()},\code{copulaReg()},\code{copulaSampleSel()},\code{SemiParTRIVProbit()} and, 
for each equation, produces predictions 
for a new set of values of the model covariates or the original values used for the model fit. 
Standard errors of predictions can be produced and are based on the posterior distribution of the model coefficients. This function is a 
wrapper for \code{predict.gam()} in \code{mgcv}. Please see the documentation of \code{predict.gam()} for full details. 
}

\usage{

\method{predict}{SemiParBIVProbit}(object, eq, ...)


}


\arguments{ 

\item{object}{A fitted code{SemiParBIVProbit}/\code{copulaReg}/\code{copulaSampleSel} object as 
              produced by \code{SemiParBIVProbit()}/\code{copulaReg()}/\code{copulaSampleSel()}.}

\item{eq}{The equation to be considered for prediction.}

\item{...}{Other arguments as in \code{predict.gam()} in \code{mgcv}.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\section{WARNINGS}{

When \code{type = "response"} (which gives predictions on the scale of the response variable). For 
the case of continuous responses this function will NOT produce correct predictions for the outcome variable (except for the Gaussian case). This
is because for all distributions (except the Gaussian) implemented in this package the distribution parameters 
determine the mean and variance
through functions of them.

} 






\seealso{
\code{\link{SemiParBIVProbit}}, \code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{plot.SemiParBIVProbit}}
}


\keyword{prediction}







