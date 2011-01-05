\name{SemiParBIVProbit-package}
\alias{SemiParBIVProbit-package}
\docType{package}

\title{Semiparametric Bivariate Probit Modelling}

\description{ \code{SemiParBIVProbit} provides a function for bivariate probit modelling where the linear predictors can be flexibly specified 
by using parametric and nonparametric components.}

\details{ \code{SemiParBIVProbit} provides a function for flexible bivariate probit modelling. The underlying representation and 
estimation of the model is based on a penalized regression spline approach, with automatic smoothness selection. The algorithm is based on 
the package \code{trust} (which carries out function minimization using a trust region algorithm) and 
an adaptation of a low level smoothness selection fitting function of the package \code{mgcv}. 
 
\code{\link{SemiParBIVProbit}} supports the use of many smoothers as extracted from \code{mgcv}. Scale invariant tensor product smooths
are not currently supported. Estimation is by penalized maximum likelihood with automatic smoothness selection by approximate 
Un-Biased Risk Estimator (UBRE) score. 

Confidence intervals for smooth components are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided. Functions \code{\link{plot.SemiParBIVProbit}} and 
\code{\link{summary.SemiParBIVProbit}} extract such information from a fitted \code{\link{SemiParBIVProbitObject}}. Model/variable selection is also 
possible via the use of shrinakge smoothers or information criteria. 

}


\author{
Giampiero Marra and Rosalba Radice

Maintainer: Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

}

\references{

Marra G. and Radice R. (in press), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. The Canadian 
Journal of Statistics.

}


\seealso{
\code{\link{SemiParBIVProbit}}
}




\keyword{package}
\keyword{semiparametric bivariate probit modelling} 
\keyword{smooth} 
\keyword{regression}
\keyword{model selection}

