\name{SemiParBIVProbit-package}
\alias{SemiParBIVProbit-package}
\docType{package}

\title{Semiparametric Bivariate Probit Modelling}

\description{\code{SemiParBIVProbit} provides a function for bivariate probit modelling with semiparametric 
predictors, including linear and nonlinear effects. The use of several bivariate copula distributions and of asymmetric link 
functions is also supported.
}

\details{\code{SemiParBIVProbit} provides a function for flexible bivariate probit modelling, in the presence of 
correlated error equations, endogeneity or non-random sample selection. The underlying representation and 
estimation of the model is based on a penalized regression spline approach, with automatic smoothness selection. The use of several
bivariate copula distributions is supported. The 
numerical routine carries out function minimization using a trust region algorithm from the package \code{trust} in combination with
an adaptation of a low level smoothness selection fitting procedure from the package \code{mgcv}. 
 
\code{\link{SemiParBIVProbit}} supports the use of many smoothers as extracted from \code{mgcv}. Scale invariant tensor product smooths
are not currently supported. Estimation is by penalized maximum likelihood with automatic smoothness selection achieved 
by using the approximate Un-Biased Risk Estimator (UBRE). 

Confidence intervals for smooth components are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided. Functions \code{\link{plot.SemiParBIVProbit}} and 
\code{\link{summary.SemiParBIVProbit}} extract such information from a fitted \code{\link{SemiParBIVProbitObject}}. Model/variable 
selection is also possible via the use of shrinakge smoothers or information criteria. 

Asymmetric link functions and tools for testing the hypothesis of absence of unobserved confounding are also available. 

}


\author{
Giampiero Marra (University College London, Department of Statistical Science) and Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics)

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}
}

\references{

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G., Radice R. and Missiroli S. (2014), Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Computational Statistics}, 29(3-4), 715-741.

Marra G., Radice R. and Filippou P. (submitted), Testing the Hypothesis of Exogeneity in Regression Spline Bivariate Probit Models.

Radice R., Marra G. and M. Wojtys (submitted), Copula Regression Spline Models for Binary Outcomes.


}


\seealso{
\code{\link{SemiParBIVProbit}}
}


\keyword{asymmetric link function} 
\keyword{endogeneity} 
\keyword{package}
\keyword{non-random sample selection} 
\keyword{semiparametric bivariate probit modelling} 
\keyword{smooth} 
\keyword{regression}
\keyword{variable selection}

