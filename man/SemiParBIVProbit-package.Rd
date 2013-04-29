\name{SemiParBIVProbit-package}
\alias{SemiParBIVProbit-package}
\docType{package}

\title{Semiparametric Bivariate Probit Modelling}

\description{\code{SemiParBIVProbit} provides a function for bivariate probit modelling with semiparametric 
predictors, including linear, nonlinear and random effects.
}

\details{\code{SemiParBIVProbit} provides a function for flexible bivariate probit modelling, in the presence of 
correlated error equations, endogeneity or non-random sample selection. The underlying representation and 
estimation of the model is based on a penalized regression spline approach, with automatic smoothness selection. The 
numerical routine carries out function minimization using a trust region algorithm from the package \code{trust} in combination with
an adaptation of a low level smoothness selection fitting procedure from the package \code{mgcv}. 
 
\code{\link{SemiParBIVProbit}} supports the use of many smoothers as extracted from \code{mgcv}. Scale invariant tensor product smooths
are not currently supported. Estimation is by penalized maximum likelihood with automatic smoothness selection achieved 
by using the approximate Un-Biased Risk Estimator (UBRE). 

Confidence intervals for smooth components are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided. Functions \code{\link{plot.SemiParBIVProbit}} and 
\code{\link{summary.SemiParBIVProbit}} extract such information from a fitted \code{\link{SemiParBIVProbitObject}}. Model/variable 
selection is also possible via the use of shrinakge smoothers or information criteria. 

The use of parametric or nonparametric random effects is also allowed for.

}


\author{
Giampiero Marra (University College London, Department of Statistical Science) and Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics) with contributions 
from Georgios Papageorgiou (Imperial College London, Department of Epidemiology and Biostatistics)

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}
}

\references{

Marra G., Papageorgiou G. and Radice R. (in press), Estimation of a Semiparametric Recursive Bivariate Probit Model with Nonparametric Mixing. \emph{Australian & New Zealand Journal of Statistics}.

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (submitted), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, invited revision.


}


\seealso{
\code{\link{SemiParBIVProbit}}
}



\keyword{endogeneity} 
\keyword{package}
\keyword{non-random sample selection} 
\keyword{semiparametric bivariate probit modelling} 
\keyword{smooth} 
\keyword{nonparametric random effects}
\keyword{regression}
\keyword{variable selection}

