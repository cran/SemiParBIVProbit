\name{SemiParBIVProbit-package}
\alias{SemiParBIVProbit-package}
\docType{package}

\title{Semiparametric Bivariate Regression Models}

\description{
\code{\link{SemiParBIVProbit}} provides a function for fitting bivariate binary models with semiparametric 
predictors, including linear and nonlinear effects. Several bivariate copula distributions are supported. The dependence parameter of 
the bivariate distribution employed can be specified as a function of a semiparametric predictor as well. Smoothness selection is
achieved automatically and interval calculations are based on a Bayesian approach.

\code{SemiParBIVProbit} has been originally designed to deal with binary responses only. However, the use of 
bivariate equation systems in which one response is binary and the other is continuous is allowed for. Other models
will be incorporated from time to time on a case-study basis.


}

\details{
\code{SemiParBIVProbit} provides a function for fitting flexible bivariate binary models, in the presence of 
associated error equations, endogeneity, non-random sample selection or partial observability. The underlying representation and 
estimation of the model is based on a penalized likelihood-based regression spline approach, with automatic smoothness selection. Several
bivariate copula distributions are available. The 
numerical routine carries out function minimization using a trust region algorithm from the package \code{trust} in combination with
an adaptation of a low level smoothness estimation fitting procedure from the package \code{mgcv}. 
 
\code{SemiParBIVProbit} supports the use of many smoothers as extracted from \code{mgcv}. Estimation is by penalized 
maximum likelihood with automatic smoothness estimation achieved 
by using the approximate Un-Biased Risk Estimator (UBRE), which can also be viewed as an approximate AIC. 

Confidence intervals for smooth components and nonlinear functions of the model
parameters are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided and based on the approach
implemented in \code{mgcv}. Functions \code{\link{plot.SemiParBIVProbit}} and 
\code{\link{summary.SemiParBIVProbit}} extract such information from a fitted \code{\link{SemiParBIVProbit}} object. Model/variable 
selection is also possible via the use of shrinakge smoothers or information criteria. 

Tools for testing the hypothesis of uncorrelated error equations/absence of unobserved 
confounding/absence of endogeneity/absence of non-random sample selection are available (see \code{\link{gt.bpm}}
and \code{\link{LM.bpm}}). 

For recursive bivariate and sample selection models \code{\link{AT}} and \code{\link{prev}} calculate 
the average effect of an endogenous covariate and prevalence, respectively. \code{\link{mb}} provides the
nonparametric (worst-case) Manski's bound which is useful to check whether the average effect from a recursive model or prevalence from
a sample selection model is included within the possibilites of the bound.

If it makes sense then the dependence parameter of the copula distribution can be specified as a function of covariates or a grouping factor.

\code{SemiParBIVProbit} has been originally designed to deal with bivariate binary responses in the contexts of associated 
error equations, endogeneity, non-random sample selection and partial observability. However, the use of 
bivariate equation systems in which one response is binary and the other is continuous is allowed for. This case
is relevant for estimating, for instance, the effect that a binary endogenous variable has on a continuous response or viceversa. 
There are many continuous distributions which can be employed within this context and we plan to include several 
options on a case-study basis. Other bivariate models will be implemented from time to time, depending on the study at hand. 





}


\author{

Giampiero Marra (University College London, Department of Statistical Science) and Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics)

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

Thanks to Bear Braumoeller (Department of Political Science, The Ohio State University) for suggesting the implementation of bivariate models with partial observability.

}

\references{

Marra G. (2013), On P-values for Semiparametric Bivariate Probit Models. \emph{Statistical Methodology}, 10(1), 23-28.

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G., Radice R. and Missiroli S. (2014), Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Computational Statistics}, 29(3-4), 715-741.

Marra G., Radice R. and Filippou P. (in press), Regression Spline Bivariate Probit Models: A Practical Approach to Testing for Exogeneity. \emph{Communications in Statistics - Simulation and Computation}.

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.

Radice R., Marra G. and M. Wojtys (in press), Copula Regression Spline Models for Binary Outcomes. \emph{Statistics and Computing}.


}


\seealso{  \code{\link{SemiParBIVProbit}}  }


\keyword{endogeneity} 
\keyword{package}
\keyword{non-random sample selection} 
\keyword{partial observability} 
\keyword{semiparametric bivariate probit modelling} 
\keyword{smooth} 
\keyword{regression}
\keyword{variable selection}

