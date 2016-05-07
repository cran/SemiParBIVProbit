\name{SemiParBIVProbit-package}
\alias{SemiParBIVProbit-package}
\docType{package}

\title{Semiparametric Copula Bivariate Regression Models}

\description{

This package provides functions for fitting copula bivariate regression models with various types of covariate effects for several situations.  
Many bivariate copula and marginal distributions are supported. The copula dependence parameter can 
also be specified as a flexible function of covariates.

This package, and in particular its main function \code{SemiParBIVProbit()}, has been originally designed 
to deal with binary responses and, in fact, the first model
that has been introduced was the semiparametric bivariate probit model (hence the name of the package). However, 
since then we have introduced more models and options, and as a consequence fitting functions (to make the 
different models more manageable for development).

The fitting functions are listed below.

\code{SemiParBIVProbit()} allows the user to fit bivariate regression models with binary responses (where the link functions
are not restricted to be just probit). This is useful to fit bivariate binary models in the presence of 
non-random sample selection or correlated/associated responses or endogeneity or partial observability. In fact, this function includes
an argument which allows the user to employ the desired model. This function also permits to fit models in which one response
is binary and the other is continuous; this is useful to deal with correlated/associated responses or endogeneity, for instance.

\code{copulaReg()} allows the user to fit bivariate models with continuous margins in the presence of correlated/associated responses. 

\code{copulaSampleSel()} allows the user to fit bivariate sample selection models with continuous response, instead of binary as it would be the case
when using SemiParBIVProbit().  

\code{SemiParTRIVProbit()} allows the user to fit trivariate probit models (with and without double sample selection).  

Other models/options will be incorporated from time to time.


}

\details{
\code{SemiParBIVProbit} provides functions for fitting flexible copula regression models, in the presence of 
associated error equations, endogeneity, non-random sample selection or partial observability. The underlying representation and 
estimation of the model is based on a penalized likelihood-based regression spline approach, with automatic 
smoothness selection. Several
marginal and bivariate copula distributions are available. The 
numerical routine carries out function minimization using a trust region algorithm in combination with
an adaptation of a smoothness estimation fitting procedure for GAMs (see \code{mgcv} for more details on the this last point). 
 
Many smoothers are supported and are extracted from \code{mgcv}. Estimation is by penalized 
maximum likelihood with automatic smoothness estimation achieved 
by using an approximate AIC. 

Confidence intervals for smooth components and nonlinear functions of the model
parameters are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided and based on the approach
implemented in \code{mgcv}. The usual plotting and summary finctions are also available. Model/variable 
selection is also possible via the use of shrinakge smoothers and/or information criteria. 


If it makes sense then the dependence parameter of the copula distribution can be specified as a function of covariates or a grouping factor, for instance.

The models implemented also allow the user to specify all marginal distribution and copula parameters as functions of
covariates, as it would be typical in a GAMLSS approach.


}


\author{

Giampiero Marra (University College London, Department of Statistical Science) and Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics)

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

Thanks to Bear Braumoeller (Department of Political Science, The Ohio State University) for suggesting the implementation of bivariate models with partial observability.

Thanks to Panagiota Filippou for the work done on \code{SemiParTRIVProbit()}.

}

\references{

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G. and Radice R. (submitted), A Bivariate Copula Additive Model for Location, Scale and Shape.

Marra G., Radice R. and Missiroli S. (2014), Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Computational Statistics}, 29(3-4), 715-741.

Marra G., Radice R. and Filippou P. (in press), Regression Spline Bivariate Probit Models: A Practical Approach to Testing for Exogeneity. \emph{Communications in Statistics - Simulation and Computation}.

Marra G., Radice R., Barnighausen T., Wood S.N. and McGovern M.E. (submitted). A Simultaneous Equation Approach to Estimating HIV Prevalence with Non-Ignorable Missing Responses.

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.

Radice R., Marra G. and M. Wojtys (in press), Copula Regression Spline Models for Binary Outcomes. \emph{Statistics and Computing}.


}


\seealso{\code{\link{SemiParBIVProbit}}, \code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{SemiParTRIVProbit}}}


\keyword{endogeneity} 
\keyword{package}
\keyword{non-random sample selection} 
\keyword{partial observability} 
\keyword{semiparametric bivariate modelling} 
\keyword{semiparametric trivariate modelling} 
\keyword{smooth} 
\keyword{regression}
\keyword{variable selection}

