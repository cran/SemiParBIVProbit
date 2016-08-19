\name{SemiParBIVProbit-package}
\alias{SemiParBIVProbit-package}
\docType{package}

\title{Semiparametric Copula Regression Models}

\description{

This package provides functions for fitting several classes of copula regression models with several types of covariate 
effects. Many copula and marginal distributions are supported. The copula dependence parameter can 
also be specified as a flexible function of covariates.

This package, and in particular its main function \code{SemiParBIVProbit()}, has been originally designed 
to deal with bivariate binary responses and in fact the first model
introduced was the semiparametric bivariate probit model (hence the name of the package). However, 
since then more models and options have been introduced.

The main fitting functions are listed below.

\code{SemiParBIVProbit()} fits bivariate regression models with binary responses (where the link functions
are not restricted to be just probit). This is useful to fit bivariate binary models in the presence of 
(i) non-random sample selection or (ii) associated responses/endogeneity or (iii) partial observability. This function includes
the \code{Model} argument which allows the user to fit the model in one of three situations mentioned above.

\code{copulaReg()} fits bivariate models with binary/discrete/continuous margins in the presence of 
associated responses/endogeneity. 

\code{copulaSampleSel()} fits bivariate sample selection models with continuous/discrete response (instead of 
binary as it would be the case when using \code{SemiParBIVProbit()}).  

\code{SemiParTRIV()} fits trivariate probit models (with and without double sample selection).  

Other models/options will be incorporated from time to time.


}

\details{
\code{SemiParBIVProbit} provides functions for fitting flexible copula regression models in various situations. The underlying 
representation and 
estimation of the modelling framework is based on a penalized regression spline approach, with automatic 
smoothness selection. Several marginal and copula distributions are available. The 
numerical routine carries out function minimization using a trust region algorithm in combination with
an adaptation of a smoothness estimation fitting procedure for GAMs (see \code{mgcv} for more details on this last point). 
 
Many smoothers are supported and are extracted from \code{mgcv}. Estimation is by penalized 
maximum likelihood with automatic smoothness estimation achieved by using an approximate AIC. 

Confidence intervals for smooth components and nonlinear functions of the model
parameters are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided and based on the approach
implemented in \code{mgcv}. The usual plotting and summary functions are also available. Model/variable 
selection is also possible via the use of shrinakge smoothers and/or information criteria. 

The dependence parameter of the copula distribution can be specified as a function of covariates or a grouping factor if it makes sense. More
generally, it is possible to specify all marginal distribution and copula parameters as functions of
covariates, as it would be typical in GAMLSS.


}


\author{

Giampiero Marra (University College London, Department of Statistical Science) and Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics)

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

Thanks to Bear Braumoeller (Department of Political Science, The Ohio State University) for suggesting the implementation of bivariate models with partial observability.


}

\references{

Key references:

Filippou P., Marra G. and Radice R. (submitted), A Penalised Estimation Approach to Additive Trivariate Probit Regression.

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G. and Radice R. (submitted), Bivariate Copula Additive Models for Location, Scale and Shape.

Marra G., Radice R., Barnighausen T., Wood S.N. and McGovern M.E. (in press), A Simultaneous Equation Approach to Estimating HIV Prevalence with Non-Ignorable Missing Responses. \emph{Journal of the American Statistical Association}.

Marra G. and Wyszynski K. (in press), Semi-Parametric Copula Sample Selection Models for Count Responses. \emph{Computational Statistics and Data Analysis}.

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.

Radice R., Marra G. and Wojtys M. (2016), Copula Regression Spline Models for Binary Outcomes. \emph{Statistics and Computing}, 26(5), 981-995. 

Wojtys M. and Marra G. (submitted). Copula based generalized additive models with non-random sample selection.

}


\seealso{\code{\link{SemiParBIVProbit}}, \code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{SemiParTRIV}}}


\keyword{endogeneity} 
\keyword{package}
\keyword{non-random sample selection} 
\keyword{partial observability} 
\keyword{flexible copula regression modelling} 
\keyword{smooth} 
\keyword{regression}


