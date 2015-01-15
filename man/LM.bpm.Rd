\name{LM.bpm}
\alias{LM.bpm}

\title{Lagrange Multiplier Test (Score Test)}

\description{ 
Before fitting a bivariate probit model, \code{LM.bpm} can be used to test the hypothesis of absence of endogeneity, correlated model equations/errors 
or non-random sample selection.}


\usage{

LM.bpm(formula, data = list(), weights = NULL, subset = NULL, Model, 
       hess = TRUE, infl.fac = 1, pPen1 = NULL, pPen2 = NULL)

}

\arguments{ 


\item{formula}{A list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms are used to specify 
               smooth smooth functions of predictors. Note that if \code{Model = "BSS"} then the first formula MUST refer 
               to the selection equation.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{Model}{It indicates the type of model to be used in the analysis. Possible values are "B" (bivariate model) and 
             "BSS" (bivariate model with sample selection).}

\item{hess}{If \code{FALSE} then the expected (rather than observed) information matrix is employed.}

\item{infl.fac}{Inflation factor for the model degrees of freedom in the UBRE score. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{pPen1, pPen2}{Optional list specifying any penalties to be applied to the parametric model terms of equations 1 and 2.}


}


\details{

This Lagrange multiplier test (also known as score test) is used here for testing the null 
hypothesis that \eqn{\rho}{\rho} is equal to 0 (i.e. no endogeneity, non-random sample selection or correlated model equations/errors, depending 
on  the model being fitted). Its main advantage is that it does 
not require an estimate of the model parameter vector under the alternative hypothesis. Asymptotically, it takes a Chi-squared distribution 
with one degree of freedom. Full details can be found in Marra et al. (2014) and Marra et al. (submitted).


}


\value{

It returns a numeric p-value corresponding to the null hypothesis that \eqn{\rho}{\rho} is equal to 0.

}


\section{WARNINGS}{

This test's implementation is only valid for models with normal errors and probit links.

} 




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G., Radice R. and Filippou P. (submitted), Regression Spline Bivariate Probit Models: A Practical Approach to Testing for Exogeneity.

Marra G., Radice R. and Missiroli S. (2014), Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Computational Statistics}, 29(3-4), 715-741.

}




\seealso{
\code{\link{SemiParBIVProbit}}
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{lagrange multiplier test} 
\keyword{score test} 
\keyword{endogeneity}
\keyword{non-random sample selection}
\keyword{correlated model equations/errors}








