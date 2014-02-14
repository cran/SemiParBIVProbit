\name{AT}
\alias{AT}

\title{Average treatment effect of a binary endogenous variable}

\description{ 
\code{AT} can be used to calculate the sample average treatment effect of a binary endogenous predictor/treatment, with corresponding confidence intervals.}



\usage{

AT(x, eq, nm.bin="", E=TRUE, treat=TRUE, delta=FALSE, sig.lev=0.05, 
   s.meth="svd", n.sim=1000)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{eq}{Equation containing the binary endogenous predictor of interest.} 

\item{nm.bin}{Name of the binary endogenous variable.} 

\item{E}{If \code{TRUE}, then \code{AT} calculates the sample ATE. If \code{FALSE}, then it calculates the sample AT for the treated individuals only.}  

\item{treat}{If \code{TRUE}, then \code{AT} calculates the AT using the treated only. If \code{FALSE}, then it calculates the effect on 
the control group. This only makes sense if used jointly with \code{E=FALSE}.}  

\item{delta}{If \code{TRUE} then an approximate delta method is used for confidence interval calculation, otherwise Bayesian posterior simulation (the most reliable option, despite a bit slower) 
             is employed. Note that for models including random effects in the linear predictors or power link functions only posterior simulation is allowed for.}  

\item{sig.lev}{Significance level.} 

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta=FALSE}.} 

\item{s.meth}{Matrix decomposition used to determine the matrix root of the covariance 
              matrix. This is used when \code{delta=FALSE}. See the documentation of the \code{mvtnorm} package for further details.}




}


\details{

AT measures the sample average difference in outcomes under 
the treatment (the binary predictor/treatment assumes value 1) and under the 
control (the binary treatment assumes value 0). See the references below for details. 

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated AT and upper confidence interval limit.}
\item{sig.lev}{Significance level used.}
\item{est.ATb}{When \code{delta=FALSE} it returns a vector containing simulated values of the average treatment effect. This is used to calculate intervals.}

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\references{

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian 
Journal of Statistics}, 39(2), 259-279.

Radice R., Marra G. and M. Wojtys (submitted), Copula Regression Spline Models for Binary Outcomes with Application in Health Care Utilization.


}


\seealso{
\code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}  
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{ATE} 
\keyword{average treatment effect} 
\keyword{bayesian posterior simulation}








