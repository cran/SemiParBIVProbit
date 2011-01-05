\name{ATE}
\alias{ATE}

\title{Average Treatment Effect}

\description{ 
\code{ATE} can be used to calculate the average treatment effect of a binary predictor/treatment, with 
corresponding `confidence' intervals calculated via posterior simulation.}



\usage{

ATE(x,eq,nm.end,sig.lev=0.05,n.sim=1000,s.meth="svd")

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{eq}{The equation containing the binary predictor of interest.} 

\item{nm.end}{The name of the binary variable.} 

\item{sig.lev}{Significance level.} 

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters.} 

\item{s.meth}{Matrix decomposition used to determine the matrix root of the covariance matrix. See the documentation of the \code{mvtnorm} package for further details.} 


}


\details{

ATE measures the impact that a binary predictor/treatment has on the probability of observing 1 for a binary response, conditional 
on the other covariates. Notice that in order to calculate this effect, all covariates, except for the treatment, are taken to 
be equal to their modes. The corresponding `confidence' intervals are calculated via simulation from the posterior distribution 
of the estimated model parameters.  

}


\value{
The function returns the ATE with corresponding `confidence' intervals.
}


\references{

Marra G. and Radice R. (in press), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. The Canadian 
Journal of Statistics.

}

\author{ 

Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

}


\seealso{
\code{\link{InfCr}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}  
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit models} 
\keyword{average treatment effect} 
\keyword{bayesian posterior simulation}








