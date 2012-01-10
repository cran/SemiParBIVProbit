\name{AT}
\alias{AT}

\title{Average Effect}

\description{ 
\code{AT} can be used to calculate the unconditional average treatment effect (ATE) and average treatment effect on the treated (ATT) of a binary predictor/treatment, with corresponding `confidence' intervals conveniently calculated via posterior simulation.}



\usage{

AT(x,eq,nm.bin="",sig.lev=0.05,n.sim=1000,s.meth="svd",E=TRUE,treat=TRUE)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{eq}{The equation containing the binary predictor of interest.} 

\item{nm.bin}{The name of the binary variable.} 

\item{sig.lev}{Significance level.} 

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters.} 

\item{s.meth}{Matrix decomposition used to determine the matrix root of the covariance matrix. See the documentation of the \code{mvtnorm} package for further details.}

\item{E}{If \code{TRUE}, then \code{AT} calculates the ATE. If \code{FALSE}, then it calculates the ATT.}  

\item{treat}{If \code{TRUE}, then \code{AT} calculates the ATT. If \code{FALSE}, then it calculates the average treatment effect on 
the control group. This only makes sense if used jointly with \code{E=FALSE}.}  


}


\details{

AT measures the average causal difference in outcomes under the treatment (the binary predictor/treatment assumes value 1) and under the control (the binary treatment assumes value 0). The corresponding `confidence' intervals are conveniently calculated via simulation from the posterior distribution 
of the estimated model parameters.  

}



\value{

\item{res}{It returns three values, namely the lower confidence interval limit, estimated AT and upper confidence interval limit, respectively.}
\item{sig.lev}{The significance level used.}

}




\author{

Maintainer: Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

}


\references{

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian 
Journal of Statistics}, 39(2), 259-279.

}


\seealso{
\code{\link{InfCr}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}  
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{ATE} 
\keyword{ATT} 
\keyword{average treatment effect} 
\keyword{average treatment effect on the treated} 
\keyword{bayesian posterior simulation}








