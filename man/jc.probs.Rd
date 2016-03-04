\name{jc.probs}
\alias{jc.probs}

\title{Joint or conditional probabilities from fitted bivariate model}

\description{ 
\code{jc.probs} can be used to calculate the joint or conditional probabilities from a fitted bivariate model with intervals obtained 
using posterior simulation.}



\usage{

jc.probs(x, y1, y2, type = "bivariate", cond = 0,
         intervals = FALSE, n.sim = 100, prob.lev = 0.05)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} or \code{copulaReg} object as produced by \code{SemiParBIVProbit()} or \code{copulaReg()}.}

\item{y1}{Value of response for first margin.}

\item{y2}{Value of response for second margin.}

\item{type}{This argument can take two: \code{"bivariate"} (the probabilities are calculated from the fitted 
            bivariate model) and \code{"independence"} (the calculation is done from univariate fits).}

\item{cond}{There are three possible values: 0 (joint probabilities are delivered), 1 (conditional probabilities are delivered and 
            conditioning is with the respect to the first margin), 2 (as before but conditioning is with the respect to 
            the second margin).}  
  
\item{intervals}{If \code{TRUE} then intervals for the probabilities are also produced.}   
  
\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used
             for interval calculations.} 

\item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.} 
  
  
}


\details{

This function calculates joint or conditional probabilities from a fitted bivariate model or a model assuming independence, with intervals 
obtained using posterior simulation.

}



\value{

\item{res}{It returns three values: estimated probabilities (\code{p12}), with lower and upper interval limits (\code{CIpr}) 
if \code{intervals = TRUE}, and \code{p1} and \code{p2} (the marginal probabilities).}
                        
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\references{

Marra G. and Radice R. (submitted), Flexible Bivariate Copula-Based Regression Models.


}


\seealso{
\code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{copulaReg}}  
}

\examples{
## see examples for SemiParBIVProbit and copulaReg
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{copula}
\keyword{bivariate probit model} 
\keyword{bayesian posterior simulation}







