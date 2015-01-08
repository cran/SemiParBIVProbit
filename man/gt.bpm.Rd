\name{gt.bpm}
\alias{gt.bpm}

\title{Gradient test}

\description{ 
\code{gt.bpm} can be used to test the hypothesis of absence of endogeneity, correlated model equations/errors or non-random sample selection.}


\usage{

gt.bpm(object)

}

\arguments{ 

\item{object}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.} 

}


\details{

The gradient test was first proposed by Terrell (2002) and it is based on classic likelihood theory. See Marra et al. (submitted) for full details. 



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

Terrell G. (2002), The Gradient Statistic. \emph{Computing Science and Statistics}, 34, 206-215.


}




\seealso{
\code{\link{SemiParBIVProbit}}
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{gradient test} 
\keyword{endogeneity}
\keyword{non-random sample selection}
\keyword{correlated model equations/errors}








