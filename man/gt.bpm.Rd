\name{gt.bpm}
\alias{gt.bpm}

\title{Gradient test}

\description{ 
\code{gt.bpm} can be used to test the hypothesis of absence of endogeneity or correlated model equations/errors or non-random sample selection.}


\usage{

gt.bpm(object)

}

\arguments{ 

\item{object}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.} 

}


\details{

The gradient test was first proposed by Terrell (2002) and it is based on classic likelihood theory results. See Marra et al. (submitted) for full details. 



}


\value{

It returns a numeric p-value corresponding to the null hypothesis \eqn{\rho}{\rho=0}.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{


Marra G., Radice R. and Filippou P. (submitted), Testing the Hypothesis of Exogeneity in Regression Spline Bivariate Probit Models.

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








