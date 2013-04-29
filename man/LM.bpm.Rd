\name{LM.bpm}
\alias{LM.bpm}

\title{Lagrange Multiplier Test (Score Test)}

\description{ 
Before fitting a bivariate probit model, \code{LM.bpm} can be used to test the hypothesis of absence of endogeneity 
or correlated model equations/errors or non-random sample selection.}


\usage{

LM.bpm(formula.eq1, formula.eq2, data, selection=FALSE, FI=FALSE)

}

\arguments{ 

\item{formula.eq1}{A formula for equation 1.} 

\item{formula.eq2}{A formula for equation 2.} 

\item{data}{An optional data frame, list or environment containing the variables in the model. If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{LM.bpm} is called.}

\item{selection}{If \code{TRUE}, then the test is performed for the sample selection model case.}

\item{FI}{If \code{TRUE}, then the Fisher (rather than the observed) information matrix is used.}

}


\details{
This Lagrange multiplier test (also known as score test) is used here for testing the null 
hypothesis \eqn{\rho}{\rho=0} (i.e. no endogeneity or non-random sample selection, depending on the model being fitted). Its main advantage is that it does 
not require an estimate of the model parameter vector under the alternative hypothesis. Asymptotically, it takes a Chi-squared distribution 
with one degree of freedom. Full details can be found in Marra et al. (submitted).

Note that, when \code{selection=FALSE} and \code{FI=TRUE}, a convenient simplification based on the result that 
the Fisher information matrix becomes block diagonal is employed (Marra et al., submitted).


}


\value{

It returns a numeric p-value corresponding to the null hypothesis \eqn{\rho}{\rho=0}.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G., Radice R. and Missiroli S., Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Submitted}.


}




\seealso{
\code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}
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








