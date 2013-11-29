\name{est.prev}
\alias{est.prev}

\title{Estimated overall prevalence from sample selection models}

\description{ 
\code{est.prev} can be used to calculate the overall estimated prevalence from a sample selection model, with corresponding confidence intervals.}



\usage{

est.prev(x, sig.lev=0.05)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{sig.lev}{Significance level.} 


}


\details{

est.prev estimates the overall prevalence of a disease (e.g. HIV) when there are missing values that are not at random. 

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated prevalence and upper confidence interval limit.}
\item{sig.lev}{Significance level used.}

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\seealso{
\code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}  
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{prevalence} 









