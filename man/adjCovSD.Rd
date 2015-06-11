\name{adjCovSD}
\alias{adjCovSD}

\title{Adjustment for the covariance matrix from a SemiParBIVProbit model fitted to complex survey data.}

\description{ 

\code{adjCovSD} can be used to adjust the covariance matrix of a fitted \code{SemiParBIVProbit} object.}


\usage{

adjCovSD(x, design)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{design}{A \code{svydesign} object as produced by \code{svydesign()} from the \code{survey} package.} 




}


\details{

This function has been extracted from the \code{survey} package and adapted to the class of bivariate models. It computes the sandwich 
variance estimator for a bivariate model fitted to data from a complex sample survey (Lumley, 2004).






}



\value{

This function returns a fitted \code{SemiParBIVProbit} object which is identical to that supplied in \code{adjCovSD} but with adjusted 
covariance matrix.


}



\section{WARNINGS}{

This correction may not be appropriate for models that require the use of penalties in fitting.

} 


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\references{

Lumley T. (2004), Analysis of Complex Survey Samples. \emph{Journal of Statistical Software}, 9(8), 1-19.


}




\seealso{
\code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}  
}


\keyword{covariance matrix adjustment} 
\keyword{complex survey design} 










