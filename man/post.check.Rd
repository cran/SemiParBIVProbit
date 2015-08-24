\name{post.check}
\alias{post.check}

\title{Diagnostic plots for continuous response margin}

\description{ 

It produces diagnostic plots based on quantile residuals. 

}

\usage{

post.check(x, main = "Density Estimate",
           xlab = "Quantile Residuals", ...)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands.}

}


\details{ 

If the model fits the response well then the plots should look normally distributed.
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{SemiParBIVProbit}}  
}

\examples{
## see example 5 for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{regression} 
\keyword{Q-Q plot} 
\keyword{density plot}
\keyword{histogram}









