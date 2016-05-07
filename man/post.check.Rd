\name{post.check}
\alias{post.check}

\title{Diagnostic plots for continuous response margin}

\description{ 

It produces diagnostic plots based on quantile residuals. 

}

\usage{

post.check(x, main = "Histogram and Density Estimate of Residuals", 
           main2 = "Histogram and Density Estimate of Residuals",
           xlab = "Quantile Residuals", xlab2 = "Quantile Residuals", ...)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{main}{Title for the plot.}

\item{main2}{Title for the plot in the second row. This comes into play only when fitting models with two continuous margins.}

\item{xlab}{Title for the x axis.}

\item{xlab2}{Title for the x axis in the second row. As above.}

\item{...}{Other graphics parameters to pass on to plotting commands.}

}


\details{ 

If the model fits the response well then the plots should look normally distributed.
When fitting models with two continuous margins, four plots will be produced. In this case,
the arguments \code{main2} and \code{xlab2} come into play and allow for different
labelling across the plots. 
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

\keyword{semiparametric bivariate modelling} 
\keyword{semiparametric trivariate probit modelling} 
\keyword{regression} 
\keyword{Q-Q plot} 
\keyword{density plot}
\keyword{histogram}









