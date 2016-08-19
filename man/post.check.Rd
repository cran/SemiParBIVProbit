\name{post.check}
\alias{post.check}

\title{Diagnostic plots for discrete/continuous response margin}

\description{ 

It produces diagnostic plots based on (randomised) quantile residuals. 

}

\usage{

post.check(x, main = "Histogram and Density Estimate of Residuals", 
           main2 = "Histogram and Density Estimate of Residuals",
           xlab = "Quantile Residuals", xlab2 = "Quantile Residuals", 
           test = FALSE, ...)

}

\arguments{ 

\item{x}{A fitted \code{copulaReg}/\code{copulaSampleSel} object.}

\item{main}{Title for the plot.}

\item{main2}{Title for the plot in the second row. This comes into play only when fitting models with two non-binary margins.}

\item{xlab}{Title for the x axis.}

\item{xlab2}{Title for the x axis in the second row. As above.}

\item{test}{If \code{TRUE} then the shapiro test for normality is performed on normalised quantile residuals.}

\item{...}{Other graphics parameters to pass on to plotting commands.}

}


\details{ 

If the model fits the response well then the plots should look normally distributed.
When fitting models with discrete and/or continuous margins, four plots will be produced. In this case,
the arguments \code{main2} and \code{xlab2} come into play and allow for different
labelling across the plots. 
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{copulaReg}}, \code{\link{copulaSampleSel}}  
}


\keyword{flexible copula regression modelling} 
\keyword{regression} 
\keyword{Q-Q plot} 
\keyword{density plot}
\keyword{histogram}









