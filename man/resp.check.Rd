\name{resp.check}
\alias{resp.check}

\title{Plots for response variable}

\description{ 

It produces a histogram of the response along with the estimated density from the assumed distribution as well as a normal Q-Q plot for 
the normalised quantile response. 

}

\usage{

resp.check(y, margin = "N", main = "Histogram and Density of Response",
               xlab = "Response", ...)

}

\arguments{ 

\item{y}{Response.}
 
\item{margin}{The distributions allowed are: normal ("N"), log-normal ("LN"), 
              Gumbel ("GU"), reverse Gumbel ("rGU"), logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA").} 

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands.}

}


\details{

Prior to fitting a model with binary-continuous bivariate response, a distribution for the continuous response
may be chosen by looking at the histogram of the response along with the estimated density from the assumed distribution and at the
normalised quantile responses. These will provide a rough guide to the adequacy of the chosen distribution.
The latter are defined as the quantile standard normal function of the cumulative distribution function of the response with scale and location
estimated by MLE. These should behave approximately as normally distributed variables (even though the original 
observations are not). Therefore, a normal Q-Q plot is appropriate here. 

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









