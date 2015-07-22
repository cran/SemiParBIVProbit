\name{prev}
\alias{prev}

\title{Estimated overall prevalence from sample selection model}

\description{ 
\code{prev} can be used to calculate the overall estimated prevalence from a sample selection model, with corresponding interval
obtained using the delta method or posterior simulation.}



\usage{

prev(x, sw = NULL, type = "bivariate", ind = NULL, delta = FALSE,  
     n.sim = 100, prob.lev = 0.05, hd.plot = FALSE, 
     main = "Histogram and Kernel Density of Simulated Prevalences", 
     xlab = "Simulated Prevalences", ...)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{sw}{Survey weights.}  

\item{type}{This argument can take three values: \code{"naive"} (the prevalence is calculated ignoring the presence of observed and unobserved confounders), \code{"univariate"} (the effect is obtained from the univariate probit model which neglects the presence of unobserved confounders) and \code{"bivariate"} (the effect is obtained from the bivariate model which accounts for observed and unobserved confounders).}
             
\item{ind}{Binary logical variable. It can be used to calculate the prevalence for a subset of the data.} 

\item{delta}{If \code{TRUE} then the delta method is used for confidence interval calculations, otherwise Bayesian posterior 
             simulation is employed.} 
             
\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.}  
             
\item{prob.lev}{Overall probability of the left and right tails of the prevalence distribution used for interval calculations.}
              

\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated prevalences is produced. This can only 
               be produced when \code{delta = FALSE}.}

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}

}


\details{

\code{prev} estimates the overall prevalence of a disease (e.g., HIV) when there are missing values that are not at random. 
An interval for the estimated prevalence can be obtained using the delta method or posterior simulation. 


}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated prevalence and upper confidence interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.prev}{If \code{delta = FALSE} then it returns a vector containing simulated values of the prevalence. This 
               is used to calculate an interval.}

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.



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
\keyword{delta method}
\keyword{non-random sample selection} 









