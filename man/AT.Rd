\name{AT}
\alias{AT}

\title{Average treatment effect of a binary or continuous endogenous variable}

\description{ 
\code{AT} can be used to calculate the sample average treatment effect of a binary or continuous endogenous predictor/treatment, with 
corresponding interval obtained using mainly posterior simulation.}



\usage{

AT(x, nm.end, E = TRUE, treat = TRUE, type = "bivariate", ind = NULL, 
   delta = FALSE, n.sim = 100, prob.lev = 0.05, hd.plot = FALSE,
   prob.plot = FALSE, main = "Histogram and Kernel Density of Simulated Average Effects", 
   xlab = "Simulated Average Effects", ...)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{nm.end}{Name of the endogenous variable.} 

\item{E}{If \code{TRUE} then \code{AT} calculates the sample ATE. If \code{FALSE} then it calculates the sample AT 
         for the treated individuals only.}  

\item{treat}{If \code{TRUE} then \code{AT} calculates the AT using the treated only. If \code{FALSE} then it calculates the effect on 
the control group. This only makes sense if \code{E = FALSE}.}  

\item{type}{This argument can take three values: \code{"naive"} (the effect is calculated ignoring the presence of observed and unobserved 
            confounders), \code{"univariate"} (the effect is obtained from the univariate model 
            which neglects the presence of unobserved confounders) and \code{"bivariate"} (the effect is obtained from 
            the bivariate model which accounts for observed and unobserved confounders).}



\item{ind}{Binary logical variable. It can be used to calculate the AT for a subset of the data. Note that it does not make sense to use \code{ind} 
           when some observations are excluded from the AT calculation (e.g., when using \code{E = FALSE}).}  
  
\item{delta}{If \code{TRUE} then the delta method is used for confidence interval calculations, otherwise Bayesian posterior 
             simulation (the most reliable option, despite a bit slower) is employed. This option is only available for binary treatment and 
             binary outcome.}  

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the AT distribution used for interval calculations.} 

\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated average effects is produced. This can 
               only be produced when \code{delta = FALSE}.}

\item{prob.plot}{For the case of continuous endogenous variable and binary outcome, if \code{TRUE} then a plot  
                 showing probability that the binary outcome is equal to 1 for each value of the endogenous variable 
                 and respective intervals is produced.}

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}


}


\details{

AT measures the sample average difference in outcomes under 
treatment (the binary predictor or treatment assumes value 1) and under 
control (the binary treatment assumes value 0). Posterior simulation 
and delta method can be used to obtain a confidence/credible interval. The former is more reliable, especially at small
sample sizes. See the references below for details. 

AT can also calculate the effect that a continuous endogenous variable has on a binary outcome. In this case the effect will depend on the 
unit increment chosen (as shown by the plot produced). An overall AT is provided as well.

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated AT and upper interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.AT}{If \code{delta = FALSE} then it returns a vector containing simulated values of the average treatment effect. This 
               is used to calculate intervals.}
\item{AT.so}{It returns a vector containing the estimated effect for each single observation. This may not be available in some cases.}            
\item{Effects}{For the case of continuous endogenous variable and binary outcome, it returns a matrix made up of 
               three columns containing the effects for each unit increase in the endogenous variable and respective intervals.}            
\item{Pr}{For the case of continuous endogenous variable and binary outcome, it returns a matrix made up of 
               three columns containing the probability that the binary outcome is equal to 1 for each value of the endogenous variable 
               and respective intervals.}            
               

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\references{

Radice R., Marra G. and M. Wojtys (in press), Copula Regression Spline Models for Binary Outcomes. \emph{Statistics and Computing}.


}


\seealso{
\code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}  
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{ATE} 
\keyword{average treatment effect} 
\keyword{bayesian posterior simulation}
\keyword{delta method}








