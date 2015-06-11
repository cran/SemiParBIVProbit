\name{AT}
\alias{AT}

\title{Average treatment effect of a binary endogenous variable}

\description{ 
\code{AT} can be used to calculate the sample average treatment effect of a binary endogenous predictor/treatment, with 
corresponding interval obtained using the delta method or posterior simulation.}



\usage{

AT(x, eq, nm.bin, E = TRUE, treat = TRUE, naive = FALSE, ind = NULL, 
   sub.l = 50, delta = FALSE, n.sim = 100, prob.lev = 0.05, hd.plot = FALSE,
   main = "Histogram and Kernel Density of Simulated Average Effects", 
   xlab = "Simulated Average Effects", ...)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{eq}{Equation containing the binary endogenous predictor of interest.} 

\item{nm.bin}{Name of the binary endogenous variable.} 

\item{E}{If \code{TRUE} then \code{AT} calculates the sample ATE. If \code{FALSE} then it calculates the sample AT 
         for the treated individuals only.}  

\item{treat}{If \code{TRUE} then \code{AT} calculates the AT using the treated only. If \code{FALSE} then it calculates the effect on 
the control group. This only makes sense if \code{E = FALSE}.}  

\item{naive}{It indicates whether AT is calculated from a univariate probit model. This option has been 
             introduced to compare adjusted (for unobserved confounding) and unadjusted estimates. Note that, although
             the unadjusted estimates do not account for unobserved confounding, they do control for observed confounding.} 

\item{ind}{Binary logical variable. It can be used to calculate the AT for a subset of the data. Note that it does not make sense to use \code{ind} 
           when some observations are excluded from the AT calculation (e.g., when using \code{E = FALSE}).}  
  
\item{sub.l}{Length of subset of observations to use in the calculation of the treatment effect when the outcome of interest is continuous.}  

\item{delta}{If \code{TRUE} then the delta method is used for confidence interval calculations, otherwise Bayesian posterior 
             simulation (the most reliable option, despite a bit slower) is employed.}  

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the AT distribution used for interval calculations.} 

\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated average effects is produced. This can 
               only be produced when \code{delta = FALSE}.}

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

For the case of continuous response, the AT is calculated via numerical integration, which can result in costly computation depending on the
size of the dataset and number of simulated coefficient vectors used for interval construction. In the model building 
process, we recommend using a subset of observations to calculate the AT (see \code{sub.l}). In 
our experience this still provides representative average effects. Once a preferred model has been found, the AT can be calculated using 
the entire  dataset by setting \code{sub.l} to the number of observations (note that this can be time consuming). 

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated AT and upper interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.AT}{If \code{delta = FALSE} then it returns a vector containing simulated values of the average treatment effect. This 
               is used to calculate intervals.}
\item{AT.so}{It returns a vector containing the estimated effect for each single observation.}            
               

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








