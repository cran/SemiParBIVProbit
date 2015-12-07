\name{OR}
\alias{OR}

\title{Average causal odds ratio of a binary or continuous endogenous variable}

\description{ 
\code{OR} can be used to calculate the sample average causal odds ratio of a binary or continuous endogenous predictor/treatment, with 
corresponding interval obtained using posterior simulation.}



\usage{

OR(x, nm.end, E = TRUE, treat = TRUE, type = "bivariate", ind = NULL, 
   n.sim = 100, prob.lev = 0.05, hd.plot = FALSE,
   prob.plot = FALSE, main = "Histogram and Kernel Density of Simulated Odds Ratios", 
   xlab = "Simulated Odds Ratios", ...)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}
 
\item{nm.end}{Name of the endogenous variable.} 

\item{E}{If \code{TRUE} then \code{OR} calculates the sample OR. If \code{FALSE} then it calculates the sample OR 
         for the treated individuals only.}  

\item{treat}{If \code{TRUE} then \code{OR} calculates the OR using the treated only. If \code{FALSE} then it calculates the ratio using  
the control group. This only makes sense if \code{E = FALSE}.}  

\item{type}{This argument can take three values: \code{"naive"} (the effect is calculated ignoring the presence of observed and 
unobserved confounders), \code{"univariate"} (the effect is obtained from the univariate model which neglects 
the presence of unobserved confounders) and \code{"bivariate"} (the effect is obtained from the bivariate model which accounts for observed and unobserved confounders).}

\item{ind}{Binary logical variable. It can be used to calculate the OR for a subset of the data. Note that it does not make sense to use \code{ind} 
           when some observations are excluded from the OR calculation (e.g., when using \code{E = FALSE}).}  

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the OR distribution used for interval calculations.} 

\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated odds ratios is produced. This can 
               only be produced when \code{delta = FALSE}.}

\item{prob.plot}{For the case of continuous endogenous variable and binary outcome, if \code{TRUE} then a plot  
                 showing probability that the binary outcome is equal to 1 for each value of the endogenous variable 
                 and respective intervals is produced.}

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}


}


\details{

OR calculates the sample average causal odds ratio for a binary or continuous treatment. Posterior simulation is used to obtain a 
confidence/credible interval. 


}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated OR and upper interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.OR}{It returns a vector containing simulated values of the average OR. This 
               is used to calculate intervals.}
\item{OR.so}{It returns a vector containing the estimated effect for each single observation. This may not be available in some cases.}            
\item{Ratios}{For the case of continuous endogenous treatment and binary outcome, it returns a matrix made up of 
               three columns containing the odds ratios for each unit increase in the endogenous variable and respective intervals.}            
\item{Pr}{For the case of continuous endogenous exposure and binary outcome, it returns a matrix made up of 
               three columns containing the probability that the binary outcome is equal to 1 for each value of the endogenous variable, 
               and respective intervals.}               

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
\keyword{OR} 
\keyword{odds ratio} 
\keyword{bayesian posterior simulation}








