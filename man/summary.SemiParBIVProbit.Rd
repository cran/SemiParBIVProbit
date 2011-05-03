\name{summary.SemiParBIVProbit}
\alias{summary.SemiParBIVProbit}

\title{SemiParBIVProbit summary}

\description{
It takes a fitted \code{SemiParBIVProbit} object produced by \code{SemiParBIVProbit()} and produces some summaries from it.
}

\usage{

\method{summary}{SemiParBIVProbit}(object,n.sim=1000,s.meth="svd",sig.lev=0.05,thrs1=0.5,thrs2=0.5,...)


}


\arguments{ 

\item{object}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate `confidence' intervals for \eqn{\rho}{\rho}.} 

\item{s.meth}{Matrix decomposition used to determine the matrix root of the covariance matrix. See the documentation of \code{mvtnorm} for further details.} 

\item{sig.lev}{Significance level used for intervals obtained via posterior simulation.} 

\item{thrs1,thrs2}{Threshold to calculate the correct prediction ratios for the two binary responses. Default value is 0.5.} 

\item{...}{Other arguments.}

}


\details{ 

The frequentist approximation for p-values of smooth terms described, for instance, in Marra and Wood (2011) is used. If \eqn{{\bf c}_i}{c_i}  
is the parameter vector for the ith smooth term and \eqn{{\bf V}_i}{V_i} its corresponding covariance matrix, then the 
statistic is \eqn{{\bf c}_i^\prime {\bf V}_i^{k-} {\bf c}_i}{c_i'V_i^{k-}c_i} where \eqn{{\bf V}^{k-}_i}{V_i^{k-}} is the \eqn{k}{k} 
pseudo-inverse of \eqn{{\bf V_i}}{V_i} and \eqn{k}{k} is the estimated rank of \eqn{{\bf V_i}}{V_i}. \eqn{k}{k} is determined heuristically as follows. It is the minimum
value between the maximum \code{edf} value allowed for the ith smooth term and the smallest integer not less than the quantity calculated
as 2*\code{edf}. P-values are obtained by comparing the chi.sq statistic to the chi-squared distribution with \eqn{k}{k} degrees of freedom. 

Notice that p-values are approximate. Simulation evidence suggests that if the data are not too noisy, then this p-value definition can lead to reliable
covariate selection, otherwise shrinkage approaches (Marra and Wood, 2011) are more reliable. 

}



\value{

\item{tableP1}{It returns a table containing parametric estimates, their standard errors, z-values and p-values for equation 1.}
\item{tableP2}{As above but for equation 2.}
\item{tableNP1}{It returns a table of nonparametric summaries for each smooth component including estimated degrees of freedom, estimated rank, 
                approximate Wald statistic for testing the null hypothesis that the smooth term is zero, and p-value for equation 1.}
\item{tableNP2}{As above but for equation 2.}
\item{n}{Sample size.}
\item{rho}{Estimated correlation parameter between the two equations.} 
\item{formula1}{The original GAM formula used for equation1.}
\item{formula2}{The original GAM formula used for equation2.}
\item{l.sp1}{Number of smooth components in equation 1.}
\item{l.sp2}{Number of smooth components in equation 2.}
\item{t.edf}{Total degrees of freedom of the estimated bivariate probit model.}
\item{CIrs}{It returns `confidence' intervals for \eqn{\rho}{\rho}.}
\item{sel}{This is used for internal calculations.}
\item{table.R}{Joint absolute frequency table of observed binary value combinations. This as well as all remaining quantities below are not provided in the sample selection case.}
\item{table.P}{Joint absolute frequency table of predicted binary value combinations. }
\item{table.F}{It is given by \code{table.P}/\code{table.R}.}
\item{MR}{Based on the comparison between \code{table.P} and \code{table.R}, it reports the matching rate as percentage.}
\item{P1,P2}{Marginal probabilities for the two binary outcomes y1=1 and y2=1.}
\item{QPS1,QPS2}{Quadratic probability scores for the two binary responses as suggested by Diebold and Rudebusch (1989). Values are on the interval [0,2], with 0 indicating a perfect fit.}
\item{CR1,CR2}{Correct prediction ratios for the two binary responses as percentage.}

}

\references{

Diebold F.X. and Rudebusch G.D. (1989), Scoring the Leading Indicators. \emph{Journal of Business}, 62, 369-391. 

Marra G. and Wood S.N. (2011), Practical Variable Selection for Generalized Additive Models. \emph{Computational Statistics and Data Analysis}, 55, 2372-2387.

}



\author{
Giampiero Marra, University College London

Rosalba Radice, London School of Hygiene and Tropical Medicine

}



\seealso{
\code{\link{AT}}, \code{\link{InfCr}}, \code{\link{SemiParBIVProbitObject}}, \code{\link{plot.SemiParBIVProbit}}
}

\examples{
## see examples for SemiParBIVProbit
}


\keyword{models}
\keyword{information criteria}
\keyword{smooth}
\keyword{regression}
\keyword{shrinkage smoothers}
\keyword{variable selection}






