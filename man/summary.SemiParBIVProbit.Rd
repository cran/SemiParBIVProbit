\name{summary.SemiParBIVProbit}
\alias{summary.SemiParBIVProbit}

\title{SemiParBIVProbit summary}

\description{
It takes a fitted \code{SemiParBIVProbit} object produced by \code{SemiParBIVProbit()} and produces some summaries from it.
}

\usage{

\method{summary}{SemiParBIVProbit}(object, n.sim = 100, s.meth = "svd", prob.lev = 0.05, 
                                   thrs1 = 0.5, thrs2 = 0.5, ...)

                         
}


\arguments{ 

\item{object}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate `confidence' intervals for the association parameter and other measures (e.g., gamma measure). It should be increased if
            more precision is required.} 

\item{s.meth}{Matrix decomposition used to determine the matrix root of the covariance matrix. See the documentation of \code{mvtnorm} for further details.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 

\item{thrs1, thrs2}{Threshold to calculate the correct prediction ratios for the two binary responses. Default value is 0.5. It only works when \code{Method = "B"}.} 

\item{...}{Other arguments.}

}


\details{ 

As in the package \code{mgcv}, based on the results of Marra and Wood (2012), `Bayesian p-values' are returned for the smooth terms. These have 
better frequentist performance than their frequentist counterpart. Let \eqn{\hat{\bf f}}
and \eqn{{\bf V}_f}{V_f} denote the vector of values of a smooth term evaluated at the original covariate values and the
corresponding Bayesian covariance matrix, and let \eqn{{\bf V}_f^{r-}}{V_f^{r-}} denote 
the rank \eqn{r}{r} pseudoinverse of \eqn{{\bf V}_f}{V_f}. The statistic used 
is \eqn{T=\hat{\bf f}^\prime {\bf V}_f^{r-} \hat{\bf f}}. This is 
compared to a chi-squared distribution with degrees of freedom given by \eqn{r}{r}, which is obtained by 
biased rounding of the estimated degrees of freedom. See Marra (2013) for further details. 

Covariate selection can also be achieved using a single penalty shrinkage approach as shown in Marra and Wood (2011). 

Posterior simulation is used to obtain `confidence' intervals of nonlinear functions of parameters, such as the odds ratio and 
gamma measure discussed by Tajar et al. (2001). 

For more details on model summaries see Radice, Marra and Wojtys (submitted). 

}



\value{

\item{tableP1}{Table containing parametric estimates, their standard errors, z-values and p-values for equation 1.}
\item{tableP2}{As above but for equation 2.}
\item{tableNP1}{Table of nonparametric summaries for each smooth component including estimated degrees of freedom, estimated rank, 
                approximate Wald statistic for testing the null hypothesis that the smooth term is zero and 
                corresponding p-value, for equation 1.}
\item{tableNP2}{As above but for equation 2.}
\item{n}{Sample size.}
\item{rho}{Estimated correlation parameter between the two equations.} 
\item{theta}{Estimated copula parameter linking the two equations.}
\item{formula1, formula2}{Formulas used for equations 1 and 2.}
\item{l.sc1, l.sc2}{Number of smooth components in equations 1 and 2.}
\item{t.edf}{Total degrees of freedom of the estimated bivariate probit model.}
\item{CIrs}{`Confidence' interval for either \eqn{\rho}{\rho} or \eqn{\theta}{\theta}, depending of the model fitted.}
\item{CIl1, CIl2}{`Confidence' intervals for the shape parameters when an asymmetric link approach is employed.}
\item{sel}{This is used for internal calculations.}
\item{n.sel}{Number of selected observations in the sample selection case.}
\item{BivD, nu, PL}{These are used for internal calculations.}
\item{good}{Indicator variable indicating the observations actually used in model fitting.}
\item{OR, CIor}{Odds ratio and related CI.The odds ratio is a measure of association between binary random variables and is defined as 
                p00p11/p10p01. In the case of independence this ratio is equal to 1. It can take values in the range (-Inf, Inf) and 
                it does not depend on the marginal probabilities (Tajar et al., 2001). Interval is calculated using posterior simulation.}
\item{GM, CIgm}{Gamma measure and related CI. This measure of association was proposed by Goodman and Kruskal (1954). It is defined as 
                (OR - 1)/(OR + 1), can take values in the range (-1, 1) and does not depend on the marginal probabilities.
                Interval is calculated using posterior simulation. }
\item{xi1, xi2}{Power or shape parameters of the link functions of the two equations if an asymmetric link approach is used.}
\item{table.R}{Joint absolute frequency table of observed binary value combinations. This as well as all remaining quantities below are not provided when \code{Method = "BSS"} or \code{Method = "BPO"}.}
\item{table.P}{Joint absolute frequency table of predicted binary value combinations. }
\item{table.F}{It is given by \code{table.P}/\code{table.R}.}
\item{MR}{Based on the comparison between \code{table.P} and \code{table.R}, it reports the matching rate as percentage.}
\item{P1, P2}{Marginal probabilities for the two binary outcomes y_1 = 1 and y_2 = 1.}
\item{QPS1, QPS2}{Quadratic probability scores for the two binary responses as suggested by Diebold and Rudebusch (1989). Values are on the interval [0,2], with 0 indicating a perfect fit.}
\item{CR1, CR2}{Correct prediction ratios for the two binary responses as percentage.}

}

\references{

Diebold F.X. and Rudebusch G.D. (1989), Scoring the Leading Indicators. \emph{Journal of Business}, 62(3), 369-391. 

Goodman L.A. and Kruskal W.H. (1954), Measures of Association for Cross Classification. \emph{Journal of the American Statistical Association}, 49, 732-764.

Marra G. (2013), On P-values for Semiparametric Bivariate Probit Models. \emph{Statistical Methodology}, 10(1), 23-28.

Radice R., Marra G. and M. Wojtys (submitted), Copula Regression Spline Models for Binary Outcomes.

Marra G. and Wood S.N. (2011), Practical Variable Selection for Generalized Additive Models. \emph{Computational Statistics and Data Analysis}, 55(7), 2372-2387.

Marra G. and Wood S.N. (2012), Coverage Properties of Confidence Intervals for Generalized Additive Model Components. \emph{Scandinavian Journal of Statistics}, 39(1), 53-74.

Tajar M., Denuit M. and Lambert P. (2001), Copula-Type Representation for Random Couples with Bernoulli Margins. Discussion Papaer 0118, Universite Catholique De Louvain.


}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{AT}}, \code{\link{est.prev}}, \code{\link{SemiParBIVProbitObject}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}
}

\examples{
## see examples for SemiParBIVProbit
}


\keyword{information criteria}
\keyword{smooth}
\keyword{regression}
\keyword{shrinkage smoother}
\keyword{variable selection}






