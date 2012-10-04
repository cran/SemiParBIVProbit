\name{SemiParBIVProbitObject}
\alias{SemiParBIVProbitObject}

\title{Fitted SemiParBIVProbit object}

\description{A fitted semiparametric bivariate probit object returned by function \code{SemiParBIVProbit} and of class.}


\value{ 


\item{fit}{A list of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$argument} and \code{fit$S.h} return
           the estimated parameters and overall penalty matrix scaled by its smoothing parameters, for the bivariate probit model. See the 
           documentation of \code{trust} for diagnostics.}

\item{gam1}{A univariate GAM object for equation 1. See the documentation of \code{mgcv} for full details.}

\item{gam2}{A univariate GAM object for equation 2.}

\item{gam2.1}{A univariate GAM object for equation 2, estimated using a Heckman-type sample selection correction procedure (\code{selection=TRUE}).}

\item{p.weights}{Prior weights used during model fitting.}

\item{sp}{Estimated smoothing parameters of the smooth components for the fitted bivariate probit model.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{iter.if}{Number of iterations performed in the initial step of the Fisher scoring algorithm.}

\item{rho}{Estimated correlation coefficient between the two equations.}

\item{n}{Sample size.}

\item{n.sel}{Number of selected observations in the sample selection case.}

\item{X1}{It returns the design matrix associated with the first linear predictor.}

\item{X2}{It returns the design matrix associated with the second linear predictor.} 

\item{X1.d2}{Number of columns of the design matrix for equation 1. This is used for internal calculations.}

\item{X2.d2}{Number of columns of the design matrix for equation 2.} 

\item{l.sp1}{Number of smooth components in equation 1.}

\item{l.sp2}{Number of smooth components in equation 2.}

\item{He}{Penalized hessian.}

\item{HeSh}{Unpenalized hessian.}

\item{Vb}{Inverse of the penalized hessian. This corresponds to the Bayesian variance-covariance matrix used for `confidence' interval calculations.}

\item{F}{This is given by \code{Vb*HeSh}.} 

\item{t.edf}{Total degrees of freedom of the estimated bivariate probit model. It is calculated as \code{sum(diag(F))}.}

\item{bs.mgfit}{A list of values and diagnostics extracted from \code{magic}.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm converged.}
                
\item{wor.c}{It contains the working model quantities given by the square root of 
             the weight matrix times the pseudo-data vector and design matrix, \code{rW.Z} and \code{rW.X}.}                

\item{p11,p10,p01,p00}{Distribution function of a bivariate normal with zero means, unit variances and correlation \code{rho} evaluated 
                       at (y1=1,y2=1), (y1=1,y2=0), (y1=0,y2=1) and (y1=0,y2=0).}

\item{eta1,eta2}{The estimated linear predictors for the two equations.}

\item{dat}{It returns the full design matrix associated with the observed binary variables and two linear predictors.}

\item{sel}{This is used for internal calculations.}

\item{masses}{Probabilities associated with the \code{K} biavariate mass points.}

\item{K}{Number of bivariate mass points.}

\item{iter.npRE}{Number of iterations performed for a starting value procedure when using nonparametric random effect.}

\item{iter.npREsp}{Number of cumulative extra-iterations performed in the 
smoothing parameter estimation loop when using nonparametric random effects.} 

\item{npRE}{This is used for internal calculations.}

\item{logL}{It returns the value of the (unpenalized) log-likelihood evaluated at the (penalized or unpenalized) parameter 
            estimates.}

\item{eb.u1,eb.u2}{Estimated nonparametric random effects for each identifier.}

\item{Eb.u1,Eb.u2}{Estimated nonparametric random effects for each observation.}

\item{id}{It represents the individual identifier.}

\item{uidf}{It returns the number of observations within each group.}

\item{T.sv}{This is used for internal calculations.}


}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{SemiParBIVProbit}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}
}











