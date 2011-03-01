\name{SemiParBIVProbitObject}
\alias{SemiParBIVProbitObject}

\title{Fitted SemiParBIVProbit object}

\description{A fitted semiparametric bivariate probit object returned by function \code{SemiParBIVProbit} and of class.}


\value{ 

\item{fit}{A list of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$argument} and \code{fit$S.h} return
           the estimated parameters and overall penalty matrix, scaled by its smoothing parameters, for the bivariate probit model. See the 
           documentation of \code{trust} for diagnostics.}

\item{gam1}{A univariate GAM object for equation 1. This is useful to compare univariate and bivariate estimates. See the documentation of 
           \code{mgcv} for full details.}

\item{gam2}{A univariate GAM object for equation 2.}

\item{sp}{Estimated smoothing parameters of the smooth components for the fitted bivariate probit model.}

\item{rho}{Estimated correlation coefficient between the two equations.}

\item{n}{Sample size.}

\item{X1}{It returns the design matrix associated with the first linear predictor.}

\item{X2}{It returns the design matrix associated with the second linear predictor.} 

\item{X1.d2}{Number of columns of the design matrix for equation 1. This is used for internal calculations.}

\item{X2.d2}{Number of columns of the design matrix for equation 2.} 

\item{l.sp1}{Number of smooth components of equation 1.}

\item{l.sp2}{Number of smooth components of equation 2.}

\item{n.e}{Number of equations. This is used for internal calculations.}

\item{He}{Penalized hessian.}

\item{HeSh}{Unpenalized hessian.}

\item{Vb}{Inverse of the penalized hessian. This corresponds to the Bayesian variance-covariance matrix used for `confidence' interval calculations.}

\item{F}{This is given by \code{Vb*HeSh}.} 

\item{t.edf}{Total degrees of freedom of the estimated bivariate probit model. It is calculated as \code{sum(diag(F))}.}

\item{bs.mgfit}{A list of values and diagnostics extracted from \code{magic}. See the documentation 
                of \code{magic} in \code{mgcv} for full details.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm converged.}
                
\item{wor.c}{It contains the working model quantities given by the square root of 
             the weight matrix times the pseudo-data vector and design matrix, \code{rW.Z} and \code{rW.X}.}                

\item{p11,p10,p01,p00}{Distribution function of a bivariate normal with zero means, unit variances and correlation \code{rho} evaluated 
                       at (y1=1,y2=1), (y1=1,y2=0), (y1=0,y2=1) and (y1=0,y2=0).}

\item{eta1,eta2}{The estimated linear predictors for the two equations.}

\item{dat}{It returns the full design matrix associated with the observed binary variables and the two linear predictors.}

}


\author{
Giampiero Marra, University College London

Rosalba Radice, London School of Hygiene and Tropical Medicine

}






\seealso{
\code{\link{SemiParBIVProbit}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}
}











