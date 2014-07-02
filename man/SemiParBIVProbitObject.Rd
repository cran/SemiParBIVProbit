\name{SemiParBIVProbitObject}
\alias{SemiParBIVProbitObject}

\title{Fitted SemiParBIVProbit object}

\description{A fitted semiparametric bivariate probit object returned by function \code{SemiParBIVProbit} and of class.}


\value{ 


\item{fit}{List of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$gradient} and \code{fit$S.h} return
           the gradient vector and overall penalty matrix scaled by its smoothing parameters, for the bivariate probit model. See the 
           documentation of \code{trust} for details on the diagnostics provided.}

\item{gam1}{Univariate fit for equation 1. See the documentation of \code{mgcv} for full details.}

\item{gam2}{Univariate fit for equation 2.}

\item{gam2.1}{Univariate fit for equation 2, estimated using an adaptation of the Heckman sample selection correction procedure when \code{selection = TRUE}.}

\item{coefficients}{The coefficients of the fitted semiparametric sample selection model. They are given in the following 
                    order: parametric and regression spline (if present)
                    coefficients for the first equation, parametric and regression spline
                    coefficients for the second equation, and association coefficient.}

\item{weights}{Prior weights used during model fitting.}

\item{sp}{Estimated smoothing parameters of the smooth components for the fitted bivariate probit model.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{iter.if}{Number of iterations performed in the initial step of the algorithm.}

\item{iter.inner}{Number of iterations performed inside smoothing parameter estimation step.}

\item{rho}{Estimated correlation coefficient between the two equations, when a bivariate normal model is employed.}

\item{theta}{Estimated copula parameter linking the two equations.}

\item{delta}{This is an additional estimated copula parameter linking the two equations and is provided when an archimedean copula with two parameters 
             is employed.}

\item{KeT}{Estimated Kendall's tau coefficient between the two equations.}

\item{n}{Sample size.}

\item{n.sel}{Number of selected observations in the sample selection case.}

\item{X1, X2}{Design matrices associated with the first and second linear predictors.}

\item{X1.d2, X2.d2}{Number of columns of the design matrices for equations 1 and 2.}

\item{l.sp1, l.sp2}{Number of smooth components in equations 1 and 2.}

\item{He}{Penalized hessian. This is the same as \code{HeSh} for unpenalized models.}

\item{HeSh}{Unpenalized hessian.}

\item{Vb}{Inverse of \code{He}. This corresponds to the Bayesian variance-covariance matrix 
          used for `confidence' interval calculations.}

\item{F}{This is given by \code{Vb*HeSh}.} 

\item{fp}{If \code{TRUE}, then a fully parametric model was fitted.} 

\item{aut.sp}{If \code{FALSE}, then a fully parametric model with regression splines was fitted.}  
 
\item{t.edf}{Total degrees of freedom of the estimated bivariate model. It is calculated as \code{sum(diag(F))}.}

\item{edf1, edf2}{Degrees of freedom for the two equations of the estimated bivariate model. It is calculated when splines are used.}

\item{bs.mgfit}{List of values and diagnostics extracted from \code{magic}.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm stopped before reaching the maximum number of iterations allowed.}
                
\item{wor.c}{Working model quantities given by the square root of 
             the weight matrix times the pseudo-data vector and design matrix, \code{rW.Z} and \code{rW.X}.}                

\item{p11, p10, p01, p00}{Distribution function of a bivariate normal with zero means, unit variances and correlation \code{rho} evaluated at (y1 = 1, y2 = 1), (y1 = 1, y2 = 0), (y1 = 0, y2 = 1) and (y1 = 0, y2 = 0).}
                       
\item{p0}{Distribution function of a standardised univariate normal evaluated at y1 = 0. This is only provided when \code{selection = TRUE}.}                       

\item{eta1, eta2}{Estimated linear predictors for the two equations.}

\item{y1, y2}{Responses of the two equations.}

\item{sel}{This is used for internal calculations.}

\item{BivD, nu}{These are used for internal calculations.}

\item{PL, eqPL, valPL, fitPL, spPL}{These are used for internal calculations.}

\item{xi1, xi2}{Shape parameters of the link functions of the two equations if an asymmetric link approach is used.}

\item{logLik}{Value of the (unpenalized) log-likelihood evaluated at the (penalized or unpenalized) parameter 
            estimates.}

\item{nC}{Copula identifier. This is used for internal calculations.}

\item{H.n}{This is used for internal calculations.}

\item{pPen1, pPen2}{This is used for internal calculations. List specifying any penalties to be applied to the parametric model terms of the model equations.}

\item{good}{Indicator variable indicating the observations actually used in model fitting.}

\item{y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1}{These are used for internal calculations.}

\item{qu.mag, gp1, gp2}{These are used for internal calculations.}

\item{X2s}{Full design matrix of second equation in sample selection case.}


}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{SemiParBIVProbit}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}
}











