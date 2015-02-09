\name{SemiParBIVProbitObject}
\alias{SemiParBIVProbitObject}

\title{Fitted SemiParBIVProbit object}

\description{A fitted semiparametric bivariate probit object returned by function \code{SemiParBIVProbit} and of class "SemiParBIVProbit".}


\value{ 


\item{fit}{List of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$gradient}, \code{fit$Fisher} and \code{fit$S.h} return
           the gradient vector, Fisher information (when used) and overall penalty matrix scaled 
           by its smoothing parameters, for the fitted bivariate probit model. See the 
           documentation of \code{trust} for details on the diagnostics provided.}

\item{gam1}{Univariate fit for equation 1. See the documentation of \code{mgcv} for full details.}

\item{gam2}{Univariate fit for equation 2.}

\item{gam2.1}{Univariate fit for equation 2, estimated using an adaptation of the Heckman sample selection correction procedure 
              when \code{Method = "BSS"}.}
              
\item{gam3}{Univariate fit for equation 3. This is available when the association parameter is modelled as functions of covariates.}              

\item{coefficients}{The coefficients of the fitted model. They are given in the following 
                    order: parametric and regression spline (if present)
                    coefficients for the first equation, parametric and regression spline
                    coefficients for the second equation, and association coefficient.}

\item{weights}{Prior weights used during model fitting.}

\item{sp}{Estimated smoothing parameters of the smooth components.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{iter.if}{Number of iterations performed in the initial step of the algorithm.}

\item{iter.inner}{Number of iterations performed within the smoothing parameter estimation step.}

\item{rho}{Estimated correlation coefficient between the two equations, when a bivariate normal model is employed.}

\item{theta}{Estimated copula parameter linking the two equations.}

\item{n}{Sample size.}

\item{n.sel}{Number of selected observations in the sample selection model case.}

\item{X1, X2}{Design matrices associated with the first and second linear predictors.}

\item{X1.d2, X2.d2}{Number of columns of \code{X1} and \code{X2}.}

\item{l.sp1, l.sp2}{Number of smooth components in equations 1 and 2.}

\item{He}{Penalized -hessian/Fisher. This is the same as \code{HeSh} for unpenalized models.}

\item{HeSh}{Unpenalized -hessian/Fisher.}

\item{Vb}{Inverse of \code{He}. This corresponds to the Bayesian variance-covariance matrix 
          used for confidence/credible interval calculations.}

\item{F}{This is given by \code{Vb*HeSh}.} 
 
\item{t.edf}{Total degrees of freedom of the estimated bivariate model. It is calculated as \code{sum(diag(F))}.}

\item{edf1, edf2, edf3}{Degrees of freedom for the two equations of the fitted bivariate model (and for the third
                        equation as well if the dependence parameter is modelled as function of covariates). They  
                        are calculated when splines are used.}

\item{bs.mgfit}{List of values and diagnostics extracted from \code{magic} in \code{mgcv}.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm stopped before reaching the maximum number of iterations allowed.}
                
\item{wor.c}{Working model quantities.}                

\item{p11, p10, p01, p00}{Model probabilities evaluated at (y_1 = 1, y_2 = 1), (y_1 = 1, y_2 = 0), (y_1 = 0, y_2 = 1) and (y_1 = 0, y_2 = 0).}
                       
\item{p0}{Distribution function of a standardised univariate normal evaluated at y_1 = 0. This is only provided when \code{Method = "BSS"}.}                       

\item{eta1, eta2, eta3}{Estimated linear predictors for the two equations (as well as the third equation if present).}

\item{y1, y2}{Responses of the two equations.}

\item{xi1, xi2}{Estimated shape parameters of the link functions of the two equations if an asymmetric link approach is used.}

\item{logLik}{Value of the (unpenalized) log-likelihood evaluated at the (penalized or unpenalized) parameter 
            estimates.}

\item{good}{Indicator variable indicating the observations actually used in model fitting.}

\item{respvec}{List containing response vectors.}

\item{X2s}{Full design matrix of outcome equation in sample selection case.}

\item{OR, GM}{Odds ratio and Gamma measure. See \code{\link{summary.SemiParBIVProbit}} for details.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{SemiParBIVProbit}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}
}











