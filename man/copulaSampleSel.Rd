\name{copulaSampleSel}
\alias{copulaSampleSel}

\title{Semiparametric Copula Bivariate Regression Models with Non-Random Sample Selection}

\description{ 
\code{copulaSampleSel} can be used to fit bivariate sample selection models where the linear predictors of the two main equations can 
be flexibly specified using a variety of covariate effects. Several bivariate copula distributions are supported. During the 
model fitting process, the 
possible presence of non-random sample selection is accounted for. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines. Note that, if it makes sense, the dependence parameter of the employed bivariate 
distribution can be specified as a function of covariates.

}



\usage{
copulaSampleSel(formula, data = list(), weights = NULL, subset = NULL,  
                 BivD = "N", margins = c("probit", "N"), gamlssfit = FALSE,
                 fp = FALSE, infl.fac = 1, 
                 rinit = 1, rmax = 100, 
                 iterlimsp = 50, tolsp = 1e-07,
                 gc.l = FALSE, parscale, extra.regI = "t")
}



\arguments{ 

\item{formula}{In the basic setup this will be a list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms 
               are used to specify smooth smooth functions of 
                  predictors. For the case of more than two equations see the example below and the documentation 
                  of \code{SemiParBIVProbit()} for more details. Note that 
                  the first formula MUST refer to the selection equation.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{copulaSampleSel} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}


\item{margins}{It indicates the distributions used for the two margins. The first is one of "probit", "logit", "cloglog" which refer 
               to the link function of the first equation whose response is always assumed to be binary. The response 
               for the second equation can be normal ("N"), 
               log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), 
               logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), 
               gamma with identity link for the location parameter ("GAi"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic distribution).}
               
\item{gamlssfit}{If \code{gamlssfit = TRUE} then 
                 a gamlss if fitted for the outcome equation. This is may used for obtaining better calibrated starting values, for instance.
               }               

\item{BivD}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F", "AMH", "FGM" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), survival Clayton,
             rotated Clayton (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), Frank, Ali-Mikhail-Haq and 
             Farlie-Gumbel-Morgenstern.}            

\item{fp}{If \code{TRUE} then a fully parametric model with unpenalised regression splines if fitted. See the example below.}              

\item{infl.fac}{Inflation factor for the model degrees of freedom in the approximate AIC. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds. See the documentation
             of \code{trust} for further details.}

\item{rmax}{Maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest 
            descent path.} 

\item{iterlimsp}{A positive integer specifying the maximum number of loops to be performed before the smoothing parameter estimation 
                 step is terminated.} 

\item{tolsp}{Tolerance to use in judging convergence of the algorithm when automatic smoothing parameter estimation is used.} 

\item{gc.l}{This is relevant when working with big datasets. If \code{TRUE} then the garbage collector is called more often than it is 
usually done. This keeps the memory footprint down but it will slow down the routine.}
      
\item{parscale}{The algorithm will operate as if optimizing objfun(x / parscale, ...) where parscale is a scalar. If missing then no 
rescaling is done. See the 
documentation of \code{trust} for more details.}    

\item{extra.regI}{If "t" then regularization as from \code{trust} is applied to the information matrix if needed. 
                  If different from "t" then extra regularization is applied via the options "pC" (pivoted Choleski - this
                  will only work when the information matrix is semi-positive or positive definite) and "sED" (symmetric eigen-decomposition).}  
     
                      
}


\details{

The underlying algorithm is based on an extension of the procedure used for \code{SemiParBIVProbit()}. For more details
see ?SemiParBIVProbit. 

This function works as \code{SemiParSampleSel()} in \code{SemiParSampleSel} and has been included in \code{SemiParBIVProbit} (which
already included sample selection models for binary outcomes) for the user's convenience (given some
requests). \code{copulaSampleSel()} allows for the use of 
many continuous distributions and different link functions for the selection equation whereas \code{SemiParSampleSel()} allows only for
a probit link and normal or gamma outcome. The latter, however, includes many discrete distributions. 

If there are factors in the model, before fitting, the user has to ensure 
that the numbers of factor variables' levels in the selected sample 
are the same as those in the complete dataset. Even if a model could be fitted in such a situation,
the model may produce fits which are 
not coherent with the nature of the correction sought. For more details see \code{?SemiParBIVProbit}. 

}


\value{
The function returns an object of class \code{copulaSampleSel} as described in \code{copulaSampleSelObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Marra G. and Radice R. (2013), Estimation of a Regression Spline Sample Selection Model. \emph{Computational Statistics and Data Analysis}, 61, 158-173.

Wojtys M. and Marra G. (submitted). Copula-Based Generalized Additive Models with Non-Random Sample Selection.

}


\section{WARNINGS}{

Convergence failure may sometimes occur. Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be 0 and the latter positive definite.
\code{SemiParBIVProbit()} will produce some warnings if there is a convergence issue.

In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). Using \code{gamlssfit = TRUE} is typically more effective than the first two options as
this will provide better calibrated starting values as compared to those obtained from the default starting value procedure.
The default option is, however, \code{gamlssfit = FALSE} only because it tends to be computationally cheaper and because the 
default starting value procedure has typically been found to do a satisfactory job in most cases. 
(The results obtained when using 
\code{gamlssfit = FALSE} and \code{gamlssfit = TRUE} could also be compared to check if starting values make any difference.)

The above suggestions may help, especially the latter option. However, the user should also consider
re-specifying the model and/or using a diferrent dependence structure and/or checking that the chosen marginal 
distribution fit the responses well.
In our experience, we found that convergence failure typically occurs 
when the model has been misspecified and/or the sample size and/or number of selected observations is low compared 
to the complexity of the model. Examples
of misspecification include using a Clayton copula rotated by 90 degrees when a positive
association between the margins is present instead, using marginal distributions that do not fit
the responses (again, this is a bit more relevant when one of the two responses is continuous), and 
employing a copula which does not accommodate the type and/or strength of
the dependence between the margins (e.g., using AMH when the association between the margins is strong).
It is also worth bearing in mind that the use of a three parameter marginal distribution requires the data
to be more informative than a situation in which a two parameter distribution is used instead.

When comparing competing models (for instance, by keeping the linear predictor specifications
fixed and changing the copula), if the computing time for a set of alternatives
is considerably higher than that of another set then it may mean that 
the more computationally demanding models are not able to fit the data very well (as a higher number of 
iterations is required to reach convergence). As a practical check, this may be verified by
fitting all competing models and, provided convergence is achieved, comparing their respective AIC and BICs, for instance.

Extra attention is required when specifying
the dependence parameter as a function of covariates. This is because in these situations the dependence parameter mainly models the 
association between the unobserved confounders in the two equations. Therefore, this option would make sense when it 
is believed that the  
strength of the association between the unobservables in the two equations varies based on some grouping factor or across geographical 
areas, for instance.

} 

\seealso{
\code{\link{copulaReg}}, \code{\link{SemiParBIVProbit}}, \code{\link{adjCov}}, \code{\link{VuongClarke}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{copulaSampleSelObject}}, \code{\link{conv.check}}, \code{\link{summary.copulaSampleSel}}, \code{\link{predict.SemiParBIVProbit}}
}

\examples{

\dontrun{ 

library(SemiParBIVProbit)

######################################################################
## Generate data
## Correlation between the two equations and covariate correlation 0.5 
## Sample size 2000 
######################################################################

set.seed(0)

n <- 2000

rh <- 0.5      

sigmau <- matrix(c(1, rh, rh, 1), 2, 2)
u      <- rMVN(n, rep(0,2), sigmau)

sigmac <- matrix(rh, 3, 3); diag(sigmac) <- 1
cov    <- rMVN(n, rep(0,3), sigmac)
cov    <- pnorm(cov)

bi <- round(cov[,1]); x1 <- cov[,2]; x2 <- cov[,3]
  
f11 <- function(x) -0.7*(4*x + 2.5*x^2 + 0.7*sin(5*x) + cos(7.5*x))
f12 <- function(x) -0.4*( -0.3 - 1.6*x + sin(5*x))  
f21 <- function(x) 0.6*(exp(x) + sin(2.9*x)) 

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + u[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +           u[, 2]
yo <- y*(ys > 0)
  
dataSim <- data.frame(ys, yo, bi, x1, x2)

## CLASSIC SAMPLE SELECTION MODEL
## the first equation MUST be the selection equation

resp.check(yo[ys > 0], "N")

out <- copulaSampleSel(list(ys ~ bi + x1 + x2, 
                            yo ~ bi + x1), 
                       data = dataSim)
conv.check(out)
post.check(out)
summary(out)

AIC(out)
BIC(out)


## SEMIPARAMETRIC SAMPLE SELECTION MODEL

## "cr" cubic regression spline basis      - "cs" shrinkage version of "cr"
## "tp" thin plate regression spline basis - "ts" shrinkage version of "tp"
## for smooths of one variable, "cr/cs" and "tp/ts" achieve similar results 
## k is the basis dimension - default is 10
## m is the order of the penalty for the specific term - default is 2

out <- copulaSampleSel(list(ys ~ bi + s(x1, bs = "tp", k = 10, m = 2) + s(x2), 
                            yo ~ bi + s(x1)), 
                       data = dataSim)
conv.check(out) 
post.check(out)
AIC(out)

## compare the two summary outputs
## the second output produces a summary of the results obtained when only 
## the outcome equation is fitted, i.e. selection bias is not accounted for

summary(out)
summary(out$gam2)

## estimated smooth function plots
## the red line is the true curve
## the blue line is the naive curve not accounting for selection bias

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)] - mean(f21(x1.s))

plot(out, eq = 2, ylim = c(-1, 0.8)); lines(x1.s, f21.x1, col = "red")
par(new = TRUE)
plot(out$gam2, se = FALSE, lty = 3, lwd = 2, ylim = c(-1, 0.8), 
     ylab = "", rug = FALSE)


## SEMIPARAMETRIC SAMPLE SELECTION MODEL with association 
## and dispersion parameters 
## depending on covariates as well

eq.mu.1   <- ys ~ bi + s(x1) + s(x2)
eq.mu.2   <- yo ~ bi + s(x1)
eq.sigma2 <-    ~ bi
eq.theta  <-    ~ bi + x1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma2, eq.theta)

out <- copulaSampleSel(fl, data = dataSim)
conv.check(out)   
post.check(out)
summary(out)
out$sigma2
out$theta

jc.probs(out, 0, 0.3, intervals = TRUE)[1:4,]

outC0 <- copulaSampleSel(fl, data = dataSim, BivD = "C0")
conv.check(outC0)
post.check(outC0)
AIC(out, outC0)
BIC(out, outC0)

#
#

#######################################################
## example using Gumbel copula and normal-gamma margins
#######################################################

y  <- rgamma(n, shape = 1/2^2, exp(-0.68 - 1.5*bi + f21(x1) + u[, 2])*2^2)
yo <- y*(ys > 0)
  
dataSim <- data.frame(ys, yo, bi, x1, x2)


out <- copulaSampleSel(list(ys ~ bi + s(x1) + s(x2), 
                            yo ~ bi + s(x1)), 
                        data = dataSim, BivD = "G0", 
                        margins = c("probit", "GA"))
conv.check(out)
post.check(out)
summary(out)


#
#
}
}


\keyword{non-random sample selection}
\keyword{semiparametric bivariate modelling} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







