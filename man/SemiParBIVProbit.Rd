\name{SemiParBIVProbit}
\alias{SemiParBIVProbit}

\title{Semiparametric Bivariate Regression Models}

\description{ 
\code{SemiParBIVProbit} can be used to fit bivariate binary models where the linear predictors of the two main equations can 
be flexibly specified using 
parametric and regression spline components. Several bivariate copula distributions are supported. During the model fitting process, the 
possible presence of associated error equations, endogeneity, non-random sample selection or partial observability is accounted for. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines. Note that, if it makes sense, the dependence parameter of the copula
distribution can be specified as a function of covariates.

\code{SemiParBIVProbit} can also be used to fit bivariate models in which one response is binary and the other is continuous. This case
is relevant for estimating, for instance, the effect that a binary endogenous variable has on a continuous response or viceversa. 
There are many continuous distributions that can be employed within this context and we plan to include several 
options on a case-study basis. Other bivariate models will be implemented from time to time, depending on the case study at hand.

}



\usage{
SemiParBIVProbit(formula, data = list(), weights = NULL, subset = NULL,  
                 start.v = NULL, Model = "B", BivD = "N", 
                 margins = c("probit","probit"), gamlssfit = FALSE,
                 fp = FALSE, hess = TRUE, infl.fac = 1, 
                 rinit = 1, rmax = 100, 
                 iterlimsp = 50, tolsp = 1e-6,
                 gc.l = FALSE, parscale, extra.regI = "t")
}



\arguments{ 

\item{formula}{In the basic setup this will be a list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms 
               are used to specify smooth smooth functions of 
                  predictors. \code{SemiParBIVProbit} supports the use shrinkage smoothers for variable
                  selection purposes. See the examples below and the documentation of \code{mgcv} for further 
                  details on formula specifications. Note that 
                  if \code{Model = "BSS"} then the first formula MUST refer to the selection equation. Furthermore,
                  if it makes sense, a third equation for the dependence parameter can be specified (see Example 1 below).
                  When one outcome is binary and the other continuous then the first equation must refer to the binary outcome whereas 
                  the second the continuous one. For the binary-continuous bivariate model, two more equations can be specified, one modelling 
                  the dispersion parameter and the other the dependence parameter.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParBIVProbit} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{start.v}{Starting values for the parameters of the two equations and association coefficient 
               can be provided here.} 

\item{Model}{It indicates the type of model to be used in the analysis. Possible values are "B" (bivariate model), 
             "BSS" (bivariate model with non-random sample selection), "BPO" (bivariate model with partial observability)
             and "BPO0" (bivariate model with partial observability and zero correlation).}

\item{margins}{It indicates the distributions used for the two margins. The first is always probit which means that the response of the 
               first equation is always assumed to be binary. The response for the second equation can be binary (case in which a 
               "probit" link is used), normal ("N"), log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), 
               logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), Dagum ("DAGUM").}
               
\item{gamlssfit}{In the case of continuous margin for the outcome equation, if \code{gamlssfit = TRUE} then 
                 a gamlss if fitted for the outcome equation. This is useful to calculate, for instance, the average
                 treatment effect from the model which does not account for endogeneity/associated equation errors.
               }               

\item{BivD}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), survival Clayton,
             rotated Clayton (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), and Frank. Note that 
             for empirical applications, options N, one of Clayton, Joe or Gumbel (with all respective rotated versions), and F may suffice.
             This is because the Clayton, Joe and Gumbel copulas are very similar; for instance, C90, J270 and G270 are very close.
             Nevertheless, there might be situations in which the use of a specific copula may result in more stable computation.}            

\item{fp}{If \code{TRUE} then a fully parametric model with unpenalised regression splines if fitted. See the example below.}              

\item{hess}{If \code{FALSE} then the expected/Fisher (rather than observed) information matrix is employed. The Fisher information matrix
            is not available for cases different from binary treatment and binary outcome.}

\item{infl.fac}{Inflation factor for the model degrees of freedom in the UBRE score. Smoother models can be obtained setting 
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
                  only makes sense when the information matrix is semi-positive or positive definite) and "sED" (symmetric eigen-decomposition).}  
     
                      
}


\details{

The bivariate models considered in this package consist of two model equations which depend on flexible linear predictors and
whose association between the responses is modelled through parameter \eqn{\theta}{\theta} of a standardised bivariate normal 
distribution or that of a bivariate copula distribution. The linear predictors of the two equations are flexibly specified using 
parametric components and smooth functions of covariates. The same can be done for the dependence parameter if it makes sense.
Replacing the smooth components with their regression spline expressions 
yields a fully parametric bivariate model. In principle, classic 
maximum likelihood estimation can be employed. However, to avoid overfitting, penalized likelihood maximization is 
employed instead. Here, the use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support 
from the data. The trade-off between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimize the approximate Un-Biased Risk Estimator (UBRE), which can also be viewed as an approximate AIC.

Automatic smoothing parameter estimation is integrated 
using a performance-oriented iteration approach (Gu, 1992; Wood, 2004). Roughly speaking, at each iteration, (i) the penalized weighted least squares 
problem is solved and (ii) the smoothing parameters of that problem estimated by approximate UBRE. Steps (i) and (ii) are iterated until 
convergence. Details of the underlying fitting methods are given in Radice, Marra and Wojtys (in press). Releases previous to 3.2-7
were based on the algorithms detailed in Marra and Radice (2011, 2013). 

}


\value{
The function returns an object of class \code{SemiParBIVProbit} as described in \code{SemiParBIVProbitObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Gu C. (1992), Cross validating non-Gaussian data. \emph{Journal of Computational and Graphical Statistics}, 1(2), 169-179.

Marra G. (2013), On P-values for Semiparametric Bivariate Probit Models. \emph{Statistical Methodology}, 10(1), 23-28.

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian 
Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G., Radice R. and Missiroli S. (2014), Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Computational Statistics}, 29(3-4), 715-741.

Marra G., Radice R. and Filippou P. (in press), Regression Spline Bivariate Probit Models: A Practical Approach to Testing for Exogeneity. \emph{Communications in Statistics - Simulation and Computation}.

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.

Radice R., Marra G. and M. Wojtys (in press), Copula Regression Spline Models for Binary Outcomes. \emph{Statistics and Computing}.

Marra G. and Wood S.N. (2011), Practical Variable Selection for Generalized Additive Models. \emph{Computational Statistics and Data Analysis}, 55(7), 2372-2387.

Marra G. and Wood S.N. (2012), Coverage Properties of Confidence Intervals for Generalized Additive Model Components. \emph{Scandinavian Journal of Statistics}, 39(1), 53-74.

Poirier D.J. (1980), Partial Observability in Bivariate Probit Models. \emph{Journal of Econometrics}, 12, 209-217.

Wood S.N. (2004), Stable and efficient multiple smoothing parameter estimation for generalized additive models. \emph{Journal of the American Statistical Association}, 99(467), 673-686.

}


\section{WARNINGS}{

In our experience, convergence failure may occur, for instance, when the number of observations is low, there are important mistakes 
in the model specification (i.e., using C90 when the model errors are positively associated), or the model is 'over-complicated'. Convergence 
failure may also mean that 
an infinite cycling between steps (i) and (ii) occurs. In this case, the smoothing parameters are set to the values 
obtained from the non-converged algorithm (\code{\link{conv.check}} will notify the user if this is the case). In such cases, we 
recommend re-specifying the model, using some extra regularisation (see \code{extra.regI}) and/or using some
rescaling (see \code{parscale}).

In the contexts of endogeneity and non-random sample selection, extra attention is required when specifying
the dependence parameter as a function of covariates. This is because in these situations the dependence coefficient mainly models the 
association between the unobserved confounders in the two equations. This option would make sense, for instance, when it is believed that the  
association coefficient varies across geographical areas.


} 

\seealso{
\code{\link{AT}}, \code{\link{OR}}, \code{\link{RR}}, \code{\link{adjCov}}, \code{\link{prev}}, \code{\link{gt.bpm}}, \code{\link{LM.bpm}}, \code{\link{VuongClarke}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbitObject}}, \code{\link{conv.check}}, \code{\link{summary.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}
}

\examples{

library(SemiParBIVProbit)

############
## EXAMPLE 1
############
## Generate data
## Correlation between the two equations 0.5 - Sample size 400 

set.seed(0)

n <- 400

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*x1 + f2(x2) + u[,2] > 0, 1, 0)

dataSim <- data.frame(y1, y2, x1, x2, x3)

#
#

## CLASSIC BIVARIATE PROBIT

out  <- SemiParBIVProbit(list(y1 ~ x1 + x2 + x3, 
                              y2 ~ x1 + x2 + x3), 
                         data = dataSim)
conv.check(out)
summary(out)
AIC(out)
BIC(out)

## SEMIPARAMETRIC BIVARIATE PROBIT

## "cr" cubic regression spline basis      - "cs" shrinkage version of "cr"
## "tp" thin plate regression spline basis - "ts" shrinkage version of "tp"
## for smooths of one variable, "cr/cs" and "tp/ts" achieve similar results 
## k is the basis dimension - default is 10
## m is the order of the penalty for the specific term - default is 2
## For COPULA models use BivD argument 

out  <- SemiParBIVProbit(list(y1 ~ x1 + s(x2, bs = "tp", k = 10, m = 2) + s(x3), 
                              y2 ~ x1 + s(x2) + s(x3)),  
                         data = dataSim)
conv.check(out)
summary(out, cm.plot = TRUE)
AIC(out)


## estimated smooth function plots - red lines are true curves

x2 <- sort(x2)
f1.x2 <- f1(x2)[order(x2)] - mean(f1(x2))
f2.x2 <- f2(x2)[order(x2)] - mean(f2(x2))
f3.x3 <- rep(0, length(x3))

par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(out, eq = 1, select = 1, seWithMean = TRUE, scale = 0)
lines(x2, f1.x2, col = "red")
plot(out, eq = 1, select = 2, seWithMean = TRUE, scale = 0)
lines(x3, f3.x3, col = "red")
plot(out, eq = 2, select = 1, seWithMean = TRUE, scale = 0)
lines(x2, f2.x2, col = "red")
plot(out, eq = 2, select = 2, seWithMean = TRUE, scale = 0)
lines(x3, f3.x3, col = "red")

## p-values suggest to drop x3 from both equations, with a stronger 
## evidence for eq. 2. This can be also achieved using shrinkage smoothers

outSS <- SemiParBIVProbit(list(y1 ~ x1 + s(x2, bs = "ts") + s(x3, bs = "cs"), 
                               y2 ~ x1 + s(x2, bs = "cs") + s(x3, bs = "ts")), 
                          data = dataSim)
conv.check(outSS)                          

plot(outSS, eq = 1, select = 1, scale = 0, shade = TRUE)
plot(outSS, eq = 1, select = 2, ylim = c(-0.1,0.1))
plot(outSS, eq = 2, select = 1, scale = 0, shade = TRUE)
plot(outSS, eq = 2, select = 2, ylim = c(-0.1,0.1))

\dontrun{ 

## SEMIPARAMETRIC BIVARIATE PROBIT with association parameter 
## depending on covariates as well

outD <- SemiParBIVProbit(list(y1 ~ x1 + s(x2), 
                              y2 ~ x1 + s(x2),
                                 ~ x1 + s(x2) ), 
                          data = dataSim)
conv.check(outD)  
summary(outD)
outD$theta

plot(outD, eq = 1, seWithMean = TRUE)
plot(outD, eq = 2, seWithMean = TRUE)
plot(outD, eq = 3, seWithMean = TRUE)
graphics.off()

#
#

############
## EXAMPLE 2
############
## Generate data with one endogenous variable 
## and exclusion restriction

set.seed(0)

n <- 400

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

cov   <- rMVN(n, rep(0,2), Sigma)
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*y1 + f2(x2) + u[,2] > 0, 1, 0)

dataSim <- data.frame(y1, y2, x1, x2)

#

## Testing the hypothesis of absence of endogeneity... 

LM.bpm(list(y1 ~ x1 + s(x2), y2 ~ y1 + s(x2)), dataSim, Model = "B")

# p-value suggests presence of endogeneity, hence fit a bivariate model


## CLASSIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIVProbit(list(y1 ~ x1 + x2, 
                             y2 ~ y1 + x2), 
                        data = dataSim)
conv.check(out)                        
summary(out)
AIC(out); BIC(out)

## SEMIPARAMETRIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIVProbit(list(y1 ~ x1 + s(x2), 
                             y2 ~ y1 + s(x2)), 
                        data = dataSim)
conv.check(out)                        
summary(out)
AIC(out); BIC(out)

#

## Testing the hypothesis of absence of endogeneity post estimation... 

gt.bpm(out)

#
## average treatment effect, risk ration and odds ratio with CIs

mb(y1, y2, Model = "B")
AT(out, nm.end = "y1", hd.plot = TRUE) 
RR(out, nm.end = "y1") 
OR(out, nm.end = "y1") 
AT(out, nm.end = "y1", type = "univariate") 


## try a Clayton copula model... 

outC <- SemiParBIVProbit(list(y1 ~ x1 + s(x2), 
                              y2 ~ y1 + s(x2)), 
                         data = dataSim, BivD = "C0")
conv.check(outC)                         
summary(outC)
AT(outC, nm.end = "y1") 

## try a Joe copula model... 

outJ <- SemiParBIVProbit(list(y1 ~ x1 + s(x2), 
                              y2 ~ y1 + s(x2)), 
                         data = dataSim, BivD = "J0")
conv.check(outJ)
summary(outJ, cm.plot = TRUE)
AT(outJ, "y1") 


VuongClarke(out, outJ)

#
## recursive bivariate probit modelling with unpenalized splines 
## can be achieved as follows

outFP <- SemiParBIVProbit(list(y1 ~ x1 + s(x2, bs = "cr", k = 5), 
                               y2 ~ y1 + s(x2, bs = "cr", k = 6)), 
                          fp = TRUE, data = dataSim)
conv.check(outFP)                            
summary(outFP)

#
#################
## See also ?meps
#################

############
## EXAMPLE 3
############
## Generate data with a non-random sample selection mechanism 
## and exclusion restriction

set.seed(0)

n <- 2000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

SigmaC <- matrix(0.5, 3, 3); diag(SigmaC) <- 1
cov    <- rMVN(n, rep(0,3), SigmaC)
cov    <- pnorm(cov)
bi <- round(cov[,1]); x1 <- cov[,2]; x2 <- cov[,3]
  
f11 <- function(x) -0.7*(4*x + 2.5*x^2 + 0.7*sin(5*x) + cos(7.5*x))
f12 <- function(x) -0.4*( -0.3 - 1.6*x + sin(5*x))  
f21 <- function(x) 0.6*(exp(x) + sin(2.9*x)) 

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + u[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +         + u[, 2] > 0
yo <- y*(ys > 0)
  
dataSim <- data.frame(y, ys, yo, bi, x1, x2)

## Testing the hypothesis of absence of non-random sample selection... 

LM.bpm(list(ys ~ bi + s(x1) + s(x2), yo ~ bi + s(x1)), dataSim, Model = "BSS")

# p-value suggests presence of sample selection, hence fit a bivariate model

#
## SEMIPARAMETRIC SAMPLE SELECTION BIVARIATE PROBIT
## the first equation MUST be the selection equation

out <- SemiParBIVProbit(list(ys ~ bi + s(x1) + s(x2), 
                             yo ~ bi + s(x1)), 
                        data = dataSim, Model = "BSS")
conv.check(out)                          
gt.bpm(out)                        

## compare the two summary outputs
## the second output produces a summary of the results obtained when
## selection bias is not accounted for

summary(out)
summary(out$gam2)

## corrected predicted probability that 'yo' is equal to 1

mb(ys, yo, Model = "BSS")
prev(out, hd.plot = TRUE)
prev(out, type = "univariate", hd.plot = TRUE)

## estimated smooth function plots
## the red line is the true curve
## the blue line is the univariate model curve not accounting for selection bias

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)]-mean(f21(x1.s))

plot(out, eq = 2, ylim = c(-1.65,0.95)); lines(x1.s, f21.x1, col="red")
par(new = TRUE)
plot(out$gam2, se = FALSE, col = "blue", ylim = c(-1.65,0.95), 
     ylab = "", rug = FALSE)

#
#
## try a Clayton copula model... 

outC <- SemiParBIVProbit(list(ys ~ bi + s(x1) + s(x2), 
                              yo ~ bi + s(x1)), 
                         data = dataSim, Model = "BSS", BivD = "C0")
conv.check(outC)
summary(outC, cm.plot = TRUE)
prev(outC)

#
################
## See also ?hiv
################

############
## EXAMPLE 4
############
## Generate data with partial observability

set.seed(0)

n <- 10000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

y1 <- ifelse(-1.55 + 2*x1 + x2 + u[,1] > 0, 1, 0)
y2 <- ifelse( 0.45 - x3        + u[,2] > 0, 1, 0)
y  <- y1*y2

dataSim <- data.frame(y, x1, x2, x3)


## BIVARIATE PROBIT with Partial Observability

out  <- SemiParBIVProbit(list(y ~ x1 + x2, 
                              y ~ x3), 
                         data = dataSim, Model = "BPO")
conv.check(out)
summary(out)

# first ten estimated probabilities for the four events from object out

cbind(out$p11, out$p10, out$p00, out$p01)[1:10,]


# case with smooth function 
# (a bit more computationally intensive)  

f1 <- function(x) cos(pi*2*x) + sin(pi*x)

y1 <- ifelse(-1.55 + 2*x1 + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse( 0.45 - x3            + u[,2] > 0, 1, 0)
y  <- y1*y2

dataSim <- data.frame(y, x1, x2, x3)

out  <- SemiParBIVProbit(list(y ~ x1 + s(x2), 
                              y ~ x3), 
                         data = dataSim, Model = "BPO")

conv.check(out)
summary(out, cm.plot = TRUE)


# plot estimated and true functions

x2 <- sort(x2); f1.x2 <- f1(x2)[order(x2)] - mean(f1(x2))
plot(out, eq = 1, scale = 0); lines(x2, f1.x2, col = "red")

#
################
## See also ?war
################

############
## EXAMPLE 5
############
## Generate data with one endogenous binary variable 
## and continuous outcome

set.seed(0)

n <- 1000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

cov   <- rMVN(n, rep(0,2), Sigma)
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <-        -0.25 - 1.25*y1 + f2(x2) + u[,2] 

dataSim <- data.frame(y1, y2, x1, x2)


## RECURSIVE Model

resp.check(y2, margin = "N") 

out <- SemiParBIVProbit(list(y1 ~ x1 + x2, 
                             y2 ~ y1 + x2), 
                        data = dataSim, margins = c("probit","N"))
conv.check(out)                        
summary(out)
post.check(out)

## SEMIPARAMETRIC RECURSIVE Model

out <- SemiParBIVProbit(list(y1 ~ x1 + s(x2), 
                             y2 ~ y1 + s(x2), ~ 1, ~ 1), data = dataSim, 
                             margins = c("probit","N"), gamlssfit = TRUE)
conv.check(out)                        
summary(out)
AT(out, nm.end = "y1")
AT(out, nm.end = "y1", type = "univariate")

#
#

############
## EXAMPLE 6
############
## Generate data with one endogenous continuous exposure 
## and binary outcome

set.seed(0)

n <- 1000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

cov   <- rMVN(n, rep(0,2), Sigma)
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2) 

y1 <-        -0.25 - 2*x1    + f2(x2) + u[,2] 
y2 <- ifelse(-0.25 - 0.25*y1 + f1(x2) + u[,1] > 0, 1, 0)

dataSim <- data.frame(y1, y2, x1, x2)


out <- SemiParBIVProbit(list(y2 ~ y1 + s(x2), 
                             y1 ~ x1 + s(x2), ~ 1, ~ 1), data = dataSim, 
                             margins = c("probit","N"))
conv.check(out)                        
summary(out)
AT(out, nm.end = "y1", hd.plot = TRUE)
AT(out, nm.end = "y1", type = "univariate")
RR(out, nm.end = "y1", prob.plot = TRUE)
RR(out, nm.end = "y1", type = "univariate")
OR(out, nm.end = "y1", hd.plot = TRUE, prob.plot = TRUE)
OR(out, nm.end = "y1", type = "univariate")

}
}

\keyword{asymmetric link function}
\keyword{non-random sample selection}
\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{partial observability} 
\keyword{endogeneity} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







