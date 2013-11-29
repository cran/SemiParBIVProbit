\name{SemiParBIVProbit}
\alias{SemiParBIVProbit}

\title{Semiparametric Bivariate Probit Modelling}

\description{ 
\code{SemiParBIVProbit} can be used to fit bivariate probit models where the linear predictors can be flexibly specified using parametric and 
regression spline components, and parametric or nonparametric random effects. The use of several 
bivariate copula distributions and of asymmetric link functions is supported. During the model fitting process, the possible presence of 
correlated error equations, endogeneity or non-random sample selection is accounted for. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines (isotropic). The current implementation does not support scale invariant tensor 
product smooths. 

}



\usage{
SemiParBIVProbit(formula.eq1, formula.eq2, data=list(), weights=NULL,  
                 start.v=NULL, BivD="N", nu=3, PL="P", eqPL="both", xi1=1, xi2=1,
                 selection=FALSE, H.n=TRUE, 
                 gamma=1, aut.sp=TRUE, fp=FALSE, RE=FALSE, RE.type="N", 
                 NGQ=10, K=2, id=NULL, e.npRE=TRUE, pPen1 = NULL, pPen2 = NULL, 
                 rinit=1, rmax=100, fterm=sqrt(.Machine$double.eps), 
                 mterm=sqrt(.Machine$double.eps), iterlimsp=50, pr.tolsp=1e-6,    
                 control.sp=list(maxit=50,tol=1e-6,step.half=25,
                                 rank.tol=sqrt(.Machine$double.eps)) )
}


\arguments{ 

\item{formula.eq1}{A formula for equation 1. \code{s} terms are used to specify smooth smooth functions of 
                  predictors. \code{SemiParBIVProbit} supports the use shrinkage smoothers for variable
                  selection purposes. See the examples below and the documentation of \code{mgcv} for further 
                  details on formula specifications. Note that 
                  if \code{selection=TRUE} then the formula MUST refer to the selection equation.} 

\item{formula.eq2}{A formula for equation 2.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParBIVProbit} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{start.v}{Starting values for the parameters of the two equations and association coefficient 
               can be provided here.} 

\item{BivD}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F", "T" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), survival Clayton,
             rotated Clayton (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), Frank, Student-t.
             Archimedean copulas with two parameters are also implemented; these are the Clayton-Gumbel and its rotated 
             versions ("BB1.0","BB1.90","BB1.180","BB1.270"), Joe-Gumbel ("BB6.0","BB6.90","BB6.180","BB6.270"), 
             Joe-Clayton ("BB7.0","BB7.90","BB7.180","BB7.270") and Joe-Frank ("BB8.0","BB8.90","BB8.180","BB8.270").} 

\item{nu}{Degrees of freedom for Student-t.}

\item{PL}{Power link function to be employed. Possible choices are "P" (classic probit), "PP" (power probit) and "RPP" (reciprocal power probit).}

\item{eqPL}{Equation(s) on which the power link function approach should be applied. Possible choices are "both", "first" and "second".}

\item{xi1,xi2}{Power or shape parameters of the link functions of the two equations if a power link approach is used.}

\item{selection}{If \code{TRUE}, then the numerical routine for bivariate probit modelling in the presence of non-random sample selection is employed.}

\item{H.n}{For \code{BivD="N"} only, if \code{FALSE} then the expected (rather than observed) information matrix 
is employed.}

\item{gamma}{Inflation factor for the model degrees of freedom in the UBRE score. Smoother models can be obtained setting 
             this parameter to a value greater than 1. Typically \code{gamma=1.4} achieves this.}

\item{aut.sp}{If \code{TRUE}, then automatic multiple smoothing parameter selection is carried out. If \code{FALSE}, then smoothing parameters 
              are set to the values obtained from the univariate fits.}

\item{fp}{If \code{TRUE}, then a fully parametric model with regression splines if fitted. See the example below.}  
 
\item{RE}{If \code{TRUE}, then the routine for bivariate probit modelling with random effects is employed. In this case, 
            data have to be provided ordered by cluster.}

\item{RE.type}{If \code{RE=TRUE}, then two choices are possible: nonparametric random effects ("NP") and bivariate normal random effects ("N").}

\item{NGQ}{Number of quadrature points.}

\item{K}{If \code{RE.type="NP"}, then \code{K} represents the number of bivariate mass points and must be supplied. Default is 3.} 

\item{id}{If \code{RE=TRUE}, then the individual random effect identifier must be supplied.} 

\item{e.npRE}{If \code{RE.type="NP"}, then \code{e.npRE} indicates whether some extra iterations are carried out.}

\item{pPen1,pPen2}{Optional list specifying any penalties to be applied to the parametric model terms of equations 1 and 2.}

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds. See the documentation
             of \code{trust} for further details.}

\item{rmax}{Maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest 
            descent path.} 

\item{fterm}{Positive scalar giving the tolerance at which the difference in objective function values in a step is considered close 
             enough to zero to terminate the algorithm.}

\item{mterm}{Positive scalar giving the tolerance at which the two-term Taylor-series approximation to the difference in objective 
             function values in a step is considered close enough to zero to terminate the algorithm.}

\item{iterlimsp}{A positive integer specifying the maximum number of loops to be performed before the smoothing parameter estimation step is terminated.} 

\item{pr.tolsp}{Tolerance to use in judging convergence of the algorithm when automatic smoothing parameter selection is used.} 

\item{control.sp}{List containing iteration control constants for the automatic smoothing parameter selection method. These 
                are: \code{maxit}: maximum number of iterations of the 
               \code{magic} algorithm; \code{tol}: tolerance to use in judging convergence; \code{step.half}: if a trial step fails then 
               the method tries halving it up to a maximum of \code{step.half} times; \code{rank.tol}: constant used to test for numerical rank 
               deficiency of the problem. See the documentation of \code{magic} in \code{mgcv} for further details.}



}


\details{

The bivariate model has a probit link for each of the two equations, and models the association
between the responses through the correlation parameter \eqn{\rho}{\rho} of a standardised bivariate normal distribution, or that 
of a bivariate copula distribution, \eqn{\theta}{\theta}. In a 
semiparametric bivariate probit 
model the linear predictors are flexibly specified using 
parametric components and smooth functions of covariates. Replacing the smooth components with their regression spline expressions 
yields a fully parametric bivariate probit model. In principle, classic 
maximum likelihood estimation can be employed. However, to avoid overfitting, penalized likelihood maximization has 
to be employed instead. Here the use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support 
from the data. The trade-off between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimize the approximate Un-Biased Risk Estimator (UBRE).

Automatic smoothing parameter selection is integrated 
using a performance-oriented iteration approach (Gu, 1992; Wood, 2004). At each iteration, (i) the penalized weighted least squares 
problem is solved and (ii) the smoothing parameters of that problem estimated by approximate UBRE. Steps (i) and (ii) are iterated until 
convergence. Details of the underlying fitting methods are given in Radice, Marra and Wojtys (submitted). Releases previous to 3.2-7
were based on the algorithms detailed in Marra and Radice (2011, 2013). 

The use of parametric and nonparametric random effects is allowed for via the option \code{RE}. The nonparametric case
is illustrated in Example 4 (Marra et al., in press), whereas the parametric one, in Example 5, will be completed soon. Asymmetric link functions
are also available; see Marra and Radice (submitted) for details.

}


\value{
The function returns an object of class \code{SemiParBIVProbit} as described in \code{SemiParBIVProbitObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Gu C. (1992), Cross validating non-Gaussian data. \emph{Journal of Computational and Graphical Statistics}, 1(2), 169-179.

Marra G., Papageorgiou G. and Radice R. (in press), Estimation of a Semiparametric Recursive Bivariate Probit Model with Nonparametric Mixing. \emph{Australian & New Zealand Journal of Statistics}.

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian 
Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G. and Radice R. (submitted), On the Use of Asymmetric Link Functions in Copula Regression Spline Models for Bivariate Binary Responses.

Marra G., Radice R. and Missiroli S. (in press), Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Computational Statistics}.

Marra G., Radice R. and Filippou P. (submitted), Testing the Hypothesis of Exogeneity in Regression Spline Bivariate Probit Models.

Radice R., Marra G. and M. Wojtys (submitted), Copula Regression Spline Models for Binary Outcomes with Application in Health Care Utilization.

Wood S.N. (2004), Stable and efficient multiple smoothing parameter estimation for generalized additive models. \emph{Journal of the American Statistical Association}, 99(467), 673-686.

}


\section{WARNINGS}{

In the current context, convergence failure may occur when the number of 
observations is low and the data have low information content. When convergence failure is associated 
with an infinite cycling between the two steps detailed above, several practical solutions are available. For instance, one might 
either (i) lower the total number of parameters to estimate by reducing the dimension of the regression spline 
bases, (ii) set the smoothing parameters to the values obtained from the univariate fits (\code{aut.sp=FALSE}), or (iii) set the smoothing 
parameters to the values obtained from the non-converged algorithm. The default option is (iii).


} 

\seealso{
\code{\link{AT}}, \code{\link{est.prev}}, \code{\link{LM.bpm}}, \code{\link{VuongClarke.bcm}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbitObject}}, \code{\link{sem.checks}}, \code{\link{summary.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}, \code{\link{residuals.SemiParBIVProbit}}
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
u     <- rmvnorm(n,rep(0,2),Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*x1 + f2(x2) + u[,2] > 0, 1, 0)

dataSim <- data.frame(y1,y2,x1,x2,x3)

#
#

## CLASSIC BIVARIATE PROBIT

out  <- SemiParBIVProbit(y1 ~ x1 + x2 + x3, 
                         y2 ~ x1 + x2 + x3, 
                         data=dataSim)
summary(out)
AIC(out)
BIC(out)

## SEMIPARAMETRIC BIVARIATE PROBIT

## "cr" cubic regression spline basis      - "cs" shrinkage version of "cr"
## "tp" thin plate regression spline basis - "ts" shrinkage version of "tp"
## for smooths of one variable, "cr/cs" and "tp/ts" achieve similar results 
## k is the basis dimension - default is 10
## m is the order of the penalty for the specific term - default is 2

out  <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cr",k=10,m=2) + s(x3,bs="cr",k=10), 
                         y2 ~ x1 + s(x2,bs="cr",k=10)     + s(x3,bs="cr",k=10), 
                         data=dataSim)
sem.checks(out)
summary(out)
AIC(out)


## estimated smooth function plots - red lines are true curves

x2 <- sort(x2)
f1.x2 <- f1(x2)[order(x2)]-mean(f1(x2))
f2.x2 <- f2(x2)[order(x2)]-mean(f2(x2))
f3.x3 <- rep(0,length(x3))

par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(out, eq=1, select=1, scale=0); lines(x2, f1.x2, col="red")
plot(out, eq=1, select=2, scale=0); lines(x3, f3.x3, col="red")
plot(out, eq=2, select=1, scale=0); lines(x2, f2.x2, col="red")
plot(out, eq=2, select=2, scale=0); lines(x3, f3.x3, col="red")

#
## same plots but CIs 'with intercept' 

par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(out, eq=1, select=1, seWithMean=TRUE, scale=0); lines(x2, f1.x2, col="red")
plot(out, eq=1, select=2, seWithMean=TRUE, scale=0); lines(x3, f3.x3, col="red")
plot(out, eq=2, select=1, seWithMean=TRUE, scale=0); lines(x2, f2.x2, col="red")
plot(out, eq=2, select=2, seWithMean=TRUE, scale=0); lines(x3, f3.x3, col="red")


## p-values suggest to drop x3 from both equations, with a stronger 
## evidence for eq. 2. This can be also achieved using shrinkage smoothers

outSS <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cs") + s(x3,bs="cs"), 
                          y2 ~ x1 + s(x2,bs="cs") + s(x3,bs="cs"), 
                          data=dataSim)

par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(outSS, eq=1, select=1, scale=0)
plot(outSS, eq=1, select=2, ylim=c(-0.1,0.1))
plot(outSS, eq=2, select=1, scale=0)
plot(outSS, eq=2, select=2, ylim=c(-0.1,0.1))

#
#

############
## EXAMPLE 2
############
## Generate data with one endogenous variable and exclusion restriction

set.seed(0)

n <- 300

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rmvnorm(n,rep(0,2),Sigma)

cov   <- rmvnorm(n,rep(0,2),Sigma, method="svd")
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*y1 + f2(x2) + u[,2] > 0, 1, 0)

dataSim <- data.frame(y1,y2,x1,x2)

#

## Testing the hypothesis of absence of endogeneity... 

LM.bpm(y1 ~ x1 + s(x2),y2 ~ y1 + s(x2),dataSim)

# p-value suggests presence of endogeneity, hence fit a bivariate model


## CLASSIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIVProbit(y1 ~ x1 + x2, 
                        y2 ~ y1 + x2, 
                        data=dataSim)
summary(out)

## SEMIPARAMETRIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIVProbit(y1 ~ x1 + s(x2), 
                        y2 ~ y1 + s(x2), 
                        data=dataSim)
summary(out)

#

## Testing the hypothesis of absence of endogeneity post estimation... 

gt.bpm(out)

#
## average treatment effect with CIs

AT(out,eq=2,nm.bin="y1") 


## try a Clayton copula model... 

outC <- SemiParBIVProbit(y1 ~ x1 + s(x2), 
                         y2 ~ y1 + s(x2), 
                         data=dataSim, BivD="C0")
summary(outC)
AT(outC,eq=2,nm.bin="y1") 

## try a Joe copula model... 

outJ <- SemiParBIVProbit(y1 ~ x1 + s(x2), 
                         y2 ~ y1 + s(x2), 
                         data=dataSim, BivD="J0")
summary(outJ)
AT(outJ,eq=2,nm.bin="y1") 


VuongClarke.bcm(out,outJ)

#
## recursive bivariate probit modelling with unpenalized splines 
## can be achieved as follows

outFP <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cr",k=5), 
                          y2 ~ y1 + s(x2,bs="cr",k=6), fp=TRUE, 
                          data=dataSim)
summary(outFP)

#
#

############
## EXAMPLE 3
############
## Generate data with a non-random sample selection mechanism and exclusion restriction

\dontrun{ 

set.seed(0)

n <- 1100

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rmvnorm(n,rep(0,2),Sigma)

SigmaC <- matrix(0.5, 3, 3); diag(SigmaC) <- 1
cov    <- rmvnorm(n,rep(0,3),SigmaC, method="svd")
cov    <- pnorm(cov)
bi <- round(cov[,1]); x1 <- cov[,2]; x2 <- cov[,3]
  
f11 <- function(x) -0.7*(4*x + 2.5*x^2 + 0.7*sin(5*x) + cos(7.5*x))
f12 <- function(x) -0.4*( -0.3 - 1.6*x + sin(5*x))  
f21 <- function(x) 0.6*(exp(x) + sin(2.9*x)) 

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + u[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +         + u[, 2] > 0
yo <- y*(ys > 0)
  
dataSim <- data.frame(y,ys,yo,bi,x1,x2)

## Testing the hypothesis of absence of non-random sample selection... 

LM.bpm(ys ~ bi + s(x1) + s(x2),yo ~ bi + s(x1),dataSim,selection=TRUE)

# p-value suggests presence of sample selection, hence fit a bivariate model

#
## SEMIPARAMETRIC SAMPLE SELECTION BIVARIATE PROBIT
## the first equation MUST be the selection equation

out <- SemiParBIVProbit(ys ~ bi + s(x1) + s(x2), 
                        yo ~ bi + s(x1), 
                        data=dataSim, selection=TRUE)
gt.bpm(out)                        

## compare the two summary outputs
## the second output produces a summary of the results obtained when
## selection bias is not accounted for

summary(out)
summary(out$gam2)

## corrected predicted probability that 'yo' is equal to 1

est.prev(out)

## estimated smooth function plots
## the red line is the true curve
## the blue line is the naive curve not accounting for selection bias

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)]-mean(f21(x1.s))

x11()
plot(out, eq=2, ylim=c(-1.65,0.95)); lines(x1.s, f21.x1, col="red")
par(new=TRUE)
plot(out$gam2, se=FALSE, col="blue",ylim=c(-1.65,0.95),ylab="",rug=FALSE)

#
#
## try a Clayton copula model... 

outC <- SemiParBIVProbit(ys ~ bi + s(x1) + s(x2), 
                         yo ~ bi + s(x1), 
                         data=dataSim, selection=TRUE, BivD="C0")
sem.checks(outC)
summary(outC)
est.prev(outC)

#
#

############
## EXAMPLE 4
############
## Generate data with subject specific features 

set.seed(0)

n.cl <- 100
n.ob.cl <- 13

sre <- round( runif(n.cl,n.ob.cl,n.ob.cl+2) ) 
t.n <- sum(sre) 
id  <- rep(seq(1,n.cl),sre) 

sigma1 <- 1
sigma2 <- 2
rho <- rho.re <- 0.5

SigmaRE <- matrix(c(sigma1^2,            sigma1*sigma2*rho.re,
                    sigma1*sigma2*rho.re,sigma2^2             ),2,2)
                    
u1 <- rmvnorm(n.cl,rep(0,2),SigmaRE,method="svd")
u  <- cbind(rep(u1[,1],sre),rep(u1[,2],sre))


Sigma  <- matrix(rho, 2, 2); diag(Sigma) <- 1
eps    <- rmvnorm(t.n,rep(0,2),Sigma)

SigmaC <- matrix(0.5, 3, 3); diag(SigmaC) <- 1
cov    <- rmvnorm(t.n,rep(0,3),SigmaC, method="svd")

x1 <- round(pnorm(cov[,1]))
x2 <- pnorm(cov[,2])
x3 <- cov[,3]

f1 <- function(x) 0.5*(exp(x) + sin(2.9*x))
f2 <- function(x) 1.8*(0.25*exp(x) - x^3)

f1.x2 <- f1(x2) - mean(f1(x2)) 
f2.x2 <- f2(x2) - mean(f2(x2))  

y1 <- ifelse(u[,1] - 0.75*x1 + f1.x2 + x3 + eps[,1] > 0, 1, 0)
y2 <- ifelse(u[,2] + 1*x1    + f2.x2      + eps[,2] > 0, 1, 0)

dataSim <- data.frame(id,y1,y2,x1,x2,x3)


out <- SemiParBIVProbit(y1 ~ x1 + s(x2) + x3,
                        y2 ~ x1 + s(x2), 
                        data=dataSim, RE=TRUE, RE.type="NP",
                        K=3, id=dataSim$id)

summary(out)

x2s <- sort(x2)
f1.x2 <- f1.x2[order(x2)]
f2.x2 <- f2.x2[order(x2)]

par(mfrow=c(1,2))
plot(out, eq=1, scale=0); lines(x2s, f1.x2, col="red")
plot(out, eq=2, scale=0); lines(x2s, f2.x2, col="red")



############
## EXAMPLE 5
############
## Generate data for bivariate normal random effect case

set.seed(0) 

# clusters; obs per cluster; total obs

N   <- 200
rps <- round(runif(N,4,7))
n   <- sum(rps)

# cluster identifier

id <- rep(seq(1,N),rps)

# normal error terms: eps

Sigma.e <- matrix(c(1,0.5,0.5,1),2,2)
eps     <- rmvnorm(n,rep(0,2),Sigma.e)

# normal random effects: u
# s1^2, rho*s1*s2,rho*s1*s2,s2^2

s1 <- 0.5; s2 <- 1.25; rho <- -0.4

Sigma.u <- matrix(c(s1^2,rho*s1*s2,rho*s1*s2,s2^2),2,2)
u       <- rmvnorm(N,rep(0,2),Sigma.u)

u1 <- rep(u[,1],rps)
u2 <- rep(u[,2],rps)

# predictors

SigmaC <- matrix(0.5, 3, 3); diag(SigmaC) <- 1
cov    <- rmvnorm(n,rep(0,3),SigmaC, method="svd")
cov    <- pnorm(cov)

x1 <- round(runif(n))
x2 <- runif(n)
x3 <- runif(n)

# responses

y1 <- ifelse(-1.15 + 2*x1 - x2 + x3 + u1 + eps[,1] > 0, 1, 0)
y2 <- ifelse(1.55 - 1.5*x1 + 1.1*x2 + 0.5*x3 + u2 + eps[,2] > 0, 1, 0)

dataSim <- data.frame(id,y1,y2,x1,x2,x3)

out <- SemiParBIVProbit(y1 ~ x1 + x2 + x3,
                        y2 ~ x1 + x2 + x3, 
                        data=dataSim, RE=TRUE, RE.type="N",
                        id=dataSim$id)
sem.checks(out)
summary(out)

}

#
#
}

\keyword{asymmetric link function}
\keyword{non-random sample selection}
\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{nonparametric random effects} 
\keyword{endogeneity} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







