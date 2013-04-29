\name{SemiParBIVProbit}
\alias{SemiParBIVProbit}

\title{Semiparametric Bivariate Probit Modelling}

\description{ 
\code{SemiParBIVProbit} can be used to fit bivariate probit models where the linear predictors can be flexibly specified using parametric and 
regression spline components, and parametric or nonparametric random effects. During the model fitting process, the possible presence of 
correlated error equations, endogeneity or non-random sample selection is accounted for. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines (isotropic). The current implementation does not support scale invariant tensor 
product smooths. 

}



\usage{
SemiParBIVProbit(formula.eq1, formula.eq2, data=list(), weights=NULL,  
                 start.v=NULL, selection=FALSE, H.pm=FALSE, gamma=1, 
                 aut.sp=TRUE, fp=FALSE, RE=FALSE, RE.type="N", NGQ=10, 
                 K=2, id=NULL, e.npRE=TRUE, rinit=1, rmax=100, 
                 fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
                 iterlimsp=50, pr.tolsp=1e-6,    
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

\item{start.v}{Starting values for the parameters of the two equations and correlation coefficient 
               can be provided here.} 

\item{selection}{If \code{TRUE}, then the numerical routine for bivariate probit modelling in the presence of non-random sample selection is employed.}

\item{H.pm}{If \code{TRUE}, then the Hessian (rather than the Fisher) information matrix is used for bivariate probit modelling. This option 
         is valid ONLY for parametric model specifications.}

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

The bivariate probit model has a probit link for each of the two equations, and models the association
between the responses through the correlation parameter \eqn{\rho}{\rho} of a standardised bivariate normal distribution. In a 
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
convergence. Details of the underlying fitting methods are given in Marra and Radice (2011,submitted). 

The use of parametric and nonparametric random effects is allowed for via the option \code{RE}. The nonparametric case
is illustrated in Example 4 (Marra et al., in press), whereas the parametric one, in Example 5, will be completed soon. 

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

Marra G. and Radice R. (submitted), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, invited revision.

Wood S.N. (2004), Stable and efficient multiple smoothing parameter estimation for generalized additive models. \emph{Journal of the American Statistical Association}, 99(467), 673-686.

}


\section{WARNINGS}{

In the current context, convergence failure may occur when \eqn{\rho}{\rho} is very high, the number of 
observations is low, and the data have low information content. When convergence failure is associated 
with an infinite cycling between the two steps detailed above, several practical solutions are available. For instance, one might 
either (i) lower the total number of parameters to estimate by reducing the dimension of the regression spline 
bases, (ii) set the smoothing parameters to the values obtained from the univariate fits (\code{aut.sp=FALSE}), or (iii) set the smoothing 
parameters to the values obtained from the non-converged algorithm. The default option is (iii).


} 

\seealso{
\code{\link{AT}}, \code{\link{InfCr}}, \code{\link{LM.bpm}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbitObject}}, \code{\link{sem.checks}}, \code{\link{summary.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}, \code{\link{residuals.SemiParBIVProbit}}
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

f1   <- function(x) (cos(pi*2*x)) + sin(pi*x)
f2   <- function(x) (x+exp(-30*(x-0.5)^2))   

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
InfCr(out)
InfCr(out,cr="BIC")

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
InfCr(out)


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

n <- 400

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rmvnorm(n,rep(0,2),Sigma)

cov   <- rmvnorm(n,rep(0,2),Sigma, method="svd")
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) (cos(pi*2*x)) + sin(pi*x)
f2   <- function(x) (x+exp(-30*(x-0.5)^2))   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*y1 + f2(x2) + u[,2] > 0, 1, 0)

dataSim <- data.frame(y1,y2,x1,x2)

#

## Testing the hypothesis of absence of endogeneity... 

LM.bpm(y1 ~ x1 + s(x2),y2 ~ y1 + s(x2),dataSim)

# p-value suggests presence of endogeneity, hence fit a bivariate model


## CLASSIC RECURSIVE BIVARIATE PROBIT

# Algorithm based on Fisher Scoring
out <- SemiParBIVProbit(y1 ~ x1 + x2, 
                        y2 ~ y1 + x2, 
                        data=dataSim)
summary(out)

# Algorithm based on Newton's method
out <- SemiParBIVProbit(y1 ~ x1 + x2, 
                        y2 ~ y1 + x2, 
                        data=dataSim,H.pm=TRUE)
summary(out)


## SEMIPARAMETRIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIVProbit(y1 ~ x1 + s(x2), 
                        y2 ~ y1 + s(x2), 
                        data=dataSim)
summary(out)

#
## average treatment effect on the treated (ATT) with CIs

AT(out,eq=2,nm.bin="y1",E=FALSE) 

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

n <- 1200

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

## compare the two summary outputs
## the second output produces a summary of the results obtained when
## selection bias is not accounted for

summary(out)
summary(out$gam2)

## corrected predicted probability that 'yo' is equal to 1
## ob.pr   is the observed proportion of 1's
## pred.pr is the mean predicted probability of observing 1, 
##                       for the sample of non-respondents 

## simulate from posterior to obtain standard 
## deviation for pred.pr (this is similar to a parametric 
## bootstrap procedure)

  pr.prS    <- NA
  eta1Sim   <- out$eta1S
  eta2Sim   <- out$eta2S 
  athrhoSim <- out$athrhoS

for(i in 1:dim(eta1Sim)[2]){

  p0  <- pnorm(-eta1Sim[,i]) 
  p01 <- pnorm2(-eta1Sim[,i],eta2Sim[,i],cov12=-tanh(athrhoSim[i]))
  pr.prS[i] <- mean( (p01/p0)[out$y1==0] )

}
  
  pred.pr <- mean( (out$p01/out$p0)[out$y1==0] )
  var.prS <- var(pr.prS)
  ob.pr   <- mean( out$y2[out$y1==1] ) 
  var.ob  <- (ob.pr*(1-ob.pr))/out$n.sel
  we      <- c( out$n.sel,(out$n-out$n.sel) )/out$n
 
  wm <- ob.pr*we[1] + pred.pr*we[2]
  sv <- sqrt(var.ob*we[1]^2 + var.prS*we[2]^2)

  qz <- qnorm(0.05/2, lower.tail = FALSE)
  lb <- wm - qz*sv # lower bound
  ub <- wm + qz*sv # upper bound 

  round(c(lb,wm,ub),2)


## estimated smooth function plots
## the red line is the true curve
## the blue line is the naive curve not accounting for selection bias

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)]-mean(f21(x1.s))

x11()
plot(out, eq=2, ylim=c(-1.2,1)); lines(x1.s, f21.x1, col="red")
par(new=TRUE)
plot(out$gam2, se=FALSE, col="blue",ylim=c(-1.2,1),ylab="",rug=FALSE)

#
#

############
## EXAMPLE 4
############
## Generate data with one endogenous variable, exclusion restriction
## and subject specific features 

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

y1 <- ifelse(u[,1] +        -0.75*x1 + f1.x2 + x3 + eps[,1] > 0, 1, 0)
y2 <- ifelse(u[,2] +  -1.5*y1 + 1*x1 + f2.x2      + eps[,2] > 0, 1, 0)

dataSim <- data.frame(id,y1,y2,x1,x2,x3)


out <- SemiParBIVProbit(y1 ~      x1 + s(x2) + x3,
                        y2 ~ y1 + x1 + s(x2), 
                        data=dataSim, RE=TRUE, RE.type="NP",
                        K=3, id=dataSim$id)

summary(out)
AT(out,eq=2,nm.bin="y1",E=FALSE)

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

N   <- 500
rps <-round(runif(N,4,7))
n   <- sum(rps)

# cluster identifier

id <- rep(seq(1,N),rps)

# normal error terms: eps

Sigma.e <- matrix(c(1,0.5,0.5,1),2,2)
eps     <- rmvnorm(n,rep(0,2),Sigma.e)

# normal random effects: u
# s1^2, rho*s1*s2,rho*s1*s2,s2^2

s1 <- 0.3; s2 <- 1.5; rho <- -0.5
Sigma.u <- matrix(c(s1^2,rho*s1*s2,rho*s1*s2,s2^2),2,2)
u       <- rmvnorm(N,rep(0,2),Sigma.u)

u1 <- rep(u[,1],rps)
u2 <- rep(u[,2],rps)

# predictors

x1 <- runif(n)-0.5

# responses

y1 <- rep(0,n)
y2 <- rep(0,n)

y1 <- replace(y1, u1 + 0.3*x1 + eps[,1] > 0, 1)
y2 <- replace(y2, u2 + 1.3*x1 + eps[,2] > 0, 1)

dataSim <- data.frame(id,y1,y2,x1)

out <- SemiParBIVProbit(y1 ~ x1,
                        y2 ~ x1, 
                        data=dataSim, RE=TRUE, RE.type="N",
                        id=dataSim$id)
sem.checks(out)
summary(out)

}

#
#
}


\keyword{non-random sample selection}
\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{nonparametric random effects} 
\keyword{endogeneity} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







