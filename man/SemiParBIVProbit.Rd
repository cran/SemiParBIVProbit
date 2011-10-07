\name{SemiParBIVProbit}
\alias{SemiParBIVProbit}

\title{Semiparametric Bivariate Probit Modelling}

\description{ 
\code{SemiParBIVProbit} can be used to fit bivariate probit models where the linear predictors are flexibly specified using parametric and 
regression spline components. During the model fitting process, the possible presence of correlated error equations, endogeneity or sample selection is accounted for. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines (isotropic). The current implementation does not support scale invariant tensor 
product smooths. 

}



\usage{
SemiParBIVProbit(formula.eq1, formula.eq2, data=list(),  
                 selection=FALSE, iterlimFS=1, iterlimSP=50, pr.tol=1e-6, 
                 gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, 
                 rinit=1, rmax=100, fterm=sqrt(.Machine$double.eps), 
                 mterm=sqrt(.Machine$double.eps),   
                 control=list(maxit=50,tol=1e-6,step.half=25,
                              rank.tol=.Machine$double.eps^0.5))
}


\arguments{ 

\item{formula.eq1}{A GAM formula for equation 1. \code{s} terms are used to specify smooth smooth functions of 
                  predictors. Notice that \code{SemiParBIVProbit} supports the use shrinkage smoothers for variable
                  selection purposes. See the examples below and the documentation of \code{mgcv} for further details on GAM formula specifications.} 

\item{formula.eq2}{A GAM formula for equation 2.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParBIVProbit} is called.} 

\item{selection}{If \code{TRUE}, then the numerical routine for bivariate probit modelling in the presence of sample selection is employed.}

\item{iterlimFS}{A positive integer specifying the maximum number of iterations of the Fisher scoring algorithm to allow in the smoothing parameter estimation step.} 

\item{iterlimSP}{A positive integer specifying the maximum number of loops to be performed before the smoothing parameter estimation step is terminated.} 

\item{pr.tol}{Tolerance to use in judging convergence of the algorithm when automatic smoothing parameter selection is used.} 

\item{gamma}{It is is an inflation factor for the model degrees of freedom in the UBRE score. Smoother models can be obtained setting 
             this parameter to a value greater than 1. Typically \code{gamma=1.4} achieves this.}

\item{aut.sp}{If \code{TRUE}, then automatic multiple smoothing parameter selection is carried out. If \code{FALSE}, then smoothing parameters 
              are set to the values obtained from the univariate fits.}

\item{fp}{If \code{TRUE}, then a fully parametric model with regression splines if fitted. This only makes sense 
          if used jointly with \code{aut.sp=FALSE}. See the example below.}  

\item{start.v}{Although not recommended, starting values for the parameters of the two equations and correlation coefficient 
               can be provided here.}  

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds. Bad initial 
             values waste a few steps while the radius is adjusted, but do not prevent the algorithm from working properly. See the documentation
             of \code{trust} for further details.}

\item{rmax}{Maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest 
            descent path.}        
        
\item{fterm}{Positive scalar giving the tolerance at which the difference in objective function values in a step is considered close 
             enough to zero to terminate the algorithm.}

\item{mterm}{Positive scalar giving the tolerance at which the two-term Taylor-series approximation to the difference in objective 
             function values in a step is considered close enough to zero to terminate the algorithm.}

\item{control}{It is a list containing iteration control constants with the following elements: \code{maxit}: maximum number of iterations of the 
               \code{magic} algorithm; \code{tol}: tolerance to use in judging convergence; \code{step.half}: if a trial step fails then 
               the method tries halving it up to a maximum of \code{step.half} times; \code{rank.tol}: constant used to test for numerical rank 
               deficiency of the problem. See the documentation of \code{magic} in \code{mgcv} for further details.}

}


\details{

The bivariate probit model has a probit link for each of the two equations, and models the association
between the responses by the correlation parameter \eqn{\rho}{\rho} of a standard bivariate normal distribution. In a semiparametric bivariate probit model the linear predictors are flexibly specified using 
parametric components and smooth functions of covariates. Replacing the smooth components with their regression spline expressions yields a fully parametric bivariate probit model. In principle, classic 
maximum likelihood estimation can be employed. However, to avoid overfitting, penalized likelihood maximization has to be employed instead. Here the use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support from the data. The trade-off between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimize the approximate Un-Biased Risk Estimator (UBRE).

The optimization problem is solved by Fisher scoring's method. Automatic smoothing parameter selection is integrated 
using a performance-oriented iteration approach (Gu, 1992; Wood, 2004). At each iteration, (i) the penalized weighted least squares 
problem is solved, then (ii) the smoothing parameters of that problem estimated by approximate UBRE. Steps (i) and (ii) are iterated until 
convergence. Details of the underlying fitting methods are given in Marra and Radice (2011). 

}


\value{
The function returns an object of class \code{SemiParBIVProbit} as described in \code{SemiParBIVProbitObject}.
}


\author{
Giampiero Marra, University College London

Rosalba Radice, London School of Hygiene and Tropical Medicine

}



\references{

Gu C. (1992), Cross validating non-Gaussian data. \emph{Journal of Computational and Graphical Statistics}, 1, 169-179.

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian 
Journal of Statistics}, 39, 259-279.

Marra G. and Radice R., A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modelling. \emph{Submitted}.

Wood S.N. (2004), Stable and efficient multiple smoothing parameter estimation for generalized additive models. \emph{Journal of the American Statistical Association}, 99, 673-686.

}


\section{WARNINGS}{

Any automatic smoothing parameter selection procedure is not likely to work well when the data have low information content. In bivariate probit models, this 
issue is especially relevant when \eqn{\rho}{\rho} is high and the number of observations low. Here, convergence failure is
typically associated with an infinite cycling between the two steps detailed above. If this occurs, as some practical solutions, one might either (i) lower the total number of parameters to estimate by reducing the dimension of the regression spline 
bases, (ii) set the smoothing parameters to the values obtained from the univariate fits (\code{aut.sp=FALSE}), or (iii) set the smoothing parameters to the values obtained from
the non-converged algorithm. The default option is (iii).

Fully parametric modelling is allowed for. However, it is not possible to specify one linear predictor as a function of parametric and smooth
components, and the other as a function of parametric terms only. If continuous covariates are available, then we should let the data determine which
effects are linear or non-linear and for which equations.


} 

\seealso{
\code{\link{AT}}, \code{\link{InfCr}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbitObject}}, \code{\link{summary.SemiParBIVProbit}}, \code{\link{residuals.SemiParBIVProbit}}
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

Sigma <- matrix(c(1,0.5,0.5,1),2,2)
u     <- rmvnorm(n,rep(0,2),Sigma)

x1 <- round(runif(n))
x2 <- runif(n)
x3 <- runif(n)

f1   <- function(x) (cos(pi*2*x)) + sin(pi*x)
f2   <- function(x) (x+exp(-30*(x-0.5)^2))   

y1 <- rep(0,n)
y1 <- replace(y1, -1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1)
y2 <- rep(0,n)
y2 <- replace(y2, -0.25 - 1.25*x1 + f2(x2) + u[,2] > 0, 1)

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
                         y2 ~ x1 + s(x2,bs="cr",k=10) + s(x3,bs="cr",k=10), 
                         data=dataSim)
summary(out)
InfCr(out)


## estimated smooth function plots - red lines are true curves

x2 <- sort(x2)
f1.x2 <- f1(x2)[order(x2)]-mean(f1(x2))
f2.x2 <- f2(x2)[order(x2)]-mean(f2(x2))
f3.x3 <- rep(0,length(x3))

par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(out, eq=1, select=1); lines(x2, f1.x2, col="red")
plot(out, eq=1, select=2); lines(x3, f3.x3, col="red")
plot(out, eq=2, select=1); lines(x2, f2.x2, col="red")
plot(out, eq=2, select=2); lines(x3, f3.x3, col="red")

#
## same plots but CIs 'with intercept' 

x11()
par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(out, eq=1, select=1, seWithMean=TRUE); lines(x2, f1.x2, col="red")
plot(out, eq=1, select=2, seWithMean=TRUE); lines(x3, f3.x3, col="red")
plot(out, eq=2, select=1, seWithMean=TRUE); lines(x2, f2.x2, col="red")
plot(out, eq=2, select=2, seWithMean=TRUE); lines(x3, f3.x3, col="red")


#
## set gamma=1.4 to obtain a smoother model

outS <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cr",k=10) + s(x3,bs="cr",k=10), 
                         y2 ~ x1 + s(x2,bs="cr",k=10) + s(x3,bs="cr",k=10), 
                         data=dataSim, gamma=1.4)
summary(outS)

## compare the two summary outputs

## p-values suggest to drop x3 from both equations, with a stronger 
## evidence for eq. 2. This can be also achieved via shrinkage smoothers

outSS <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cs",k=10) + s(x3,bs="cs",k=10), 
                          y2 ~ x1 + s(x2,bs="cs",k=10) + s(x3,bs="cs",k=10), 
                          data=dataSim, gamma=1.4)

par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(outSS, eq=1, select=1)
plot(outSS, eq=1, select=2, ylim=c(-0.1,0.1))
plot(outSS, eq=2, select=1)
plot(outSS, eq=2, select=2, ylim=c(-0.1,0.1))

#
#

############
## EXAMPLE 2
############
## Generate data with one endogenous variable and exclusion restriction

set.seed(0)

n <- 400

Sigma <- matrix(c(1,0.5,0.5,1),2,2)
u     <- rmvnorm(n,rep(0,2),Sigma)

x1 <- round(runif(n))
x2 <- runif(n)

f1   <- function(x) (cos(pi*2*x)) + sin(pi*x)
f2   <- function(x) (x+exp(-30*(x-0.5)^2))   

y1 <- rep(0,n)
y1 <- replace(y1, -1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1)
y2 <- rep(0,n)
y2 <- replace(y2, -0.25 - 1.25*y1 + f2(x2) + u[,2] > 0, 1)

dataSim <- data.frame(y1,y2,x1,x2)

#
#

## CLASSIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIVProbit(y1 ~ x1 + x2, 
                        y2 ~ y1 + x2, 
                        data=dataSim)
summary(out)

## SEMIPARAMETRIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cr",k=10), 
                        y2 ~ y1 + s(x2,bs="cr",k=10), 
                        data=dataSim)
summary(out)

AT(out,eq=2,nm.bin="y1") ## average treatment effect (ATE) of y1 with CIs
AT(out,eq=2,nm.bin="y1",E=FALSE) ## average treatment effect on the treated (ATT)

#
## recursive bivariate probit modelling with unpenalized splines 
## can be achieved as follows

outFP <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cr",k=5), 
                          y2 ~ y1 + s(x2,bs="cr",k=6), aut.sp=FALSE, fp=TRUE, 
                          data=dataSim)
summary(outFP)


#
## example using a 2D smooth 

f.biv <- function(x,z) {-0.7*exp(((3*x+3) + 0.7*(3*z-3)^2)/5)}

y1 <- rep(0,n)
y1 <- replace(y1, -1.55 + 2*x1    + f1(x2)     + u[,1] > 0, 1)
y2 <- rep(0,n)
y2 <- replace(y2, 3.5 - 1.25*y1 + f.biv(x2,x3) + u[,2] > 0, 1)

dataSim <- data.frame(y1,y2,x1,x2,x3)


outSb <- SemiParBIVProbit(y1 ~ x1 + s(x2,bs="cr"), 
                          y2 ~ y1 + s(x2,x3,bs="tp"), 
                          data=dataSim)
summary(outSb)

x11()
par(mfrow=c(1,2),mar=c(4.5,4.5,2,2))

x2 <- x3 <- seq(0, 1, length=40)
z <- outer(x2, x3, f.biv)

persp(x2, x3, z, theta = 30, phi = 25,
      xlab="x2", ylab="x3", zlab="f.biv")
plot(outSb, eq=2, select=1,theta = 30, phi = 25) 

#
#

############
## EXAMPLE 3
############
## Generate data with sample selection mechanism and exclusion restriction

set.seed(0)

n <- 2000

Sigma <- matrix(c(1,0.5,0.5,1),2,2)
u     <- rmvnorm(n,rep(0,2),Sigma)

SigmaC <- matrix( c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1), 3 , 3)

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


## SEMIPARAMETRIC SAMPLE SELECTION BIVARIATE PROBIT
## the first equation MUST be the selection equation

out <- SemiParBIVProbit(ys ~ bi + s(x1,bs="cr",k=10) + s(x2,bs="cr",k=10), 
                        yo ~ bi + s(x1,bs="cr",k=10), 
                        data=dataSim, selection=TRUE)

AT(out,eq=2,nm.bin="bi",E=FALSE) ## average treatment effect on the treated (ATT) 

## compare the two summary outputs
## the second output produces a summary of the results obtained when only 
## the outcome equation is fitted, i.e. selection bias is not accounted for

summary(out)
summary(out$gam2)

## mean predicted probabilities that 'yo' is equal to 1
## the second predicted probability does not account for selection bias

mean(round(pnorm(out$eta2[out$eta2!=0])))
mean(round(pnorm(predict(out$gam2))))

## estimated smooth function plots
## the red line is the true curve
## the blue line is the naive curve not accounting for selection bias

x1 <- sort(x1)
f21.x1 <- f21(x1)[order(x1)]-mean(f21(x1))

plot(out, eq=2, select=1,ylim=c(-1.2,1)); lines(x1, f21.x1, col="red")
par(new=TRUE)
plot(out$gam2, select=1, se=FALSE, col="blue",ylim=c(-1.2,1),ylab="",rug=FALSE)

#
#

}


\keyword{sample selection}
\keyword{semiparametric bivariate probit modelling} 
\keyword{bivariate probit model} 
\keyword{endogeneity} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







