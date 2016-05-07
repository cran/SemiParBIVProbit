\name{copulaReg}
\alias{copulaReg}

\title{Semiparametric Copula Bivariate Models with Continuous Margins}

\description{ 
\code{copulaReg} can be used to fit bivariate models with continuous margins where the linear predictors of the two main equations can 
be flexibly specified using parametric and regression spline components. Several bivariate copula distributions are supported. 
During the model fitting process, the 
possible presence of associated error equations is accounted for. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines. Note that, if it makes sense, the dependence parameter of the employed bivariate 
distribution, or more generally all distribution specific parameters can be specified as functions of covariates.

}



\usage{
copulaReg(formula, data = list(), weights = NULL, subset = NULL,  
                 BivD = "N", margins = c("N","N"), 
                 gamlssfit = FALSE, fp = FALSE, infl.fac = 1, 
                 rinit = 1, rmax = 100, 
                 iterlimsp = 50, tolsp = 1e-07,
                 gc.l = FALSE, parscale, extra.regI = "t")
}



\arguments{ 

\item{formula}{In the basic setup this will be a list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms 
               are used to specify smooth smooth functions of 
                  predictors. For the case of more than two equations see the example below and the documentation 
                  of \code{SemiParBIVProbit()} for more details.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{copulaReg} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}


\item{margins}{It indicates the distributions used for the two margins. Possible choices are: normal ("N"), log-normal ("LN"), 
Gumbel ("GU"), reverse Gumbel ("rGU"), 
               logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), 
               gamma with identity link for the location parameter ("GAi"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic distribution).}
               
\item{gamlssfit}{If \code{gamlssfit = TRUE} then gamlss models are fitted. This is useful for obtaining 
                 starting values, for instance.}               

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


}


\value{
The function returns an object of class \code{copulaReg} as described in \code{copulaRegObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G. and Radice R. (submitted), A Bivariate Copula Additive Model for Location, Scale and Shape.

}



\section{WARNINGS}{

Convergence failure may sometimes occur. Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be 0 and the latter positive definite.
\code{copulaReg()} will produce some warnings if there is a convergence issue.

In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). Using \code{gamlssfit = TRUE} is typically more effective than the first two options as
this will provide better calibrated starting values as compared to those obtained from the default starting value procedure.
The default option is, however, \code{gamlssfit = FALSE} only because it tends to be computationally cheaper and because the 
default starting value procedure has typically been found to do a satisfactory job in most cases. 
(The results obtained when using 
\code{gamlssfit = FALSE} and \code{gamlssfit = TRUE} could also be compared to check if starting values make any difference.)

The above suggestions may help, especially the latter option. However, the user should also consider
re-specifying the model, and/or using a diferrent dependence structure and/or checking that the chosen marginal 
distributions fit the responses.
In our experience, we found that convergence failure typically occurs 
when the model has been misspecified and/or the sample size is low compared to the complexity of the model. Examples
of misspecification include using a Clayton copula rotated by 90 degrees when a positive
association between the margins is present instead, using marginal distributions that do not fit
the responses, and 
employing a copula which does not accommodate the type and/or strength of
the dependence between the margins (e.g., using AMH when the association between the margins is strong).
It is also worth bearing in mind that the use of three parameter marginal distributions requires the data
to be more informative than a situation in which two parameter distributions are used instead.

When comparing competing models (for instance, by keeping the linear predictor specifications
fixed and changing the copula), if the computing time for a set of alternatives
is considerably higher than that of another set then it may mean that 
the more computationally demanding models are not able to fit the data very well (as a higher number of 
iterations is required to reach convergence). As a practical check, this may be verified by
fitting all competing models and, provided convergence is achieved, comparing their respective AIC and BICs, for instance.

} 

\seealso{
\code{\link{adjCov}}, \code{\link{VuongClarke}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{copulaRegObject}}, \code{\link{conv.check}}, \code{\link{summary.copulaReg}}, \code{\link{predict.SemiParBIVProbit}}
}

\examples{

library(SemiParBIVProbit)

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

y1 <- -1.55 + 2*x1    + f1(x2) + u[,1]
y2 <- -0.25 - 1.25*x1 + f2(x2) + u[,2]

dataSim <- data.frame(y1, y2, x1, x2, x3)

resp.check(y1, "N")
resp.check(y2, "N")

eq.mu.1     <- y1 ~ x1 + s(x2) + s(x3)
eq.mu.2     <- y2 ~ x1 + s(x2) + s(x3)
eq.sigma2.1 <-    ~ 1
eq.sigma2.2 <-    ~ 1
eq.theta    <-    ~ x1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma2.1, eq.sigma2.2, eq.theta)

# the order above is the one to follow when
# using more than two equations

out  <- copulaReg(fl, data = dataSim)
conv.check(out)
post.check(out)
summary(out)
AIC(out)
BIC(out)
jc.probs(out, 1.4, 2.3, intervals = TRUE)[1:4,]


}

\keyword{semiparametric bivariate modelling} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}
\keyword{copula}


