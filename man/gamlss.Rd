\name{gamlss}
\alias{gamlss}

\title{Generalised Additive Models for Location, Scale and Shape}

\description{ 
\code{gamlss} fits flexible univariate regression models with several continuous and discrete distributions, and types of covariate 
effects. The purpose of this function was only to provide, in some cases, starting values 
for the simultaneous models in the package, but it has now been made available in the form of a proper function should the user wish to fit 
univariate models using the general estimation approach of this package. The distributions implemented here 
have been parametrised according to Rigby and Stasinopoulos (2005).
 
}



\usage{
gamlss(formula, data = list(), weights = NULL, subset = NULL,  
       margin = "N", robust = FALSE, rc = 3, rsim = 100, infl.fac = 1, 
       rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
       gc.l = FALSE, parscale, extra.regI = "t")
}



\arguments{ 

\item{formula}{List of equations. This should contain one or more equations.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{gamlss} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}


\item{margin}{Possible distributions are normal ("N"), normal where 
               sigma2 corresponds to 
               the standard deviation instead of the variance ("N2"),  
               log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), 
               logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), 
               gamma with identity link for the location parameter ("GAi"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic distribution), Poisson ("PO"), zero truncated 
               Poisson ("ZTP"), negative binomial - type I ("NBI"), negative 
               binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG").}
 
\item{robust}{If \code{TRUE} then the robust version of the model is fitted.}
             
\item{rc}{Robust constant.}       

\item{rsim}{Number of simulation draws to use in robust estimation.}  
 
\item{infl.fac}{Inflation factor for the model degrees of freedom in the approximate AIC. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds.}

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

The underlying algorithm is based on the procedure used for \code{SemiParBIVProbit()}. For more details
see ?SemiParBIVProbit. 

There are many continuous/discrete distributions to choose from and we plan to include more 
options. Get in touch if you are interested in a particular distribution.


}


\value{
The function returns an object of class \code{gamlss} as described in \code{gamlssObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G. and Radice R. (submitted), Bivariate Copula Additive Models for Location, Scale and Shape.

Rigby R.A. and Stasinopoulos D.M. (2005). Generalized additive models for location, scale and shape (with discussion). \emph{Journal of the Royal Statistical Society, Series C}, 54(3), 507-554.

}



\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{gamlss()} will produce some warnings if there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model. In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). However, the user should especially consider
re-specifying/simplifying the model, and/or checking that the chosen distribution fits the response well.
In our experience, we found that convergence failure typically occurs 
when the model has been misspecified and/or the sample size is low compared to the complexity of the model. 
It is also worth bearing in mind that the use of three parameter distributions requires the data
to be more informative than a situation in which two parameter distributions are used instead.


} 

\seealso{
\code{\link{plot.SemiParBIVProbit}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{gamlssObject}}, \code{\link{conv.check}}, \code{\link{summary.gamlss}}, \code{\link{predict.SemiParBIVProbit}}
}

\examples{

\dontrun{ 

library(SemiParBIVProbit)

set.seed(0)

n <- 400

x1 <- round(runif(n))
x2 <- runif(n)
x3 <- runif(n)
f1 <- function(x) cos(pi*2*x) + sin(pi*x)
y1 <- -1.55 + 2*x1 + f1(x2) + rnorm(n)

dataSim <- data.frame(y1, x1, x2, x3)
resp.check(y1, "N")

eq.mu <- y1 ~ x1 + s(x2) + s(x3)
eq.s2 <-    ~ s(x3)
fl    <- list(eq.mu, eq.s2)

out <- gamlss(fl, data = dataSim)

conv.check(out)
post.check(out)

plot(out, eq = 1, scale = 0, pages = 1, seWithMean = TRUE)
plot(out, eq = 2, seWithMean = TRUE)

summary(out)

AIC(out)
BIC(out)
}

}

\keyword{flexible regression modelling} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{distribution}


