\name{SemiParTRIVProbit}
\alias{SemiParTRIVProbit}

\title{Semiparametric Trivariate Probit Regression}

\description{ 
\code{SemiParTRIVProbit} can be used to fit trivariate binary models where the linear predictors of the two main equations can 
be flexibly specified using 
parametric and regression spline components. 

}



\usage{
SemiParTRIVProbit(formula, data = list(), weights = NULL, subset = NULL,  
                 Model = "T", penCor = "unpen", sp.penCor = 3,
                 fp = FALSE, infl.fac = 1, 
                 rinit = 1, rmax = 100, 
                 iterlimsp = 50, tolsp = 1e-07,
                 gc.l = FALSE, parscale, extra.regI = "t")
}


\arguments{ 

\item{formula}{In the basic setup this will be a list of three formulas. \code{s} terms 
               are used to specify smooth smooth functions of 
                  predictors. \code{SemiParBIVProbit} supports the use shrinkage smoothers for variable
                  selection purposes and more. See the examples below and the documentation of \code{SemiParBIVProbit} for further 
                  details.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParBIVProbit} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{Model}{It indicates the type of model to be used in the analysis. Possible values are "T" (trivariate model), 
             "TSS" (trivariate model with double sample selection).}

\item{penCor}{Type of penalty for correlation coefficients. Possible values are "unpen", "lasso", "ridge".}

\item{sp.penCor}{Starting value for smoothing parameter of \code{penCor}.}
                    
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

This function fits trivariate models.


For sample selection models, if there are factors in the model, before fitting, the user has to ensure 
that the numbers of factor variables' levels in the selected sample 
are the same as those in the complete dataset. Even if a model could be fitted in such a situation,
the model may produce fits which are 
not coherent with the nature of the correction sought. For more details see \code{?SemiParBIVProbit}. 


}


\value{
The function returns an object of class \code{SemiParTRIVProbit} as described in \code{SemiParTRIVProbitObject}.
}


\author{

Authors: Panagiota Filippou and Giampiero Marra

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\section{WARNINGS}{

Convergence failure may sometimes occur. Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be 0 and the latter positive definite.
\code{SemiParTRIVProbit()} will produce some warnings if there is a convergence issue.

In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). Penalising the correlations using argument \code{penCor} may help a lot in some situations as
in our experience in hard situations the correlation coefficients are typically the most difficult to estimate. 

The above suggestions may help, especially the latter option. However, the user should also consider
looking into the proportions of 1 and 0 available for each event of the trivariate model. It may 
be the case that certain events do not have many observations associated to them 
which will make estimation very difficult.
As usual, model complexity plays a role.

} 

\seealso{
\code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParTRIVProbitObject}}, \code{\link{conv.check}}, \code{\link{summary.SemiParTRIVProbit}}, \code{\link{predict.SemiParBIVProbit}}
}

\examples{

\dontrun{ 

library(SemiParBIVProbit)

############
## EXAMPLE 1
############
## Generate data
## Correlation between the two equations 0.5 - Sample size 400 

set.seed(0)

n <- 400

Sigma <- matrix(0.5, 3, 3); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,3), Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

y1 <- ifelse(-1.55 +    2*x1 - 0.6*x2 + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*x1 +     x2 + u[,2] > 0, 1, 0)
y3 <- ifelse(-0.75 + 0.25*x1 +        + u[,3] > 0, 1, 0)

dataSim <- data.frame(y1, y2, y3, x1, x2, x3)

out  <- SemiParTRIVProbit(list(y1 ~ x1 + x2, 
                               y2 ~ x1 + x2,
                               y3 ~ x1), 
                         data = dataSim)
conv.check(out)
summary(out)
AIC(out)
BIC(out)


############
## EXAMPLE 2
############
## Generate data
## with double sample selection

set.seed(0)

n <- 5000

Sigma <- matrix(c(1,   0.5, 0.4,
                  0.5,   1, 0.6,
                  0.4, 0.6,   1 ), 3, 3)

u <- rMVN(n, rep(0,3), Sigma)

x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- runif(n)
  
y1 <-  1    + 1.5*x1 -     x2 + 0.8*x3 - 1.5*x4 + u[, 1] > 0
y2 <-  1    - 2.5*x1 + 1.2*x2 +     x3          + u[, 2] > 0
y3 <-  1.58 + 1.5*x1 - 2.5*x2                   + u[, 3] > 0

dataSim <- data.frame(y1, y2, y3, x1, x2, x3, x4)

f.l <- list(y1 ~ x1 + x2 + x3 + x4, 
            y2 ~ x1 + x2 + x3, 
            y3 ~ x1 + x2)   
          
out <- SemiParTRIVProbit(f.l, data = dataSim, Model = "TSS")

prev(out)
prev(out, type = "univariate")
prev(out, type = "naive")

}

}

\keyword{semiparametric trivariate probit modelling}  
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







