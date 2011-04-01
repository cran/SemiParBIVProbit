\name{residuals.SemiParBIVProbit}
\alias{residuals.SemiParBIVProbit}

\title{Residuals}

\description{ 
\code{residuals} efficiently calculates Pearson and working residuals.}



\usage{


\method{residuals}{SemiParBIVProbit}(object,...)

}

\arguments{ 

\item{object}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{...}{Other arguments.}
 

}


\details{

Justification for these residuals is from the penalized IRLS algorithm used to fit the model (Marra and Radice, in press). 

}


\value{

\item{r.p}{It contains three columns corresponding to the Pearson residuals from the fitted model. The first two columns are of interest; they 
           refer to the two linear predictors in the bivariate probit model.} 

\item{r.w}{It contains three columns corresponding to the working residuals from the fitted model.} 



}



\author{

Giampiero Marra, University College London

Rosalba Radice, London School of Hygiene and Tropical Medicine

}



\references{

Marra G. and Radice R. (in press), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{The Canadian 
Journal of Statistics}.

}



\section{WARNINGS}{

Unfortunately residual plots from models fitted to binary data do not contain much information about the goodness-of-fit of the model. For instance, if one wants to check distributional assumptions from residuals, QQ-plots and/or half-normal plots with simulated envelope provide no useful information (this has also been checked in simulation). For binary data, it is necessary to be able to group the residuals according to some criterion.     

} 




\seealso{
\code{\link{InfCr}}, \code{\link{SemiParBIVProbit}}, \code{\link{plot.SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}  
}



\keyword{Pearson residuals} 
\keyword{working residuals} 









