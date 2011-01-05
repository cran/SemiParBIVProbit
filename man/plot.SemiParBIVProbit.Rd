\name{plot.SemiParBIVProbit}
\alias{plot.SemiParBIVProbit}

\title{SemiParBIVProbit plotting}


\description{It takes a fitted \code{SemiParBIVProbit} object produced by \code{SemiParBIVProbit()} and plots the 
component smooth functions that make it up on the scale of the linear predictor. 
}

\usage{
\method{plot}{SemiParBIVProbit}(x, eq, select, rug=TRUE, se=TRUE, se.l=1.95996, seWithMean=FALSE, 
     n=100, xlab = NULL, ylab=NULL, xlim=NULL, ylim = NULL, main=NULL, 
     trans = I, n2 = 40, theta = 30, phi = 30, too.far = 0.1, ...)
}



\arguments{ 
  
\item{x}{A fitted \code{SemiParBIVProbit} object as produced by \code{SemiParBIVProbit()}.}

\item{eq}{The equation from which smooth terms should be considered for printing.}

\item{rug}{If \code{TRUE} (default), then the covariate to which the plot applies is displayed as a rug plot
 at the foot of each plot of a 1-d smooth.} 
  
\item{se}{If \code{TRUE} (default), then Bayesian `confidence' intervals are added to 1-d plots at \code{se.l} standard errors above and below the estimate 
of the smooth being plotted. }

\item{se.l}{Significance level.}

\item{select}{The smooth component to be considered for printing.}

\item{n}{Number of points used for each 1-d plot. As a rule of thumb, this number has to be several times the estimated 
degrees of freedom for the smooth.}

\item{n2}{Square root of number of points used to grid estimates of 2-d functions.}

\item{theta}{One of the perspective plot angles.}

\item{phi}{The other perspective plot angle.}

\item{xlab}{If supplied, then this is used as the x label.}

\item{ylab}{If supplied, then this is used as the y label.}

\item{main}{If supplied, then this is used as title (or z axis label).}

\item{ylim}{If supplied, then this pair of numbers are used as the y limits.}

\item{xlim}{If supplied, then this pair of numbers are used as the x limits.}

\item{too.far}{When plotting 2-D smooths, if greater than 0, then it is used to determine when a location is too
far from data to be plotted. The data are scaled into the unit square before deciding what to exclude, and \code{too.far} is a distance 
within the unit square.}

\item{trans}{Function to apply to each smooth before plotting.}

\item{seWithMean}{If \code{TRUE}, then the component smooth is shown with confidence 
intervals that include the uncertainty about the overall mean. If \code{FALSE}, then the 
uncertainty relates purely to the centred smooth itself.}

\item{...}{Other graphics parameters to pass on to plotting commands.}

}


\details{This function produces plot showing the smooth terms of a fitted semiparametric bivariate probit model.

For plots of 1-d smooths, the x axis of each plot is labelled using the name of the regressor, while the y axis is labelled as \code{s(regr,edf)} 
where \code{regr} is the regressor name, and \code{edf} the estimated degrees of freedom of the smooth.

As for 2-d smooths, perspective plots are produced with the x-axes labelled with the first and second variable names and the y axis 
is labelled as \code{s(var1,var2,edf)}, which indicates the variables of which the term is a function and the \code{edf} for the term. 

Smooths of more than 2 variables are not dealt with, but simply generate a warning.

If \code{seWithMean=TRUE}, then the confidence intervals include the uncertainty about the overall mean. That is, 
although each smooth is shown centred, the confidence intervals are obtained as if every other term in the model was 
constrained to have average 0 (average taken over the covariate values) except for the smooth being plotted. The theoretical arguments 
and simulation study of Marra and Wood (submitted paper) suggests that \code{TRUE} results in intervals with
close to nominal (frequentist) coverage probabilities.
}



\value{
The function generates plots.
}



\references{

Marra G. and Wood S. N., Coverage Properties of Confidence Intervals for Generalized Additive Model Components. Submitted.

}


\author{ 

Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

This function is pretty much the same as \code{plot.gam} in \code{mgcv}.

}



\section{WARNING}{ 

The function can not deal with smooths of more than 2 variables.

} 



\seealso{
\code{\link{ATE}}, \code{\link{InfCr}}, \code{\link{SemiParBIVProbit}}, \code{\link{summary.SemiParBIVProbit}}
}



\examples{
## see examples for SemiParBIVProbit
}

\keyword{smooth} 
\keyword{regression} 
\keyword{hplot}






