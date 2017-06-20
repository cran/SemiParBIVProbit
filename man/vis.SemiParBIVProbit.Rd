\name{vis.SemiParBIVProbit}
\alias{vis.SemiParBIVProbit}

\title{Visualization function}


\description{

It takes a fitted \code{SemiParBIVProbit}/\code{copulaReg}/\code{copulaSampleSel}/\code{SemiParTRIV} object produced 
by \code{SemiParBIVProbit()}, \code{copulaReg()}, \code{copulaSampleSel()}, \code{SemiParTRIV()} and 
produces perspective or contour plot views of model predictions. This function is a 
wrapper of \code{vis.gam()} in \code{mgcv}. Please see 
the documentation of \code{vis.gam()} for full details. 

}

\usage{
vis.SemiParBIVProbit(x, eq, ...)
}


\arguments{ 
  
\item{x}{A fitted \code{SemiParBIVProbit}/\code{copulaReg}/\code{copulaSampleSel}/\code{SemiParTRIV()} object.}

\item{eq}{The equation from which smooth terms should be considered for printing.}

\item{...}{Other graphics parameters to pass on to plotting commands, as described for \code{vis.gam()} in \code{mgcv}.}


}


\value{
The function generates plots.
}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{SemiParBIVProbit}}, \code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{SemiParTRIV}}, \code{\link{plot.SemiParBIVProbit}}
}




\keyword{smooth} 
\keyword{regression} 
\keyword{hplot}






