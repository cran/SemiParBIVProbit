\name{mb}
\alias{mb}

\title{Nonparametric (worst-case) Manski's bound}

\description{ 
\code{mb} can be used to calculate the worst-case Manski's bound.}



\usage{

mb(treat, outc)

}

\arguments{ 

\item{treat}{Binary treatment variable.}
 
\item{outc}{Binary outcome variable.} 


}



\value{

It returns the nonparametric lower and upper (worst-case) Manski's bounds. For comparison, it also returns the average 
treatment effect assuming random assignment (i.e., the treatment received relies on ignorable selection). 

}


\section{WARNINGS}{

NOTE: Since the average treatment effect can be at most 1 and the worst-case Manski's bound has width 1 then
0 is always included within the possibilites of this bound. 

}






\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\references{

Manski C.F. (1990), Nonparametric Bounds on Treatment Effects. \emph{American Economic Review, Papers and Proceedings}, 80, 319-323.

}


\seealso{
\code{\link{SemiParBIVProbit-package}}, \code{\link{SemiParBIVProbit}}  
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{Manski's bound} 
\keyword{Nonparametric bound} 
\keyword{Worst-case bound} 










