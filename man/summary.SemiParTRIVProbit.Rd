\name{summary.SemiParTRIVProbit}
\alias{summary.SemiParTRIVProbit}
\alias{print.summary.SemiParTRIVProbit}

\title{SemiParTRIVProbit summary}

\description{
It takes a fitted \code{SemiParTRIVProbit} object produced by \code{SemiParTRIVProbit()} and produces some summaries from it.

}


\usage{

\method{summary}{SemiParTRIVProbit}(object, n.sim = 100, prob.lev = 0.05, ...)

\method{print}{summary.SemiParTRIVProbit}(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)                      
}


\arguments{ 

\item{object}{A fitted \code{SemiParTRIVProbit} object as produced by \code{SemiParTRIVProbit()}.}

\item{x}{\code{summary.SemiParTRIVProbit} object produced by \code{summary.SemiParTRIVProbit()}.}


\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals for the association parameter and other measures. It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 
            
\item{digits}{Number of digits printed in output.} 

\item{signif.stars}{By default significance stars are printed alongside output.} 

\item{...}{Other arguments.}            
            

}


\details{ 

This function is very similar to \code{summary.SemiParBIVProbit()}.

\code{print.summary.SemiParTRIVProbit} prints model term summaries. 


}



\value{

Very similar to what obtained when using \code{summary.SemiParBIVProbit()}.


}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{plot.SemiParBIVProbit}}, \code{\link{predict.SemiParBIVProbit}}
}

\examples{
## see examples for SemiParTRIVProbit
}


\keyword{information criteria}
\keyword{smooth}
\keyword{regression}
\keyword{shrinkage smoother}
\keyword{variable selection}






