\name{summary.copulaSampleSel}
\alias{summary.copulaSampleSel}
\alias{print.summary.copulaSampleSel}

\title{copulaSampleSel summary}

\description{
It takes a fitted \code{copulaSampleSel} object produced by \code{copulaSampleSel()} and produces some summaries from it.

}


\usage{

\method{summary}{copulaSampleSel}(object, n.sim = 100, prob.lev = 0.05, cm.plot = FALSE, 
        xlim = c(-3, 3), ylim = c(-3, 3), 
        ylab = "Margin 2", xlab = "Margin 1", 
        n.grid = 1000, n.dig = 2, ...)
  
\method{print}{summary.copulaSampleSel}(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)  
  
}


\arguments{ 

\item{object}{A fitted \code{copulaSampleSel} object as produced by \code{copulaSampleSel()}.}

\item{x}{\code{summary.copulaSampleSel} object produced by \code{summary.copulaSampleSel()}.}


\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals for the association parameter, dispersion coefficient, for instance It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 
                    
\item{cm.plot}{If \code{TRUE} then a filled bivariate contour meta plot corresponding to the assumed (estimated) bivariate model is produced.}   
             
\item{xlim, ylim}{Limits of the bivariate contour meta plot.}    

\item{ylab, xlab}{Labels for the bivariate contour meta plot.} 
                     
\item{n.grid}{Number of grid points used in contour plot construction. This is relevant for the continuous margin.} 

\item{n.dig}{Number of digit points used in rounding bivariate pdf values for contour plot construction. This is relevant for the continuous margin.} 
          
\item{digits}{Number of digits printed in output.} 

\item{signif.stars}{By default significance stars are printed alongside output.}           
          
\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{cm.plot=TRUE}.}

}


\details{ 

This function is very similar to \code{summary.SemiParBIVProbit()}.

\code{print.summary.copulaSampleSel} prints model term summaries. 


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
## see examples for copulaSampleSel
}


\keyword{information criteria}
\keyword{smooth}
\keyword{regression}
\keyword{non-random sample selection}






