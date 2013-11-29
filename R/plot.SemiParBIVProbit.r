plot.SemiParBIVProbit <- function(x, eq, ...){

if(x$RE==FALSE) Kk1 <- k.n <- 0
if(x$RE==TRUE && x$RE.type=="NP"){ Kk1 <- x$K - 1; k.n <- 0}
if(x$RE==TRUE && x$RE.type=="N"){  Kk1 <- 0; k.n <- 1}

 if(eq==1){ ss.plot <- x$gam1
            ind <- (1:length(ss.plot$coefficients))+Kk1 
            } 
 if(eq==2){ ss.plot <- x$gam2
            fir <- length(x$gam1$coefficients) + 1 + 2*Kk1 + k.n
            sec <- fir - 1 + length(ss.plot$coefficients)
            ind <- fir:sec }
                                   
           ss.plot$coefficients <- x$coefficients[ind]
           ss.plot$Vp <- x$Vb[ind,ind]
           ss.plot$edf <- diag(x$F)[ind]

  plot.gam(ss.plot, ...)  # recall that with PL this function is not good when plotting on the response scale
       
}

