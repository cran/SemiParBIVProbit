plot.SemiParBIVProbit <- function(x, eq, ...){

 if(eq==1){ ss.plot <- x$gam1
            ind <- (1:length(ss.plot$coefficients)) 
            } 
 if(eq==2){ ss.plot <- x$gam2
            fir <- length(x$gam1$coefficients) + 1 
            sec <- fir - 1 + length(ss.plot$coefficients)
            ind <- fir:sec }
                                   
           ss.plot$coefficients <- x$coefficients[ind]
           ss.plot$Vp <- x$Vb[ind,ind]
           ss.plot$edf <- diag(x$F)[ind]

  plot.gam(ss.plot, ...)  # recall that with PL this function is not good when plotting on the response scale
       
}

