plot.SemiParBIVProbit <- function(x, eq, ...){

 if(eq==1){ ss.plot <- x$gam1
            ind <- 1:x$X1.d2 
            } 
 if(eq==2){ ss.plot <- x$gam2
            ind <- (x$X1.d2+1):(x$X1.d2+x$X2.d2) 
            }
 if(eq==3){ ss.plot <- x$gam3
            ind <- (x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2) }
                                    
           ss.plot$coefficients <- x$coefficients[ind]
           ss.plot$Vp <- x$Vb[ind,ind]
           ss.plot$edf <- diag(x$F)[ind]

  plot.gam(ss.plot, ...)  # with asymmetric links this function is not good when plotting on the response scale
       
}

