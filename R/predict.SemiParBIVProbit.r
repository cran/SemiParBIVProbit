predict.SemiParBIVProbit <- function(object, eq, ...){

if(object$RE==FALSE) Kk1 <- k.n <- 0
if(object$RE==TRUE && object$RE.type=="NP"){ Kk1 <- object$K - 1; k.n <- 0}
if(object$RE==TRUE && object$RE.type=="N"){  Kk1 <- 0; k.n <- 1}

 if(eq==1){ ss.pred <- object$gam1; ind <- (1:length(ss.pred$coefficients))+Kk1 } 
      else{ ss.pred <- object$gam2; fir <- length(object$gam1$coefficients) + 1 + 2*Kk1 + k.n;
                                    sec <- fir - 1 + length(ss.pred$coefficients)
                                    ind <- fir:sec }
                                   
           ss.pred$coefficients <- object$coefficients[ind]
           ss.pred$Vp <- object$Vb[ind,ind] # recall that with GEV this function is not good when predicting on the response scale

  predict.gam(ss.pred, ...)

}
