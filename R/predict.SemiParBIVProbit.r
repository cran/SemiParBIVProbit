predict.SemiParBIVProbit <- function(object, eq, ...){

 if(eq==1){ ss.pred <- object$gam1; ind <- (1:length(ss.pred$coefficients)) } 
      else{ ss.pred <- object$gam2; fir <- length(object$gam1$coefficients) + 1
                                    sec <- fir - 1 + length(ss.pred$coefficients)
                                    ind <- fir:sec }
                                   
           ss.pred$coefficients <- object$coefficients[ind]
           ss.pred$Vp <- object$Vb[ind,ind] # recall that with PL this function is not good when predicting on the response scale

  predict.gam(ss.pred, ...)

}
