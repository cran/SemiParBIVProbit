predict.SemiParBIVProbit <- function(object, eq, ...){

                         
 if(eq==1){ ss.pred <- object$gam1
            ind <- 1:object$X1.d2 
            } 
 if(eq==2){ ss.pred <- object$gam2
            ind <- (object$X1.d2+1):(object$X1.d2+object$X2.d2) 
            }
 if(eq==3){ ss.pred <- object$gam3
            ind <- (object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2) }                                    
                                    
                                   
           ss.pred$coefficients <- object$coefficients[ind]
           ss.pred$Vp <- object$Vb[ind,ind] # with asymmetric links this function is not good when predicting on the response scale

  predict.gam(ss.pred, ...)

}
