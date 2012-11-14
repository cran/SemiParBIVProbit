predict.SemiParBIVProbit <- function(object, eq, ...){

 if(eq==1){ss.pred <- object$gam1; ind <- 1:length(ss.pred$coefficients)} 
      else{ss.pred <- object$gam2; ind <- ( length(object$gam1$coefficients)+1 ):( length(object$gam1$coefficients) + length(ss.pred$coefficients) ) }

           ss.pred$coefficients <- object$fit$argument[ind]
           ss.pred$Vp <- object$Vb[ind,ind]

  R <- predict.gam(ss.pred, ...)

R

}

