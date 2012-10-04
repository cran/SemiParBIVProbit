predict.SemiParBIVProbit <- function(object, eq, ...){

   if(eq==1) {ss.pred <- object$gam1
              fl <- length(object$gam1$coefficients) 
              ss.pred$coefficients <- object$fit$argument[1:fl]
              ss.pred$Vp <- object$Vb[1:fl,1:fl]
              
              }
              
   if(eq==2) {ss.pred <- object$gam2
              il <- length(object$gam1$coefficients)+1
              fl <- length(object$gam1$coefficients) + length(object$gam2$coefficients)
              ss.pred$coefficients <- object$fit$argument[il:fl]
              ss.pred$Vp <- object$Vb[il:fl,il:fl]
             
              }

  R <- predict.gam(ss.pred, ...)

R

}

