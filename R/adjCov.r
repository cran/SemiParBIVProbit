adjCov <- function(object, id){

if(object$PL!="P") stop("This correction does not currently work for models with asymmetric links.")

Vb     <- object$Vb 
scores <- cbind( c(object$fit$dl.dbe1)*object$X1, c(object$fit$dl.dbe2)*object$X2, object$fit$dl.drho )
scores <- aggregate.data.frame(scores,by=list(id),FUN=sum)[,-1]
nclusters <- dim(scores)[1]
meat   <- (nclusters-1)*var(scores)
covsan <- Vb %*% meat %*% Vb
object$Vb <- covsan

object

}