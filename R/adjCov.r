adjCov <- function(object,id){

# this does not work with power links and two-parameter copulas

Vb     <- object$Vb 
scores <- cbind( c(object$fit$dl.dbe1)*object$X1, c(object$fit$dl.dbe2)*object$X2, object$fit$dl.drho )
scores <- aggregate.data.frame(scores,by=list(id),FUN=sum)[,-1]
nclusters <- dim(scores)[1]
meat   <- (nclusters-1)*var(scores)
covsan <- Vb %*% meat %*% Vb
object$Vb <- covsan

object

}