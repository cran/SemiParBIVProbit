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


#\item{qu.mag, gp1, gp2, VC, Model, ig, method, magpp}{These are used for internal calculations.}
#\item{BivD, nu, sel, fp, nC}{These are used for internal calculations.}
#\item{PL, eqPL, valPL, fitPL, spPL}{These are related to the asymmetric links and used for internal calculations.}
#\item{extra.regI}{If "t" then regularization as from \code{trust} is applied to the information matrix if needed. 
#                  If different from "t" then extra regularization is applied via the options "pC" (pivoted Choleski) and
#                  "sED" (symmetric eigen-decomposition).}  
#\item{BivD, nu, PL, sel}{These are used for internal calculations.}
#\item{IV}{IV = NULL  An instrumental variable can be used if available. This can be either a binary or categorical variable. The use of a 
#          continuous IV is not supported.} 