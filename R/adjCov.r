adjCov <- function(x, id){

cont2par <- c("N","GU","rGU","LO","LN","WEI","iG","GA","iGA")  

Vb     <- x$Vb 

if( x$margins[2]=="probit" && x$Model != "BPO0"){

if(is.null(x$X3) )  mul <- 1
if(!is.null(x$X3) ) mul <- x$X3[x$good,]

scores <- cbind( c(x$fit$dl.dbe1)*x$X1[x$good,], c(x$fit$dl.dbe2)*x$X2[x$good,], c(x$fit$dl.drho)*mul )

}

if( x$Model == "BPO0" ){

scores <- cbind( c(x$fit$dl.dbe1)*x$X1[x$good,], c(x$fit$dl.dbe2)*x$X2[x$good,] )

}



if( x$margins[2] %in% cont2par ){

if( !is.null(x$X3) && !is.null(x$X4) ) mul1 <- x$X3; mul2 <- x$X4 
if(  is.null(x$X3) &&  is.null(x$X4) ) mul1 <- mul2 <- 1 
                                       
scores <- cbind( c(x$fit$dl.dbe1)*x$X1, 
                 c(x$fit$dl.dbe2)*x$X2, 
                 c(x$fit$dl.dsigma.st)*mul1,
                 c(x$fit$dl.dteta.st)*mul2       )                                           
}





scores <- aggregate.data.frame(scores,by=list(id),FUN=sum)[,-1]
nclusters <- dim(scores)[1]
meat   <- (nclusters-1)*var(scores)
covsan <- Vb %*% meat %*% Vb
x$Vb <- covsan

rm(scores, nclusters, meat, covsan, Vb)

x

}


#\item{qu.mag, gp1, gp2, VC, Model, ig, method, magpp}{These are used for internal calculations.}
#\item{BivD, nu, sel, fp, nC}{These are used for internal calculations.}
#\item{PL, eqPL, valPL, fitPL, spPL}{These are related to the asymmetric links and used for internal calculations.}
#\item{extra.regI}{If "t" then regularization as from \code{trust} is applied to the information matrix if needed. 
#                  If different from "t" then extra regularization is applied via the options "pC" (pivoted Choleski) and
#                  "sED" (symmetric eigen-decomposition).}  
#\item{BivD, nu, PL, sel}{These are used for internal calculations.}
#\item{IV}{IV = NULL  An instrumental variable can be used if available. This can be either a binary or categorical variable. The use of a 
#
#          continuous IV is not supported.} 
# importFrom(CDVine,BiCopCDF,BiCopHfunc,BiCopPDF,BiCopMetaContour)