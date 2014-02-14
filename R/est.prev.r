est.prev <- function(x, sig.lev=0.05, sw=NULL, naive=FALSE){

if(x$sel!=TRUE) stop("This function is suitable for sample selection models only.")

good <- x$good
eta2 <- x$eta2[good]
X2sg <- x$X2s[good,]


if(naive==TRUE) eta2 <- X2sg%*%coef(x$gam2) 

if(is.null(sw)) sw <- rep(1,length(eta2)) 
core <- apply( c(dnorm(eta2))*X2sg, 2, weighted.mean,  w = sw)

if(naive==FALSE) G <- c( rep(0,x$X1.d2), core ,0) else G <- c( core )  

  
  wm <- weighted.mean(pnorm(eta2), w=sw)
  
  if(naive==FALSE) Vv <- x$Vb else Vv <- x$gam2$Vp  

  sv <- sqrt( t(G)%*%Vv%*%G ) 

  qz <- qnorm(sig.lev/2, lower.tail = FALSE)
  lb <- wm - qz*sv # lower bound
  ub <- wm + qz*sv # upper bound 

  res <- c(lb,wm,ub)

  out <- list(res=res, sig.lev=sig.lev)
 
  class(out) <- "est.prev"

  out

}



