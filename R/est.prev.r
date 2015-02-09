est.prev <- function(x, sig.lev = 0.05, sw = NULL, naive = FALSE){

if(x$Model=="B" || x$Model=="BPO") stop("This function is suitable for sample selection models only.")

X2sg <- x$X2s

if(naive==TRUE)  eta2 <- X2sg%*%coef(x$gam2) 
if(naive==FALSE) eta2 <- x$eta2

if(is.null(sw)) sw <- rep(1,length(eta2)) 

core <- apply( c(dnorm(eta2))*X2sg, 2, weighted.mean,  w = sw)


if( is.null(x$X3) ) zerod <- 0
if(!is.null(x$X3) ) zerod <- rep(0, x$X3.d2)

if(naive==FALSE) G <- c( rep(0,x$X1.d2), core, zerod) 
if(naive==TRUE)  G <- c( core )  

  
  wm <- weighted.mean(pnorm(eta2), w=sw)
  
  if(naive==FALSE) Vv <- x$Vb 
  if(naive==TRUE)  Vv <- x$gam2$Vp  

  sv <- sqrt( t(G)%*%Vv%*%G ) 

  qz <- qnorm(sig.lev/2, lower.tail = FALSE)
  lb <- wm - qz*sv 
  ub <- wm + qz*sv 

  res <- c(lb,wm,ub)

  out <- list(res=res, sig.lev=sig.lev)
 
  class(out) <- "est.prev"

  out

}



