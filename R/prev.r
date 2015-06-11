prev <- function(x, sw = NULL, naive = FALSE, ind = NULL, delta = FALSE, n.sim = 100, prob.lev = 0.05, 
                      hd.plot = FALSE, main = "Histogram and Kernel Density of Simulated Prevalences", 
       xlab="Simulated Prevalences", ...){

lb <- wm <- ub <- qz <- sv <- Vv <- G <- X2sg <- 1
wms <- NA

if(x$Model=="B" || x$Model=="BPO") stop("This function is suitable for sample selection models only.")


X2sg <- x$X2s



if(naive==TRUE)  {eta2 <- X2sg%*%coef(x$gam2); Vv <- x$gam2$Vp}  
if(naive==FALSE) {eta2 <- x$eta2;              Vv <- x$Vb}


if( !is.null(ind) ){ 

    if(is.logical(ind) == FALSE) stop("ind must be a logical variable.")
    if(length(table(ind))!=2 ) stop("ind must be a logical binary variable.")
    if( length(ind)!=length(eta2) ) stop("ind must have the same length as the number of observations used in fitting.")   

    eta2 <- eta2[ind]
    X2sg <- X2sg[ind,]
    
if(!is.null(sw)) sw <- sw[ind]     

}


if(is.null(sw)) sw <- rep(1,length(eta2)) 

wm <- weighted.mean(pnorm(eta2), w=sw)



if(delta == TRUE){

core <- apply( c(dnorm(eta2))*X2sg, 2, weighted.mean,  w = sw)


if( is.null(x$X3) ) zerod <- 0
if(!is.null(x$X3) ) zerod <- rep(0, x$X3.d2)

if(naive==FALSE) G <- c( rep(0,x$X1.d2), core, zerod) 
if(naive==TRUE)  G <- c( core )  

 
  sv <- sqrt( t(G)%*%Vv%*%G ) 
  
  qz <- qnorm(prob.lev/2, lower.tail = FALSE)

  lb <- wm - qz*sv 
  ub <- wm + qz*sv 

}



if(delta == FALSE){

  if(naive==FALSE) coefm <- x$coefficients       
  if(naive==TRUE)  coefm <- x$gam2$coefficients    

  # bs <- rmvnorm(n.sim, mean = coefm, sigma=Vv, method=s.meth)
  
   bs <- rMVN(n.sim, mean = coefm, sigma=Vv)
  
  
  if(naive==FALSE) bs <- bs[, x$X1.d2 + (1 : x$X2.d2) ]
 
  p2s <- pnorm(X2sg%*%t(bs)) 
  wms <- apply(p2s, MARGIN=2, FUN=weighted.mean, w=sw)  
  bb <- quantile(wms, probs=c(prob.lev/2,1-prob.lev/2), na.rm=TRUE )

  lb <- bb[1]
  ub <- bb[2] 
  
  if(hd.plot == TRUE){
  
  hist(wms*100, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(wms*100)$y,hist(wms*100, plot = FALSE)$density)), ...)
  lines(density(wms*100))

 }

}


  res <- c(lb, wm, ub)
  

  rm(lb,wm,ub,qz,sv,Vv,G,X2sg)

  out <- list(res=res, prob.lev=prob.lev, sim.prev = wms)
 
  class(out) <- "prev"

  out

}



