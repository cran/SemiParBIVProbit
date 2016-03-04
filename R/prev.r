prev <- function(x, sw = NULL, type = "bivariate", ind = NULL, delta = FALSE, n.sim = 100, prob.lev = 0.05, 
                      hd.plot = FALSE, main = "Histogram and Kernel Density of Simulated Prevalences", 
       xlab="Simulated Prevalences", ...){


if(x$Cont == "YES") stop("This function is not suitable for bivariate models with continuous margins.")

lb <- wm <- ub <- qz <- sv <- Vv <- G <- X2sg <- 1
wms <- NA

if(x$Model=="B" || x$Model=="BPO") stop("This function is suitable for sample selection models only.")

if( !( type %in% c("naive","univariate","bivariate") ) ) stop("Error in parameter type value. It should be one of: naive, univariate or bivariate.")


X2sg <- x$X2s


if(type == "univariate")                    {eta2 <- X2sg%*%coef(x$gam2); Vv <- x$gam2$Vp}  
if(type == "bivariate" || type == "naive")  {eta2 <- x$eta2;              Vv <- x$Vb}

if( !is.null(ind) ){ 

    if(is.logical(ind) == FALSE) stop("ind must be a logical variable.")
    if(length(table(ind))!=2 ) stop("ind must be a logical binary variable.")
    if( length(ind)!=length(eta2) ) stop("ind must have the same length as the number of observations used in fitting.")   

    eta2 <- eta2[ind]
    X2sg <- X2sg[ind,]
    
if(!is.null(sw)) sw <- sw[ind]     

}


if(is.null(sw)) sw <- rep(1,length(eta2)) 

if( length(sw)!=length(eta2) ) stop("sw must have the same length as the number of observations used in fitting.") 



#######

wm <- weighted.mean(probm(eta2, x$margins[2])$pr, w=sw)

#######



if(type != "naive"){




if(delta == TRUE){

core <- colWeightedMeans( c( probm(eta2, x$margins[2], only.pr = FALSE)$d.n )*X2sg, w = sw, na.rm = FALSE) 

if( is.null(x$X3) ) zerod <- 0
if(!is.null(x$X3) ) zerod <- rep(0, x$X3.d2)

if(type == "bivariate")   G <- c( rep(0,x$X1.d2), core, zerod) 
if(type == "univariate")  G <- c( core )  

 
  sv <- sqrt( t(G)%*%Vv%*%G ) 
  
  qz <- qnorm(prob.lev/2, lower.tail = FALSE)

  lb <- wm - qz*sv 
  ub <- wm + qz*sv 

}







if(delta == FALSE){

  if(type == "bivariate")   coefm <- x$coefficients       
  if(type == "univariate")  coefm <- x$gam2$coefficients    
  
   bs <- rMVN(n.sim, mean = coefm, sigma=Vv)
  
  
  if(type == "bivariate") bs <- bs[, x$X1.d2 + (1 : x$X2.d2) ]
 
  p2s <- probm( X2sg%*%t(bs) , x$margins[2])$pr 
  wms <- colWeightedMeans( p2s, w = sw, na.rm = FALSE)
  bb <- quantile(wms, probs = c(prob.lev/2,1-prob.lev/2), na.rm=TRUE )

  lb <- bb[1]
  ub <- bb[2] 
  
  if(hd.plot == TRUE){
  
  hist(wms*100, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(wms*100)$y,hist(wms*100, plot = FALSE)$density)), ...)
  lines(density(wms*100))

                     }

}






} # end type







if( type == "naive"){

inde <- x$inde
resp <- x$y2

if( !is.null(ind) ){ 

inde <- inde[ind]
resp <- resp[ind]
                    }

resp <- resp[inde]
sw <- sw[inde]

qz <- qnorm(prob.lev/2, lower.tail = FALSE)
wm <- weighted.mean(resp, w = sw)
sv <- sqrt( (wm*(1 - wm))/length(resp) )
lb <- wm - qz*sv 
ub <- wm + qz*sv 
  
}
  










  res <- c(lb, wm, ub)
  

  rm(lb,wm,ub,qz,sv,Vv,G,X2sg)

  out <- list(res=res, prob.lev=prob.lev, sim.prev = wms)
 
  class(out) <- "prev"

  out

}

























