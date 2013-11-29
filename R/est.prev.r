est.prev <- function(x, sig.lev=0.05){

if(x$sel!=TRUE) stop("This function is suitable for sample selection models only.")

good <- x$good

pr.prS    <- NA
eta1Sim   <- x$eta1S[good,]
eta2Sim   <- x$eta2S[good,] 
s.thetSim <- x$ass.pS
y1 <- x$y1[good]
y2 <- x$y2[good]

n.sel <- table(y1)[2]
n     <- sum(table(y1))

ct <- data.frame( c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","T"),
                     c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,2) 
                   )
nC <- ct[which(ct[,1]==x$BivD),2]

epsilon <- .Machine$double.eps*10^6
dimeta1Sim <- dim(eta1Sim)[2]
pr.prS  <- matrix(NA,sum(as.numeric(y1==0)),dimeta1Sim)


for(i in 1:dimeta1Sim){

if(x$BivD %in% c("T","N"))      {parS <- tanh(s.thetSim[i]); if(parS %in% c(-1,1)) parS <- sign(parS)*0.9999999}
if(x$BivD=="F")                  parS <-   s.thetSim[i] + epsilon

if(x$BivD %in% c("C0", "C180") ) parS <-   exp(s.thetSim[i]) + epsilon  
if(x$BivD %in% c("C90","C270") ) parS <- -(exp(s.thetSim[i]) + epsilon) 

if(x$BivD %in% c("J0", "J180") ) parS <-   1+exp(s.thetSim[i]) + epsilon  
if(x$BivD %in% c("J90","J270") ) parS <- -(1+exp(s.thetSim[i]) + epsilon) 

if(x$BivD %in% c("G0", "G180") ) parS <-   1+exp(s.thetSim[i])  
if(x$BivD %in% c("G90","G270") ) parS <- -(1+exp(s.thetSim[i])) 

  p0  <- pnorm(-eta1Sim[,i]) 
  if(x$BivD=="N") p01 <- pbinorm(-eta1Sim[,i],eta2Sim[,i],cov12=-parS) else p01 <- pnorm(eta2Sim[,i]) - BiCopCDF(pnorm(eta1Sim[,i]),pnorm(eta2Sim[,i]), nC, par=parS) 
  pr.prS[,i] <- (p01/p0)[y1==0] 

}

  pred.pr <- (x$p01/x$p0)[y1==0] 
  mean.pp <- mean(pred.pr)
  var.pp  <- sum(rowVars(pr.prS))/(n-n.sel)^2
  ob.pr   <- mean( y2[y1==1] ) 
  var.ob  <- (ob.pr*(1-ob.pr))/n.sel
  we      <- c( n.sel,(n-n.sel) )/n
 
  wm <- ob.pr*we[1] + mean.pp*we[2]
  sv <- sqrt(var.ob*we[1]^2 + var.pp*we[2]^2)

  qz <- qnorm(sig.lev/2, lower.tail = FALSE)
  lb <- wm - qz*sv # lower bound
  ub <- wm + qz*sv # upper bound 

  res <- c(lb,wm,ub)

  out <- list(res=res, sig.lev=sig.lev)
 
  class(out) <- "est.prev"

  out

}



