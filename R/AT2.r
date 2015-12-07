AT2 <- function(x1, x2, index1, index2, n.sim = 100, prob.lev = 0.05, 
                  hd.plot = FALSE, 
                  main = "Histogram and Kernel Density of Simulated Average Effects", 
                  xlab = "Simulated Average Effects", ...){


etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- 1
diffEf <- fy1.y2 <- est.ATso <- y2 <- CIF <- Pr <- Effects <- C.11 <- C.10 <- C.11s <- C.10s <- v0s <- v1s <- ATs <- NULL


if( missing(index1) || missing(index2) ) stop("You must provide values for index1 and index2.")  

m2 <- c("N","GU","rGU","LO","LN","WEI","WEI2","iG","GA")
m3 <- c("DAGUM")

end <- 0
epsilon <- 0.0000001 # 0.9999999 sqrt(.Machine$double.eps)
max.p   <- 0.9999999
est.ATb <- NA
indD <- list()


if(x2$margins[2] == "DAGUM") { if( min(sqrt(x2$sigma2[index2])) <= 1) stop("Sigma parameter has value(s) smaller than 1, hence the mean is indeterminate.")}
 

########################################################
# model x1 - binary binary
########################################################

ind1  <- 1:x1$X1.d2 
ind2  <- x1$X1.d2+(1:x1$X2.d2)

X1 <- as.matrix(x1$X1[index1,])
X2 <- as.matrix(x1$X2[index1,])
eta1 <- x1$eta1[index1]
eta2 <- x1$eta2[index1] 

p.1 <- pnorm(eta1)
p.2 <- pnorm(eta2)

if( is.null(x1$X3)  ) ass.p <- x1$theta      
if( !is.null(x1$X3) ) ass.p <- x1$theta[index1]  


p.1 <-   pmax(p.1, epsilon ) 
p.1 <- ifelse(p.1 > max.p, max.p, p.1)
p.2 <-   pmax(p.2, epsilon )
p.2 <- ifelse(p.2 > max.p, max.p, p.2)

C.11  <- (pmax( BiCDF(p.1, p.2, x1$nC, ass.p) , epsilon )) / p.1
C.10  <- (p.2 - C.11) / (1 - p.1)                                  ## second equation is of interest
              
              
######
# CIs
######
              
bs <- rMVN(n.sim, mean = coef(x1), sigma=x1$Vb)

eta1s <- t(X1)%*%t(bs[,ind1])
eta2s <- t(X2)%*%t(bs[,ind2])  

p.1s <-   pmax(pnorm(eta1s), epsilon ) 
p.1s <- ifelse(p.1s > max.p, max.p, p.1s)
p.2s <-   pmax(pnorm(eta2s), epsilon )
p.2s <- ifelse(p.2s > max.p, max.p, p.2s)

if( !is.null(x1$X3) ) etds <- t(x1$X3[index1,])%*%t(bs[,(x1$X1.d2+x1$X2.d2+1):(x1$X1.d2+x1$X2.d2+x1$X3.d2)])
if(  is.null(x1$X3) ) etds <- bs[, length(coef(x1))]
   

if(x1$BivD=="F")                   ass.ps <- etds + epsilon
if(x1$BivD %in% c("N"))           {ass.ps <- tanh(etds); ass.ps <- ifelse(ass.ps < -max.p, -max.p, ass.ps) 
                                                         ass.ps <- ifelse(ass.ps >  max.p,  max.p, ass.ps) }

if(x1$BivD %in% c("C0", "C180") )  ass.ps <-   exp(etds) + epsilon
if(x1$BivD %in% c("C90","C270") )  ass.ps <- -(exp(etds) + epsilon)

if(x1$BivD %in% c("J0", "J180") )  ass.ps <-    1 + exp(etds) + epsilon
if(x1$BivD %in% c("J90","J270") )  ass.ps <- -( 1 + exp(etds) + epsilon)

if(x1$BivD %in% c("G0", "G180") )  ass.ps <-    1 + exp(etds)
if(x1$BivD %in% c("G90","G270") )  ass.ps <- -( 1 + exp(etds) )

ass.ps <- ifelse(ass.ps == Inf ,  8.218407e+307, ass.ps )
ass.ps <- ifelse(ass.ps == -Inf, -8.218407e+307, ass.ps )


if(x1$BivD == "N") {

	for(i in 1:n.sim){ 
		 C.11s[i] <- pmax( BiCDF(p.1s[i], p.2s[i], x1$nC, ass.ps[i], test = FALSE), epsilon )  / p.1s[i]
		 C.10s[i] <- (p.2s[i] - C.11s[i]) / (1 - p.1s[i])  
	                  }                 
}

if(x1$BivD != "N") {

 C.11s <- pmax( BiCDF(p.1s, p.2s, x1$nC, ass.ps, test = FALSE), epsilon ) / p.1s 
 C.10s <- (p.2s - C.11s) / (1 - p.1s)
             
}





########################################################
# model x2
########################################################


if( x2$margins[2] %in% c("N","GU","rGU","LO") )                     { lil <- -Inf;    uil <- Inf}
if( x2$margins[2] %in% c("LN","WEI","WEI2","iG","GA","DAGUM")  )    { lil <- epsilon; uil <- Inf}


ConExp0 <- function(y2, eta2, sigma2, nu, margin2, p.1, ass.p, BivD){   

	pp0 <- distrHsAT(y2, eta2, sigma2, nu, margin2) 
	p2.0   <- pp0$p2
	pdf2.0 <- pp0$pdf2 
	dc0 <- copgHsAT(p1 = p.1, p2=p2.0, teta=ass.p, BivD=BivD)$c.copula.be2	
	cond0 <- y2*dc0*pdf2.0/p.1
	cond0
                                                                   }

ConExp1 <- function(y2, eta2, sigma2, nu, margin2, p.1, ass.p, BivD){   

	pp1 <- distrHsAT(y2, eta2, sigma2, nu, margin2)
	p2.1   <- pp1$p2
	pdf2.1 <- pp1$pdf2 	
	dc1 <- copgHsAT(p1 = p.1, p2=p2.1, teta=ass.p, BivD=BivD)$c.copula.be2	
	cond1 <- y2 * ( (1 - dc1)*pdf2.1/(1-p.1))
	cond1
                                                                   }
                                                                   

integr0 <- function(eta2, sigma2, nu, margin2, p.1, ass.p, BivD, lil, uil){

  integrate(ConExp0, lower=lil, upper=uil, eta2=eta2, 
            sigma2=sigma2, nu = nu, margin2=margin2, p.1=p.1, 
            ass.p=ass.p, BivD=BivD)$value
                         
                       }

integr1 <- function(eta2, sigma2, nu, margin2, p.1, ass.p, BivD, lil, uil){

  integrate(ConExp1, lower=lil, upper=uil, eta2=eta2,  
            sigma2=sigma2, nu= nu, margin2=margin2, p.1 = p.1, 
            ass.p=ass.p, BivD=BivD)$value
                         
                       }

v.integr0 <- Vectorize(integr0)  
v.integr1 <- Vectorize(integr1) 


#######################################################


ind1  <- 1:x2$X1.d2 
ind2  <- x2$X1.d2+(1:x2$X2.d2)

X1 <- as.matrix(x2$X1[index2,])
X2 <- as.matrix(x2$X2[index2,])


eta1 <- x2$eta1[index2]
eta2 <- x2$eta2[index2] 


ass.p  <- x2$theta
sigma2 <- x2$sigma2
nu     <- x2$nu



if( x2$margins[2] %in% m2){

if( !is.null(x2$X3) && !is.null(x2$X4) ) { ass.p <- ass.p[index2]; sigma2 <- sigma2[index2] }

}



if( x2$margins[2] %in% m3){

if( !is.null(x2$X3) && !is.null(x2$X4) && !is.null(x2$X5)) { ass.p <- ass.p[index2]; sigma2 <- sigma2[index2]; nu <- nu[index2] }

}




p.1  <- pnorm(-x2$eta1[index2])
p.1 <- pmax(p.1, epsilon ) 
p.1 <- ifelse(p.1 > max.p, max.p, p.1)


v0 <- v.integr0(eta2 = x2$eta2[index2], sigma2=sigma2, nu = nu, margin2=x2$VC$margins[2], p.1=p.1, 
                ass.p=ass.p, BivD=x2$BivD, lil=lil, uil=uil)
v1 <- v.integr1(eta2=x2$eta2[index2], sigma2=sigma2, nu = nu, margin2=x2$VC$margins[2], p.1=p.1, 
                ass.p=ass.p, BivD=x2$BivD, lil=lil, uil=uil)








####################
# Effect of interest
####################

AT <- C.11*v1 - C.10*v0  
   

#######################



bs <- rMVN(n.sim, mean = coef(x2), sigma=x2$Vb)

eta1s <- t(X1)%*%t(bs[,ind1])
eta2s <- t(X2)%*%t(bs[,ind2])  

p.1s <-   pmax(pnorm(-eta1s), epsilon ) 
p.1s <- ifelse(p.1s > max.p, max.p, p.1s)




   
   if(x2$VC$margins[2] %in% m2 ){
   
   if( !is.null(x2$X3) ) sigma2.st <- x2$X3[index2,]%*%t(bs[,(x2$X1.d2+x2$X2.d2+1):(x2$X1.d2+x2$X2.d2+x2$X3.d2)]) 
   if(  is.null(x2$X3) ) sigma2.st <- bs[, length(coef(x2)) - 1]
   
   if( !is.null(x2$X4) ) etds <- x2$X4[index2,]%*%t(bs[,(x2$X1.d2+x2$X2.d2+x2$X3.d2 + 1):(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2)])
   if(  is.null(x2$X4) ) etds <- bs[, length(coef(x2))]  
   
   sigma2.st <- ifelse( sigma2.st > 20, 20, sigma2.st )  
   sigma2.st <- ifelse( sigma2.st < -17, -17, sigma2.st ) 
   sigma2s <- exp(sigma2.st)
  
  }
  
   if(x2$VC$margins[2] %in% m3 ){
   
   if( !is.null(x2$X3) ) sigma2.st <- x2$X3[index2,]%*%t(bs[,(x2$X1.d2+x2$X2.d2+1):(x2$X1.d2+x2$X2.d2+x2$X3.d2)]) 
   if(  is.null(x2$X3) ) sigma2.st <- bs[,length(coef(x2)) - 2]
   
   if( !is.null(x2$X4) ) nu.st <- x2$X4[index2,]%*%t(bs[,(x2$X1.d2+x2$X2.d2+x2$X3.d2+1):(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2)]) 
   if(  is.null(x2$X4) ) nu.st <- bs[, length(coef(x2)) - 1]   
   
   if( !is.null(x2$X5) ) etds <- x2$X5[index2, ]%*%t(bs[,(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2 + 1):(x2$X1.d2+x2$X2.d2+x2$X3.d2+x2$X4.d2+x2$X5.d2)])
   if(  is.null(x2$X5) ) etds <- bs[,length(coef(x2))]  
   
   sigma2.st <- ifelse( sigma2.st > 20, 20, sigma2.st )  
   sigma2.st <- ifelse( sigma2.st < -17, -17, sigma2.st ) 
   sigma2s <- exp(sigma2.st)
   
   
   #if(x2$VC$margins[2] == "DAGUM"){
   nu.st <- ifelse( nu.st > 20, 20, nu.st )  
   nu.st <- ifelse( nu.st < -17, -17, nu.st ) 
   nus <- exp(nu.st)   
   #}
   

  
  }  
 
 
 
 
    if(x2$BivD=="F")                   ass.ps <- etds + epsilon
    if(x2$BivD %in% c("N"))           {ass.ps <- tanh(etds); ass.ps <- ifelse(ass.ps < -max.p, -max.p, ass.ps) 
                                                            ass.ps <- ifelse(ass.ps >  max.p,  max.p, ass.ps) }
   
    if(x2$BivD %in% c("C0", "C180") )  ass.ps <-   exp(etds) + epsilon
    if(x2$BivD %in% c("C90","C270") )  ass.ps <- -(exp(etds) + epsilon)
 
    if(x2$BivD %in% c("J0", "J180") )  ass.ps <-    1 + exp(etds) + epsilon
    if(x2$BivD %in% c("J90","J270") )  ass.ps <- -( 1 + exp(etds) + epsilon)
  
    if(x2$BivD %in% c("G0", "G180") )  ass.ps <-    1 + exp(etds)
    if(x2$BivD %in% c("G90","G270") )  ass.ps <- -( 1 + exp(etds) )
    
    ass.ps <- ifelse(ass.ps == Inf ,  8.218407e+307, ass.ps )
    ass.ps <- ifelse(ass.ps == -Inf, -8.218407e+307, ass.ps )



for(i in 1:n.sim){ 

v0s[i] <- v.integr0(eta2 = eta2s[i], sigma2=sigma2s[i], nu=nus[i], margin2=x2$VC$margins[2], p.1 = p.1s[i], ass.p=ass.ps[i], BivD=x2$BivD, lil=lil, uil=uil)
v1s[i] <- v.integr1(eta2 = eta2s[i], sigma2=sigma2s[i], nu=nus[i], margin2=x2$VC$margins[2], p.1 = p.1s[i], ass.p=ass.ps[i], BivD=x2$BivD, lil=lil, uil=uil)


                 }




 # if(hd.plot == TRUE){
 # 
 # if(x$margins[2] %in% m2) mult <- 1 else mult <- 100
 # 
 # hist(est.ATb*mult, freq = FALSE, main=main, 
 #      xlab=xlab, 
 #      ylim=c(0,max(density(est.ATb*mult)$y,hist(est.ATb*mult, plot = FALSE)$density)), ...)
 # lines(density(est.ATb*mult))
 #
 #}
#res <- c(CIs[1], est.AT, CIs[2])
#
#out <- list(res=res, prob.lev=prob.lev, sim.AT=est.ATb, AT.so = est.ATso, mar2=x$margins[2], type = type, 
#            Effects = Effects, Pr = Pr, treat = y2, eq = eq, p11 = C.11, p10 = C.10)
# 							 
#rm(etap.noi, X.int, X.noi, eti1, eti0, etno, indS, bs, ind.excl, p.int1, p.int0, d.int1, d.int0,
#   p.etn, d.etn, ass.p, ass.pst, C.11, C.10, sig2, peti1s, peti0s, sigma2.st, sigma2s, eti1s, eti0s, d0, d1,
#   p.etns, etnos, etds, ass.ps)    
# 
#class(out) <- "AT"
#
#out


for(i in 1:n.sim) ATs[i] <- C.11s[i]*v1s[i] - C.10s[i]*v0s[i]



 if(hd.plot == TRUE){
  
  hist(ATs, freq = FALSE, main = main, xlab = xlab, 
       ylim = c(0, max(density(ATs)$y, hist(ATs, plot = FALSE)$density)), ...)
  lines(density(ATs))
 
 }




CIs <- as.numeric(quantile(ATs, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))

res <- c(CIs[1], AT, CIs[2])

out <- list(res = res, prob.lev = prob.lev) 

class(out) <- "AT2"

out

}

