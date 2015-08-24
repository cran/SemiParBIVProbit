summary.SemiParBIVProbit <- function(object, n.sim = 100, prob.lev = 0.05, 
                                     cm.plot = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), 
                                     ylab = "Margin 2", xlab = "Margin 1", gm = FALSE, ...){


  testStat <- getFromNamespace("testStat", "mgcv")
  liu2   <- getFromNamespace("liu2", "mgcv") 
  filled.contour <- getFromNamespace("filled.contour", "graphics")      

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- et1s <- et2s <- p1s <- p2s <- p11s <- p10s <- p00s <- p01s <- ORs <- GMs <- XX <- Xt <- V <- 1
  
  cont2par <- c("N","GU","rGU","LO","LN","WEI","WEI2","iG","GA","iGA")  
  cont3par <- c("DAGUM")  
  

  n <- object$n; n.sel <- object$n.sel
  
  
  
  tableN <- table <- list(NULL, NULL, NULL, NULL, NULL, NULL)
  CIkt <- CIor <- CIgm <- CIsig2 <- CInu <- NULL  
  epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)
  max.p   <- 0.9999999
  
  
  est.RHOb <- rep(NA,n.sim) 

 
  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 

  
  if(object$VC$Model != "BPO0") bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb)  


  if(object$VC$Model == "BPO0") epds <- rep(0, 10 )

  if(object$VC$margins[2]=="probit" && object$VC$Model != "BPO0"){
  
  if( !is.null(object$X3) ) epds <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)])
  if(  is.null(object$X3) ) epds <- bs[,lf]
  
  }
  
  
  
  
  
  if(object$VC$margins[2] %in% cont2par ){
  
  if( !is.null(object$X3) ) sigma2.st <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
  if(  is.null(object$X3) ) sigma2.st <- bs[,lf-1]
  
   sigma2.st <- ifelse( sigma2.st > 20, 20, sigma2.st )  
   sigma2.st <- ifelse( sigma2.st < -17, -17, sigma2.st ) 
   sigma2 <- exp(sigma2.st)
   if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
   
   CIsig2 <- t(apply(sigma2, MARGIN=1, FUN=quantile, probs=c(prob.lev/2,1-prob.lev/2), na.rm=TRUE ))  
  
  if( !is.null(object$X4) ) epds <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)])
  if(  is.null(object$X4) ) epds <- bs[,lf]  
   
  } 
  
  
  
  
    if(object$VC$margins[2] %in% cont3par ){
    
    if( !is.null(object$X3) ) sigma2.st <- object$X3%*%t(bs[,(object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)]) 
    if(  is.null(object$X3) ) sigma2.st <- bs[,lf-2]
    
    if( !is.null(object$X4) ) nu.st     <- object$X4%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)]) 
    if(  is.null(object$X4) ) nu.st     <- bs[,lf-1]    
    
     sigma2.st <- ifelse( sigma2.st > 20, 20, sigma2.st )  
     sigma2.st <- ifelse( sigma2.st < -17, -17, sigma2.st ) 
     sigma2 <- exp(sigma2.st)
     if(  is.null(object$X3) ) sigma2 <- t(as.matrix(sigma2))
     
     CIsig2 <- t(apply(sigma2, MARGIN = 1, FUN = quantile, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE ))  
     
     nu.st <- ifelse( nu.st > 20, 20, nu.st )  
     nu.st <- ifelse( nu.st < -17, -17, nu.st ) 
     nu <- exp(nu.st)
     if(  is.null(object$X4) ) nu <- t(as.matrix(nu))
     
     CInu <- t(apply(nu, MARGIN = 1, FUN = quantile, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE ))      
     
     
    
    if( !is.null(object$X5) ) epds <- object$X5%*%t(bs[,(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2 + 1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)])
    if(  is.null(object$X5) ) epds <- bs[,lf]  
     
  }
  
  
  
        
   if(object$BivD=="N")                 {est.RHOb <- tanh(epds); est.RHOb <- ifelse(est.RHOb < -max.p, -max.p, est.RHOb)
                                                                 est.RHOb <- ifelse(est.RHOb >  max.p , max.p, est.RHOb)}
   if(object$BivD=="F")                  est.RHOb <- epds + epsilon

   if(object$BivD %in% c("C0", "C180") ) est.RHOb <-   exp(epds) + epsilon  
   if(object$BivD %in% c("C90","C270") ) est.RHOb <- -(exp(epds) + epsilon) 

   if(object$BivD %in% c("J0", "J180") ) est.RHOb <-   1+exp(epds) + epsilon  
   if(object$BivD %in% c("J90","J270") ) est.RHOb <- -(1+exp(epds) + epsilon) 

   if(object$BivD %in% c("G0", "G180") ) est.RHOb <-   1+exp(epds)  
   if(object$BivD %in% c("G90","G270") ) est.RHOb <- -(1+exp(epds))
   
   est.RHOb <- ifelse(est.RHOb ==  Inf,  8.218407e+307, est.RHOb) 
   est.RHOb <- ifelse(est.RHOb == -Inf, -8.218407e+307, est.RHOb)    
   
   if(  is.null(object$X3) ) est.RHOb <- t(as.matrix(est.RHOb))
   
   
   CIrs <- t(apply(est.RHOb, MARGIN = 1, FUN = quantile, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE ))
   
   
   
   
   
########################   
   
   
if (gm == TRUE){   


   ####
   # for OR and GM
   ##   
   
if(object$VC$margins[2]=="probit" && object$VC$Model != "BPO0"){   
     
   et1s <- object$X1%*%t(bs[,1:object$X1.d2])    
   et2s <- object$X2%*%t(bs[,(object$X1.d2+1):(object$X1.d2+object$X2.d2)]) 
     
   p1s <- pnorm(et1s)
   p1s <- pmax(p1s, epsilon )
   p1s <- ifelse(p1s > max.p,max.p,p1s)  
   p2s <- pnorm(et2s)
   p2s <- pmax(p2s, epsilon )
   p2s <- ifelse(p2s > max.p,max.p,p2s) 
   
   p11s <- matrix(NA,dim(p1s)[1],dim(p1s)[2])
     
   if( !is.null(object$X3) ) { for(i in 1:n.sim) p11s[,i] <- BiCDF(p1s[,i], p2s[,i], object$nC, est.RHOb[,i]) }
   if(  is.null(object$X3) ) { for(i in 1:n.sim) p11s[,i] <- BiCDF(p1s[,i], p2s[,i], object$nC, est.RHOb[i])  }
    
 p11s <- pmax(p11s, epsilon )
 p11s <- ifelse(p11s > max.p, max.p, p11s)  
 p10s <- p1s - p11s 
 p00s <- (1 - p2s) - ( p1s - p11s )
 p01s <- p2s - p11s
 

ORs <- (p00s*p11s)/(p01s*p10s)

ORs  <- ifelse(ORs  ==  Inf,  8.218407e+307, ORs ) 
ORs  <- ifelse(ORs  == -Inf, -8.218407e+307, ORs ) 


GMs <- colMeans((ORs - 1)/(ORs + 1))
ORs <- colMeans(ORs)


  CIor <- as.numeric(quantile(ORs,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))
  CIgm <- as.numeric(quantile(GMs,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))


}


if(object$VC$gc.l == TRUE) gc()

}
 
#########################
         

  
  
  if(object$VC$gc.l == TRUE) gc()

  index <- 1:2
  ind1 <- 1:object$gp1
  ind2 <- object$X1.d2 + (1:object$gp2)
  ind3 <- ind4 <- ind5 <- ind6 <- NULL 
  
  if(!is.null(object$X3) ) {
  
       ind3 <- object$X1.d2 + object$X2.d2 + (1:object$gp3)
       index <- 1:3
       
       if(!is.null(object$X4) ) {
       ind4 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + (1:object$gp4)
       index <- 1:4
       }
                                
       if(!is.null(object$X5) ) {
       ind5 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + (1:object$gp5)
       index <- 1:5  
       }     
                                
       if(!is.null(object$X6) ) {
       ind6 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + (1:object$gp6)
       index <- 1:6  
       }                                 
                            
  }
                            
  ind <- list( ind1 = ind1,
               ind2 = ind2,
               ind3 = ind3, 
               ind4 = ind4,
               ind5 = ind5,
               ind6 = ind6)
                

  for(i in index){
  estimate <- coef(object)[ind[[i]]]
  se       <- SE[ind[[i]]]
  ratio    <- estimate/se
  pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
  table[[i]] <- cbind(estimate,se,ratio,pv)
  dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }

  
  if( object$l.sp1!=0 || object$l.sp2!=0 || object$l.sp3!=0 || object$l.sp4!=0 || object$l.sp5!=0 || object$l.sp6!=0){

  	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- tableN <- list(0, 0, 0, 0, 0, 0)
        XX <- object$R
        
           for(i in index){

             if(i==1) {mm <- object$l.sp1; if(mm==0) next}
             if(i==2) {mm <- object$l.sp2; if(mm==0) next} 
             if(i==3) {mm <- object$l.sp3; if(mm==0) next} 
             if(i==4) {mm <- object$l.sp4; if(mm==0) next} 
             if(i==5) {mm <- object$l.sp5; if(mm==0) next} 
             if(i==6) {mm <- object$l.sp6; if(mm==0) break} 
  
		for(k in 1:mm){

                        if(i==1){ gam <- object$gam1; ind <-  gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para                                } 
                        if(i==2){ gam <- object$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2                } 
                        if(i==3){ gam <- object$gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 }
                        if(i==4){ gam <- object$gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 }
                        if(i==5){ gam <- object$gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 }
                        if(i==6){ gam <- object$gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 }
                          
                        gam$sig2            <- 1
                        gam$scale.estimated <- FALSE                          
                          
                        if(gam$smooth[[k]]$null.space.dim == 0){
                        
                        LRB <- rbind(XX, t(mroot(object$fit$S.h)))
			LRB <- cbind(LRB[, -ind], LRB[, ind])
			ind1 <- (ncol(LRB) - length(ind) + 1):ncol(LRB)
			Rm <- qr.R(qr(LRB, tol = 0, LAPACK = FALSE))[ind1, ind1]
                        B <- mroot(object$Ve[ind, ind, drop = FALSE])
                          
			b.hat <- coef(object)[ind]
			d <- Rm %*% b.hat
			stat <- sum(d^2)
			ev <- eigen(crossprod(Rm %*% B), symmetric = TRUE, only.values = TRUE)$values
			ev[ev < 0] <- 0
			rank <- sum(ev > max(ev) * .Machine$double.eps^0.8)
			pval <- liu2(stat, ev)                          
                        Tp <- list(stat = stat, pval = pval, rank = rank)  
                          
                        }
                          
			if(gam$smooth[[k]]$null.space.dim != 0){
			
			b  <- coef(object)[ind]
			V  <- Vb[ind,ind, drop = FALSE]
			Xt <- XX[, ind, drop = FALSE] 
			pTerms.df[[i]][k] <- min(ncol(Xt), object$edf11[[i]][k])
			Tp <- testStat(b, Xt, V, pTerms.df[[i]][k], type = 0, res.df = -1)
			
			}
			
			
			pTerms.chi.sq[[i]][k] <- Tp$stat 
			pTerms.df[[i]][k] <- Tp$rank
                        pTerms.pv[[i]][k] <- Tp$pval
			                 
                }
                
              tableN[[i]] <- cbind(object$edf[[i]], pTerms.df[[i]], pTerms.chi.sq[[i]], pTerms.pv[[i]])
              dimnames(tableN[[i]])[[2]] <- c("edf", "Ref.df", "Chi.sq", "p-value")
              
            }

  if(object$VC$gc.l == TRUE) gc()

  }
  










 if(cm.plot == TRUE && object$VC$Model != "BPO0"){
 
 
 
 
 m2 <- sqs2 <- nu <- 0 
 
 if(object$margins[2] %in% cont2par ){ 
 m2 <- mean(object$eta2) 
 sqs2 <- sqrt(object$sigma2.a)
 nu <- 0
 } 
 
 if(object$margins[2] %in% cont3par ){ 
 m2 <- mean(object$eta2) 
 sqs2 <- sqrt(object$sigma2.a)
 nu <- object$nu.a
 } 

 par1 <- object$theta.a 


 if(object$BivD=="N")    {cop <- bquote(paste("Gaussian (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="F")    {cop <- bquote(paste("Frank (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="C0")   {cop <- bquote(paste("Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="C90")  {cop <- bquote(paste("90",degree," Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="C180") {cop <- bquote(paste("180",degree, " Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="C270") {cop <- bquote(paste("270",degree, " Clayton (",hat(theta)," = ",.(round(par1,2)),")",sep=""))} 
 if(object$BivD=="J0")   {cop <- bquote(paste("Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="J90")  {cop <- bquote(paste("90",degree," Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="J180") {cop <- bquote(paste("180",degree," Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="J270") {cop <- bquote(paste("270",degree," Joe (",hat(theta)," = ",.(round(par1,2)),")",sep=""))} 
 if(object$BivD=="G0")   {cop <- bquote(paste("Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="G90")  {cop <- bquote(paste("90",degree," Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="G180") {cop <- bquote(paste("180",degree," Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))}
 if(object$BivD=="G270") {cop <- bquote(paste("270",degree," Gumbel (",hat(theta)," = ",.(round(par1,2)),")",sep=""))} 
 
 
 funcm2 <- function(x22, m2, sqs2, nu, mar2){
           
           
           if(mar2=="probit")      {d.x2 <-  dnorm(x22, mean = 0, sd = 1)
                                    p2 <-  pnorm(x22, mean = 0, sd = 1)
           }
           if(mar2=="N")           {d.x2 <-  dnorm(x22, mean = m2, sd = sqs2)
           			   p2 <-  pnorm(x22, mean = m2, sd = sqs2)
           }
           if(mar2=="LN")          {d.x2 <- dlnorm(x22, meanlog = m2, sdlog = sqs2) 
           			   p2 <- plnorm(x22, meanlog = m2, sdlog = sqs2)
           }
           if(mar2=="GU")          {d.x2 <- exp(-exp((x22 - m2)/sqs2)) * (exp((x22 - m2)/sqs2)*(1/sqs2))   
           			   p2 <- 1 - exp(-exp((x22 - m2)/sqs2)) 
           }
           if(mar2=="rGU")         {d.x2 <- 1/sqs2 * exp( -( (x22-m2)/sqs2 + exp( -( (x22-m2)/sqs2) ) ) ) 
           			   p2 <- exp(-(exp(-(x22 - m2)/sqs2))) 
           }
           if(mar2=="LO")          {d.x2 <- dlogis(x22, m2, sqs2 )       
                                    p2 <- plogis(x22, m2, sqs2 ) 
           }
           if(mar2=="WEI")         {d.x2 <- sqs2/exp(m2)*(x22/exp(m2))^(sqs2-1) * exp(-(x22/exp(m2))^sqs2)            
           			   p2 <- 1 - exp(-(x22/exp(m2))^sqs2) 
           }
           if(mar2=="WEI2")         {d.x2 <- sqs2*exp(m2)*x22^(sqs2-1)*exp(-exp(m2)*x22^(sqs2))              
           			   p2 <- 1  - exp( - ( x22/(exp(m2)^( -1/sqs2 ) ) )^sqs2)  
           }          
           if(mar2=="iG")          {d.x2 <- exp(-0.5 * log(2 * pi) - log(sqs2) - (3/2) * log(x22) - ((x22 - exp(m2))^2)/(2 * sqs2^2 * (exp(m2)^2) * x22))          
                                    p2 <- pnorm(((x22/exp(m2)) - 1)/(sqs2 * sqrt(x22))) + exp(2/(exp(m2)*sqs2^2))* pnorm(-((x22/exp(m2)) + 1)/(sqs2 * sqrt(x22)))
           }
           if(mar2=="GA")          {d.x2 <- dgamma(x22, shape = 1/sqs2, scale = exp(m2) * sqs2)          
                                    p2 <- pgamma(x22, shape = 1/sqs2, scale = exp(m2) * sqs2)
           }                        
           if(mar2=="iGA")         {d.x2 <- exp(1/sqs2 * m2 + 1/sqs2 * log(1/sqs2 + 1) - lgamma(1/sqs2) - (1/sqs2 + 1) * log(x22) - ((exp(m2) * (1/sqs2 + 1))/x22))          
                                    p2 <- 1-pgamma(((exp(m2) * (1/sqs2 + 1))/x22), shape = 1/sqs2, scale=1)
           }   
           if(mar2=="DAGUM")       {d.x2 <- sqs2*nu/x22 *( ( (x22/exp(m2)^(-1/sqs2))^(sqs2*nu) )/  (  ((x22/exp(m2)^(-1/sqs2))^sqs2 + 1)^(nu+1) )  )
                                    p2 <- (1 + (x22/exp(m2)^(-1/sqs2))^-sqs2)^-nu
           }          
           
           
           list(p2 = p2, d.x2 = d.x2)
           
           
           }
          
 
 
 Cop.pdf <- function (u1, u2, par1, fam) {
 
 if(fam == 1) { # N
     t1 = qnorm(u1)
     t2 = qnorm(u2)
     res <- 1/sqrt(1 - par1^2) * exp(-(par1^2 * (t1^2 + t2^2) - 2 * par1 * t1 * t2)/(2 * (1 - par1^2)))
 } 
 if(fam == 14) { # F
     res <- (par1 * (exp(par1) - 1) * exp(par1 * u2 + par1 * u1 + par1))/(exp(par1 * u2 + par1 * u1) - exp(par1 * u2 + par1) - exp(par1 * u1 + par1) + exp(par1))^2
     }  
     
 if(fam %in% c(2, 6, 10)  ){ # 0
     theta = par1
     d1 = u1
     d2 = u2}
 if(fam %in% c(3, 7, 11)){# 90
     theta = -par1
     d1 = 1 - u1
     d2 = u2}  
 if(fam %in% c(4, 8, 12)){# 180
     theta = par1
     d1 = 1 - u1
     d2 = 1 - u2}     
 if(fam %in% c(5, 9, 13)){# 270
     theta = -par1
     d1 = u1
     d2 = 1 - u2}    
     
 if(fam %in% c(2:5)) { # C
     res <- (1 + theta) * (d1 * d2)^(-1 - theta) * (d1^(-theta) + d2^(-theta) - 1)^(-2 - 1/theta)
 }
 if(fam %in% c(6:9)) { # G    
     t1 <- (-log(d1))^(theta) + (-log(d2))^(theta)
     t2 <- exp(-t1^(1/theta))
     res <- t2/(d1 * d2) * t1^(-2 + 2/theta) * (log(d1) * log(d2))^(theta - 1) * (1 + (theta - 1) * t1^(-1/theta))
 }
 if(fam %in% c(10:13)) { # J 
     res <- ((1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))^(1/(theta) - 2) * (1 - d1)^(theta - 1) * (1 - d2)^(theta - 1) * (theta - 1 + (1 - d1)^(theta) + (1 - d2)^(theta) - (1 - d1)^(theta) * (1 - d2)^(theta))
 }
 
 res
 
 }
 
 

  Cplot <- function (fam, par1, mar2, m2, sqs2, nu, resp, ...){   
  
  
          
          size <- 100

          
          
          
                             x1 <- seq(from = xlim[1], to = xlim[2], length.out = size)               
          if(mar2=="probit") x2 <- seq(from = ylim[1], to = ylim[2], length.out = size)              #; levels <- c(0.01,0.05,0.1,0.15,0.2) } 
        if(mar2 != "probit") x2 <- seq(from = min(resp), to = max(resp), length.out = size)          #; levels <- c(0.01,0.05,0.1,0.15,0.2) } 
          
                            x11 <- rep(x1, each = size)
                            x22 <- rep(x2, times = size)
          
        if(mar2 != "probit") x2 <- x2 - mean(x2) # seq(from = ylim[1], to = ylim[2], length.out = size)
          
                                   d.x1 <- dnorm(x11)
                                   p1   <- pnorm(x11) 
          
          
          
          resf <- funcm2(x22, m2, sqs2, nu, mar2)
          
          
          
          
          md <- Cop.pdf(p1, resf$p2, par1, fam)*d.x1*resf$d.x2    
          z  <- matrix(data = md, nrow = size, byrow = TRUE)
          
                    
          filled.contour(x1, x2, z, color = topo.colors, nlevels = 16, ...) 
          

          
 }
 


  c1   <- c("N","GU","rGU","LO")
  c2   <- c("LN","WEI","WEI2","iG","GA","iGA","DAGUM") 

  if(object$margins[2] %in% c(c1,c2)){

  if(object$margins[2] %in% c2) test.resp <- seq(from = 0.0001, to = max(object$y2), length.out = 1000)
  if(object$margins[2] %in% c1) test.resp <- seq(from = min(object$y2), to = max(object$y2), length.out = 1000)
  
  test.dens <- round( funcm2(test.resp, m2, sqs2, nu, object$VC$margins[2])$d.x2, 2)
  resp <- test.resp[which(test.dens != 0)]
 
  } else resp <- object$y2
  
  
  
 
Cplot(fam = object$nC, par1 = par1, main = cop, ylab = ylab, xlab = xlab, mar2 = object$margins[2], 
      m2 = m2, sqs2 = sqs2, nu = nu, resp = resp, ...)  
      
      
      
      
  
 }
 
 
 
 
 
 
 
 
 
 
 
 
 
rm(bs, SE, Vb, epds, sigma2.st, sigma2, est.RHOb, et1s, et2s, p1s, p2s, p11s, p10s, p00s, p01s, ORs, GMs, XX, Xt, V) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], 
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], 
              n=n, theta=object$theta.a, sigma2=object$sigma2.a, nu=object$nu.a, OR = object$OR, GM = object$GM, 
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula,
              t.edf=object$t.edf, CItheta=CIrs, CIsig2=CIsig2, CInu=CInu,  
              n.sel=n.sel, CIor = CIor, CIgm = CIgm, 
              BivD=object$BivD, margins = object$margins, 
              Model=object$Model,
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6
              )
  class(res) <- "summary.SemiParBIVProbit"
      
                                        

res

}

          #cp <- c("#8E063B", "#A63945", "#BC584D" ,"#CF7355", "#DF8B5B" ,"#EAA162", "#F2B468" ,"#F6C56F", "#F6D277", "#F2DD80", "#EBE48B" ,"#E2E6BD")
          #cp <- c("#3A58B6", "#6458B9" ,"#805ABB", "#965EBB", "#A763BB" ,"#B66AB9", "#C273B8", "#CC7CB6" ,"#D586B4", "#DD91B3", "#E39CB3" ,"#E8A8B4")
          #cp=c("#023FA5", "#495DA8", "#6B77B2" ,"#868FBD", "#9EA4C6" ,"#B1B5CE", "#C1C4D5","#CED0DA", "#D8D8DE", "#DEDEE1", "#E1E1E2", "#E2E2E2")
          #cp=c("#023FA5" ,"#7D87B9" ,"#BEC1D4", "#E2E2E2", "#D6BCC0" ,"#BB7784" ,"#8E063B")
          #cp=c("#3A58B6" ,"#5B58B9" ,"#7259BA", "#855BBB" ,"#945EBB" ,"#A261BB" ,"#AD66BA", "#B76BB9", "#C072B8", "#C878B7", "#D07FB5", "#D687B4", "#DB8FB3", "#E097B3","#E49FB3", "#E8A8B4")
          #col=cp
           #heat.colors(0:n.col/n.col)color = gray(0:1) cm.colors topo.colors  terrain.colors rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
           # rainbow(n, start=.7, end=.1)   rainbow(10, start=.7, end=.1)
          #color.palette = gray(0:n.col/n.col)