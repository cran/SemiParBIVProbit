SemiParBIVProbit.fit.post <- function(SemiParFit, formula.eq2, data, 
                                      Model, VC, 
                                      PL, eqPL, valPL, fitPL,  
                                      qu.mag=NULL, gam1, gam2, gam3){

non.sel.dd <- lambda1s <- lambda2s <- eta1S <- eta2S <- athrhoS <- rho <- theta <- edf <- rho.a <- etad <- theta.a <- NULL

xi1 <- xi2 <- 1
   
He <- SemiParFit$fit$hessian
logLik <- -SemiParFit$fit$l

epsilon <- sqrt(.Machine$double.eps)
                                                                                                         
    He.eig <- eigen(He, symmetric=TRUE)
    if(min(He.eig$values) < epsilon) He.eig$values[which(He.eig$values < epsilon)] <- 0.0000001
    Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)   

                                     
if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0) && VC$fp==FALSE){ 
                                          HeSh <- He - SemiParFit$fit$S.h
                                          F <- Vb%*%HeSh
                                        }else{ HeSh <- He; F <- diag(rep(1,dim(Vb)[1])) } 
t.edf <- sum(diag(F))


dimnames(SemiParFit$fit$hessian)[[1]] <- dimnames(SemiParFit$fit$hessian)[[2]] <- dimnames(Vb)[[1]] <- dimnames(Vb)[[2]] <- dimnames(HeSh)[[1]] <- dimnames(HeSh)[[2]] <- dimnames(F)[[1]] <- dimnames(F)[[2]] <- dimnames(He)[[1]] <- dimnames(He)[[2]] <- names(SemiParFit$fit$argument)   

  
if(VC$hess == FALSE) SemiParFit$fit$Fisher <- SemiParFit$fit$hessian



  if( is.null(VC$X3) ){ if( VC$BivD == "N" ) dep <- SemiParFit$fit$argument["athrho"] else dep <- SemiParFit$fit$argument["theta.star"]  } 

  if(!is.null(VC$X3) ) dep <- SemiParFit$fit$etad  



  if( VC$BivD == "N" ) {rho <- tanh(dep); rho <- ifelse(rho == -1, -0.9999999, rho)
                                          rho <- ifelse(rho == 1 ,  0.9999999, rho) }#; names(rho) <- rep("rho",length(rho))} 
  
  if( VC$BivD != "N" ) {

                    th.st <- dep 

   if(VC$BivD=="F") theta <- th.st + epsilon 
  
   if(VC$BivD %in% c("C0", "C180") ) theta <- exp(th.st) + epsilon
   if(VC$BivD %in% c("C90","C270") ) theta <- -(exp(th.st) + epsilon)

   if(VC$BivD %in% c("J0", "J180","G0", "G180") ) theta <-    1 + exp(th.st) + epsilon
   if(VC$BivD %in% c("J90","J270","G90","G270") ) theta <- -( 1 + exp(th.st) + epsilon )
       
   theta <- ifelse(theta == Inf ,  8.218407e+307, theta )
   theta <- ifelse(theta == -Inf, -8.218407e+307, theta )    
       
   #names(theta) <- rep("theta",length(theta))
   
        }
        
        
if( VC$BivD == "N" ) rho.a <- mean(rho) else theta.a <- mean(theta) 


  
  if((PL=="PP" || PL=="RPP")){
  
       if(fitPL=="fixed"){
                     xi1 <- exp(valPL[1])
                     xi2 <- exp(valPL[2])
                         }else{
  
              if(eqPL=="both"){
              xi1 <- exp(SemiParFit$fit$argument["xi1.star"])
              xi2 <- exp(SemiParFit$fit$argument["xi2.star"])
                              }
              if(eqPL=="first"){
              xi1 <- exp(SemiParFit$fit$argument["xi1.star"])
              xi2 <- 1
                               }
              if(eqPL=="second"){
              xi1 <- 1
              xi2 <- exp(SemiParFit$fit$argument["xi2.star"])
                                }   
                               }        
              }

  if(PL=="SN"){
  
     if(fitPL=="fixed"){
                   xi1 <- valPL[1]
                   xi2 <- valPL[2]
                       }else{
  
  
              if(eqPL=="both"){
              xi1 <- SemiParFit$fit$argument["xi1.star"]
              xi2 <- SemiParFit$fit$argument["xi2.star"]
                              }
              if(eqPL=="first"){
              xi1 <- SemiParFit$fit$argument["xi1.star"]
              xi2 <- 0
                               }
              if(eqPL=="second"){
              xi1 <- 0
              xi2 <- SemiParFit$fit$argument["xi2.star"]
                                } 
                             }
              }
              
  names(xi1) <- "xi1"; names(xi2) <- "xi2"               




###
## can be made more efficient by avoiding eta1 etc in some models?
###


  if(Model=="BSS"){

  eta1      <- SemiParFit$fit$eta1
  resp      <- rep(1,length(gam1$y))
  fs        <- as.formula( paste("resp","~",formula.eq2[3],sep="") ) 

  non.sel.dd <- gam(fs, data=data, fit = FALSE)$X[SemiParFit$fit$good,] 
  
  if(VC$gc.l == TRUE) gc()  
  
  if( is.null(VC$X3) ) ll <- length(SemiParFit$fit$argument)
  if(!is.null(VC$X3) ) ll <- (VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)
 
  param <- SemiParFit$fit$argument[-c(1:VC$X1.d2,ll)] 
  
  if(length(param)!=dim(non.sel.dd)[2]){
  posit <- which(names(VC$X2[1,])%in%names(param))
  non.sel.dd <- non.sel.dd[,posit]

  }
  
  eta2 <- non.sel.dd%*%param

  SemiParFit$fit$eta2 <- eta2
   
}



if(Model=="BSS" || Model=="BPO"){

  if(Model=="BPO") {eta1 <- SemiParFit$fit$eta1; eta2 <- SemiParFit$fit$eta2} 

  p1 <- pmax(pnorm(eta1), epsilon) 
  p2 <- pmax(pnorm(eta2), epsilon ) 
  p1 <- ifelse(p1==1,0.9999999,p1)
  p2 <- ifelse(p2==1,0.9999999,p2)
  
  
  if(VC$BivD=="N") theta <- rho 
  
   p11 <- pmax( BiCDF(p1, p2, VC$nC, theta), epsilon )
  
   SemiParFit$fit$p10 <- p1 - p11
   SemiParFit$fit$p11 <- p11
   SemiParFit$fit$p00 <- (1 - p2) - ( p1 - p11 )
   SemiParFit$fit$p01 <- p2 - p11
   
   SemiParFit$fit$p1 <- p1
   SemiParFit$fit$p2 <- p2

}


######################
# Association measures
######################


p00 <- SemiParFit$fit$p00 
p01 <- SemiParFit$fit$p01 
p11 <- SemiParFit$fit$p11 
p10 <- SemiParFit$fit$p10 

p1 <- SemiParFit$fit$p1
p2 <- SemiParFit$fit$p2

OR <- (p00*p11)/(p01*p10)

OR  <- ifelse(OR  ==  Inf,  8.218407e+307, OR ) 
OR  <- ifelse(OR  == -Inf, -8.218407e+307, OR ) 

GM <- mean((OR - 1)/(OR + 1))
OR <- mean(OR)


rm(p00,p01,p10,p11,p1,p2)

if(VC$gc.l == TRUE) gc()  




if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0) ){

  edf <- list(0,0,0)
        
     for(i in 1:3){

       if(i==1) {mm <- VC$l.sp1; if(mm==0) next}
       if(i==2) {mm <- VC$l.sp2; if(mm==0) next} 
       if(i==3) {mm <- VC$l.sp3; if(mm==0) break} 

          for(k in 1:mm){

              if(i==1){ gam <- gam1; ind <- (gam$smooth[[k]]$first.para):(gam$smooth[[k]]$last.para) } 
              if(i==2){ gam <- gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 } 
              if(i==3){ gam <- gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + VC$X1.d2 + VC$X2.d2 } 
              
	      edf[[i]][k] <- sum(diag(F)[ind])
                        }
                  }
         
  if(VC$l.sp1!=0) names(edf[[1]]) <- names(gam1$sp)  
  if(VC$l.sp2!=0) names(edf[[2]]) <- names(gam2$sp)   
  if(VC$l.sp3!=0) names(edf[[3]]) <- names(gam3$sp)   

}
  
 
  sp <- SemiParFit$sp 
  
  
  if( VC$fp==FALSE ){ 
  
  if(VC$l.sp1!=0 && VC$l.sp2!=0 && VC$l.sp3!=0) names(sp) <- c(c(paste(names(gam1$sp),".eq1",sep="")),c(paste(names(gam2$sp),".eq2",sep="")),c(paste(names(gam3$sp),".eq3",sep="")))
  if(VC$l.sp1!=0 && VC$l.sp2!=0 && VC$l.sp3==0) names(sp) <- c(c(paste(names(gam1$sp),".eq1",sep="")),c(paste(names(gam2$sp),".eq2",sep="")))
  if(VC$l.sp1!=0 && VC$l.sp2==0 && VC$l.sp3==0) names(sp) <- paste(names(gam1$sp),".eq1",sep="")
  if(VC$l.sp1==0 && VC$l.sp2!=0 && VC$l.sp3!=0) names(sp) <- c(c(paste(names(gam2$sp),".eq2",sep="")),c(paste(names(gam3$sp),".eq3",sep="")))
  if(VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3!=0) names(sp) <- paste(names(gam3$sp),".eq3",sep="")
  if(VC$l.sp1==0 && VC$l.sp2!=0 && VC$l.sp3==0) names(sp) <- paste(names(gam2$sp),".eq2",sep="")
  if(VC$l.sp1!=0 && VC$l.sp2==0 && VC$l.sp3!=0) names(sp) <- c(c(paste(names(gam1$sp),".eq1",sep="")),c(paste(names(gam3$sp),".eq3",sep="")))
  
  }
  
  
  
  
  
  if(PL!="P"){

    wna <- which(is.na(names(sp)))

	if(eqPL=="both")   names(sp)[wna] <- c("xi1","xi2")      
        if(eqPL=="first")  names(sp)[wna] <- c("xi1")     
        if(eqPL=="second") names(sp)[wna] <- c("xi2")   

  }
  
  
  if( VC$BivD == "N" ){ if( length(rho)==1 ) rho.a <- rho} else { if( length(theta)==1 ) theta.a <- theta }  
  
  
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, F = F, t.edf = t.edf,
                      edf1 = edf[[1]], edf2 = edf[[2]], edf3 = edf[[3]], rho = rho, theta = theta, rho.a = rho.a, theta.a = theta.a,
                      xi1 = xi1, xi2 = xi2, sp = sp, OR = OR, GM = GM, 
                      X2s = non.sel.dd) 

}

















