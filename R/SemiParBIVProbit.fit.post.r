SemiParBIVProbit.fit.post <- function(SemiParFit, formula.eq2, data, 
                                      Model, VC, 
                                      PL, eqPL, valPL, fitPL,  
                                      qu.mag=NULL, gam1, gam2){

non.sel.dd <- lambda1s <- lambda2s <- eta1S <- eta2S <- athrhoS <- rho <- theta <- edf <- NULL
xi1 <- xi2 <- 1
   
He <- SemiParFit$fit$hessian
logLik <- -SemiParFit$fit$l
                                           
    He.eig <- eigen(He, symmetric=TRUE)
    
    if(min(He.eig$values) < .Machine$double.eps){ He <- as.matrix( nearPD( He, ensureSymmetry = FALSE )$mat )
                                                  He.eig <- eigen(He, symmetric=TRUE) }
     
Vb <- He.eig$vec%*%tcrossprod(diag(1/He.eig$val),He.eig$vec)      

                                     
if( (VC$l.sp1!=0 || VC$l.sp2!=0) && VC$fp==FALSE){ 
                                          HeSh <- He - SemiParFit$fit$S.h
                                          F <- Vb%*%HeSh
                                        }else{ HeSh <- He; F <- diag(rep(1,dim(Vb)[1])) } 
t.edf <- sum(diag(F))
  

  if( VC$BivD %in% c("N","T") ) {rho <- tanh(SemiParFit$fit$argument["athrho"]); names(rho) <- "rho"; KeT <- BiCopPar2Tau(VC$nC,rho,par2=VC$nu)} 
  else{

   th.st <- SemiParFit$fit$argument["theta.star"]  

   if(VC$BivD=="F") theta <- th.st
  
   if(VC$BivD %in% c("C0", "C180") ) theta <- exp(th.st)
   if(VC$BivD %in% c("C90","C270") ) theta <- -exp(th.st)

   if(VC$BivD %in% c("J0", "J180","G0", "G180") ) theta <-    1 + exp(th.st)
   if(VC$BivD %in% c("J90","J270","G90","G270") ) theta <- -( 1 + exp(th.st) )
       
   names(theta) <- "theta"; KeT <- BiCopPar2Tau(VC$nC,theta)
   
        }

  names(KeT) <- "k.tau"
  
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







  if(Model=="BSS"){

  eta1      <- SemiParFit$fit$eta1
  resp      <- rep(1,length(gam1$y))
  fs        <- as.formula( paste("resp","~",formula.eq2[3],sep="") ) 

  non.sel.dd <- gam(fs, data=data, fit = FALSE)$X[SemiParFit$fit$good,] 
  non.sel.d  <- gam(fs, data=data, fit = FALSE)$X 
  
  ll <- length(SemiParFit$fit$argument)
  
  param <- SemiParFit$fit$argument[-c(1:length(gam1$coef),ll)] 
  
  if(length(param)!=dim(non.sel.dd)[2]){
  posit <- which(names(VC$X2[1,])%in%names(param))
  non.sel.dd <- non.sel.dd[,posit]
  non.sel.d  <- non.sel.d[,posit]
  }
  
  eta2 <- non.sel.dd%*%param

  SemiParFit$fit$eta2 <- eta2
   
}



if(Model=="BSS" || Model=="BPO"){

  if(Model=="BPO") {eta1 <- SemiParFit$fit$eta1; eta2 <- SemiParFit$fit$eta2} 

  p1 <- pmax(pnorm(eta1), 1000*.Machine$double.eps ) 
  p2 <- pmax(pnorm(eta2), 1000*.Machine$double.eps ) 
  p1 <- ifelse(p1==1,0.9999999999999999,p1)
  p2 <- ifelse(p2==1,0.9999999999999999,p2)
  
  if(VC$BivD=="N") p11 <- pmax( pbinorm( eta1, eta2, cov12=rho), 1000*.Machine$double.eps ) 
  else{ 
   if(VC$BivD=="T") theta <- rho
   p11 <- pmax(BiCopCDF(p1,p2, VC$nC, par=theta, par2=VC$nu), 1000*.Machine$double.eps ) 
   }

   SemiParFit$fit$p10 <- p1 - p11
   SemiParFit$fit$p11 <- p11
   SemiParFit$fit$p00 <- (1 - p2) - ( p1 - p11 )
   SemiParFit$fit$p01 <- p2 - p11

}


  l.sp11 <- length(gam1$smooth)
  l.sp22 <- length(gam2$smooth) 

  if( (l.sp11!=0 || l.sp22!=0) ){

  edf <- list(0,0)
        
     for(i in 1:2){

       if(i==1) {mm <- l.sp11; if(mm==0) next}
       if(i==2) {mm <- l.sp22; if(mm==0) break} 

          for(k in 1:mm){

              if(i==1){gam <- gam1; ind <- (gam$smooth[[k]]$first.para):(gam$smooth[[k]]$last.para)} 
                  else{gam <- gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para)+VC$X1.d2} 
	      edf[[i]][k] <- sum(diag(F)[ind])
                        }
                  }
                  
  if(length(gam1$paraPen)!=0 && VC$l.sp1>1)  names(edf[[1]]) <- names(gam1$sp)[-1]
  if(length(gam1$paraPen)==0 && VC$l.sp1!=0) names(edf[[1]]) <- names(gam1$sp)  

  if(length(gam2$paraPen)!=0 && VC$l.sp2>1)  names(edf[[2]]) <- names(gam2$sp)[-1] 
  if(length(gam2$paraPen)==0 && VC$l.sp2!=0) names(edf[[2]]) <- names(gam2$sp) 
  
  }
 

  sp <- SemiParFit$sp 
  if(VC$l.sp1!=0 && VC$l.sp2!=0 && VC$fp==FALSE) names(sp) <- c(c(paste(names(gam1$sp),".eq1",sep="")),c(paste(names(gam2$sp),".eq2",sep="")))
  if(VC$l.sp1==0 && VC$l.sp2!=0 && VC$fp==FALSE) names(sp) <- paste(names(gam2$sp),".eq2",sep="")
  if(VC$l.sp1!=0 && VC$l.sp2==0 && VC$fp==FALSE) names(sp) <- paste(names(gam1$sp),".eq1",sep="")
  
  if(PL!="P"){

    wna <- which(is.na(names(sp)))

	if(eqPL=="both")   names(sp)[wna] <- c("xi1","xi2")      
        if(eqPL=="first")  names(sp)[wna] <- c("xi1")     
        if(eqPL=="second") names(sp)[wna] <- c("xi2")   

  }

   
                 list(SemiParFit = SemiParFit, He = He, logLik = logLik, Vb = Vb, HeSh = HeSh, F = F, t.edf = t.edf,
                      edf1 = edf[[1]], edf2 = edf[[2]], rho = rho, theta = theta, KeT = KeT,
                      xi1 = xi1, xi2 = xi2, sp = sp, 
                      X2s = non.sel.dd)

}

















