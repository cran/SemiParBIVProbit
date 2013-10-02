SemiParBIVProbit.fit.post <- function(SemiParFit, formula.eq2, selection, data, RE, RE.type, BivD, nu, gev.eq1, gev.eq2, shape1, shape2, nC, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, pPen1, pPen2, l.sp1, l.sp2,
                                      qu.mag=NULL, gam1, gam2, gp1, gp2, fp, weights, K=NULL, n, N=NULL, cuid=NULL, uidf=NULL){


 NP.qu.int <- eta1S <- eta2S <- athrhoS <- sigma1 <- rho.u <- sigma2 <- rho <- theta <- NULL

    if(selection==FALSE && RE==TRUE && RE.type=="NP"){
    
    	NP.qu.int <- NP.qu(SemiParFit, BivD, nC, nu, gev.eq1, gev.eq2, shape1, shape2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, pPen1, pPen2, qu.mag, gp1, gp2, fp, l.sp1, l.sp2, weights, K, n, N, cuid, uidf)
    	He <- NP.qu.int$He; logLik <- NP.qu.int$logLik 

    	SemiParFit$fit$eta1 <- NP.qu.int$Eb.u1 + SemiParFit$fit$eta1 
    	SemiParFit$fit$eta2 <- NP.qu.int$Eb.u2 + SemiParFit$fit$eta2

        X1 <- model.matrix(gam1); X2 <- model.matrix(gam2)

                                      }else{ 
                                               He <- SemiParFit$fit$hessian; logLik <- -SemiParFit$fit$l; K <- NULL
                                           }
    

    if(selection==FALSE && RE==TRUE && RE.type=="N") NP.qu.int <- EBNormal(SemiParFit$fit$argument,X1,X2,X1.d2,X2.d2,N,cuid,uidf,tol=1e-6,M=50)


    He.eig <- eigen(He,symmetric=TRUE)
    k.e <- sum(as.numeric(He.eig$val<sqrt(.Machine$double.eps)))
    
     if(k.e!=0){
      ind.e <- (length(He.eig$val)-(k.e-1)):length(He.eig$val)
      min.e <- min(He.eig$val[1:(ind.e[1]-1)])
      for(i in 1:k.e) He.eig$val[ind.e[i]] <- min.e/10^i  
      Vb <- He.eig$vec%*%tcrossprod(diag(1/He.eig$val),He.eig$vec)      
     }else{
      Vb <- He.eig$vec%*%tcrossprod(diag(1/He.eig$val),He.eig$vec) 
    }
          

                                    
  if( (l.sp1!=0 || l.sp2!=0) && fp==FALSE){ if(RE==TRUE && RE.type=="NP") SemiParFit$fit$S.h <- adiag(SemiParFit$fit$S.h,matrix(0,K-1,K-1)) 

                                                        HeSh <- He - SemiParFit$fit$S.h; F <- Vb%*%HeSh

                                       }else{ HeSh <- He; F <- diag(rep(1,dim(Vb)[1])) } 


    
  t.edf <- sum(diag(F))


  if( BivD %in% c("N","T") ) {rho <- tanh(SemiParFit$fit$argument["athrho"]); names(rho) <- "rho"; KeT <- BiCopPar2Tau(nC,rho,par2=nu)} 
  else{

   if(BivD=="F") theta <- SemiParFit$fit$argument["theta.star"]  
  
   if(BivD %in% c("C0", "C180") ) theta <- exp(SemiParFit$fit$argument["theta.star"])
   if(BivD %in% c("C90","C270") ) theta <- -exp(SemiParFit$fit$argument["theta.star"])

   if(BivD %in% c("J0", "J180","G0", "G180") ) theta <-    1 + exp(SemiParFit$fit$argument["theta.star"])
   if(BivD %in% c("J90","J270","G90","G270") ) theta <- -( 1 + exp(SemiParFit$fit$argument["theta.star"]) )
  
   names(theta) <- "theta"; KeT <- BiCopPar2Tau(nC,theta,par2=nu)

   }

  names(KeT) <- "k.tau"


  if(RE==TRUE && RE.type=="N"){

          sigma1=exp(SemiParFit$fit$argument["log.sigma1"]); names(sigma1) <- "sigma1"
          rho.u=tanh(SemiParFit$fit$argument["athrho.u"]); names(rho) <- "rho.u"
          sigma2=exp(SemiParFit$fit$argument["log.sigma2"]); names(sigma2) <- "sigma2"

             }

  if(selection==TRUE){

  eta1      <- SemiParFit$fit$eta1
  resp      <- rep(1,length(eta1))
  fs        <- as.formula( paste("resp","~",formula.eq2[3],sep="") ) 
  non.sel.d <- gam(fs, data=data[SemiParFit$fit$good,], fit = FALSE)$X
  param     <- SemiParFit$fit$argument[-c(1:length(gam1$coef),length(SemiParFit$fit$argument))]
  eta2      <- non.sel.d%*%param

  if(gev.eq1==FALSE && gev.eq2==TRUE) {p1 <- pnorm(eta1);             p2 <- 1-pgev(-eta2,shape=shape2)}   
  if(gev.eq1==TRUE  && gev.eq2==FALSE){p1 <- 1-pgev(-eta1,shape=shape1); p2 <- pnorm(eta2)            }
  if(gev.eq1==TRUE  && gev.eq2==TRUE) {p1 <- 1-pgev(-eta1,shape=shape1); p2 <- 1-pgev(-eta2,shape=shape2)}
  if(gev.eq1==FALSE && gev.eq2==FALSE){p1 <- pnorm(eta1);             p2 <- pnorm(eta2)            } 

  if(BivD=="N") p11 <- pmax( abs(pbinorm( qnorm(p1), qnorm(p2), cov12=rho)), 1000*.Machine$double.eps ) 
  else{ 
   if(BivD=="T") theta <- rho
   p11 <- BiCopCDF(p1,p2, nC, par=theta,par2=nu)
   }

   SemiParFit$fit$p00 <- (1-p2) - ( p1 - p11 )
   SemiParFit$fit$p01 <- p2 - p11

  bs    <- rmvnorm(1000, mean = SemiParFit$fit$argument, sigma=Vb, method="svd")
  eta1S <- X1%*%t(bs[,1:length(gam1$coef)]) 
  eta2S <- non.sel.d%*%t(bs[,-c(1:length(gam1$coef),length(SemiParFit$fit$argument))]) 
  athrhoS  <- bs[,length(SemiParFit$fit$argument)] 


}

  edf <- NULL
  l.sp11 <- length(gam1$smooth)
  l.sp22 <- length(gam2$smooth) 

  if( (l.sp11!=0 || l.sp22!=0) ){

  Kk1 <- Kk2 <- 0; if(RE==TRUE && RE.type=="NP"){ Kk1 <- K - 1; Kk2 <- Kk1 + 1}  
  Ks1 <- Ks2 <- 0; if(RE==TRUE && RE.type=="N"){  Ks1 <- 1; Ks2 <- 2}

  edf <- list(0,0)
        
     for(i in 1:2){

       if(i==1) {mm <- l.sp11; if(mm==0) next}
       if(i==2) {mm <- l.sp22; if(mm==0) break} 

          for(k in 1:mm){

              if(i==1){gam <- gam1; ind <- (gam$smooth[[k]]$first.para+Kk1+Ks1):(gam$smooth[[k]]$last.para+Kk1+Ks1)} 
                  else{gam <- gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para+Kk1+Ks1)+X1.d2+Kk2+Ks2} 
	      edf[[i]][k] <- sum(diag(F)[ind])
                        }
                  }
                  
  if(length(gam1$paraPen)!=0 && l.sp1>1)  names(edf[[1]]) <- names(gam1$sp)[-1]
  if(length(gam1$paraPen)==0 && l.sp1!=0) names(edf[[1]]) <- names(gam1$sp)  

  if(length(gam2$paraPen)!=0 && l.sp2>1)  names(edf[[2]]) <- names(gam2$sp)[-1] 
  if(length(gam2$paraPen)==0 && l.sp2!=0) names(edf[[2]]) <- names(gam2$sp) 
  
  }
 

  sp <- SemiParFit$sp 
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE) names(sp) <- c(c(paste(names(gam1$sp),".eq1",sep="")),c(paste(names(gam2$sp),".eq2",sep="")))
  if(l.sp1==0 && l.sp2!=0 && fp==FALSE) names(sp) <- paste(names(gam2$sp),".eq2",sep="")
  if(l.sp1!=0 && l.sp2==0 && fp==FALSE) names(sp) <- paste(names(gam1$sp),".eq1",sep="")
        

   
                 list(SemiParFit=SemiParFit,NP.qu.int=NP.qu.int,He=He,
                      logLik=logLik,X1=X1,X2=X2,K=K,Vb=Vb,HeSh=HeSh,F=F,t.edf=t.edf,edf1=edf[[1]],edf2=edf[[2]],rho=rho,theta=theta,KeT=KeT,
                      sigma1=sigma1,rho.u=rho.u,sigma2=sigma2, 
                      eta1S=eta1S,eta2S=eta2S,ass.pS=athrhoS,nC=nC,sp=sp)



}

















