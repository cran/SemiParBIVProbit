SemiParBIVProbit.fit.post <- function(SemiParFit, formula.eq2, selection, data, RE, RE.type, y1, y2, q1, q2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, l.sp1, l.sp2,
                                      qu.mag=NULL, gam1, gam2, gp1, gp2, fp, weights, K=NULL, n, N=NULL, cuid=NULL, uidf=NULL){


 NP.qu.int <- eta1S <- eta2S <- athrhoS <- sigma1 <- rho.u <- sigma2 <- NULL

    if(selection==FALSE && RE==TRUE && RE.type=="NP"){
    
    	NP.qu.int <- NP.qu(SemiParFit, y1, y2, q1, q2, y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1, X1, X2, X1.d2, X2.d2, qu.mag, gp1, gp2, fp, l.sp1, l.sp2, weights, K, n, N, cuid, uidf)
    	He <- NP.qu.int$He; logL <- NP.qu.int$logL

    	SemiParFit$fit$eta1 <- NP.qu.int$Eb.u1 + SemiParFit$fit$eta1 
    	SemiParFit$fit$eta2 <- NP.qu.int$Eb.u2 + SemiParFit$fit$eta2

        X1 <- model.matrix(gam1); X2 <- model.matrix(gam2)

                                      }else{ 
                                               He <- SemiParFit$fit$hessian; logL <- SemiParFit$fit$l; K <- NULL
                                           }
    



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
  rho <- tanh(SemiParFit$fit$argument["athrho"]); names(rho) <- "rho"
  
  if(RE==TRUE && RE.type=="N"){

          sigma1=exp(SemiParFit$fit$argument["log.sigma1"]); names(sigma1) <- "sigma1"
          rho.u=tanh(SemiParFit$fit$argument["athrho.u"]); names(rho) <- "rho.u"
          sigma2=exp(SemiParFit$fit$argument["log.sigma2"]); names(sigma2) <- "sigma2"

             }

  if(selection==TRUE){

  eta1      <- SemiParFit$fit$eta1
  fs <- as.formula( paste("y2","~",formula.eq2[3],sep="") ) 
  non.sel.d <- gam(fs, data=data, fit = FALSE)$X
  param     <- SemiParFit$fit$argument[-c(1:length(gam1$coef),length(SemiParFit$fit$argument))]
  eta2      <- non.sel.d%*%param

  SemiParFit$fit$p00 <- pnorm2(-eta1,-eta2,cov12=rho) 
  SemiParFit$fit$p01 <- pnorm2(-eta1,eta2,cov12=-rho)

  bs    <- rmvnorm(1000, mean = SemiParFit$fit$argument, sigma=Vb, method="svd")
  eta1S <- X1%*%t(bs[,1:length(gam1$coef)]) 
  eta2S <- non.sel.d%*%t(bs[,-c(1:length(gam1$coef),length(SemiParFit$fit$argument))]) 
  athrhoS  <- bs[,length(SemiParFit$fit$argument)] 

  }


   
                 list(SemiParFit=SemiParFit,NP.qu.int=NP.qu.int,He=He,
                      logL=logL,X1=X1,X2=X2,K=K,Vb=Vb,HeSh=HeSh,F=F,t.edf=t.edf,rho=rho,sigma1=sigma1,rho.u=rho.u,sigma2=sigma2, 
                      eta1S=eta1S,eta2S=eta2S,athrhoS=athrhoS)



}

















