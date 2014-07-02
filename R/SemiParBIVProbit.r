SemiParBIVProbit <- function(formula.eq1, formula.eq2, data=list(), weights=NULL,  
                             start.v=NULL, BivD="N", nu=3, PL="P", eqPL="both", valPL=c(0,0), fitPL="pLiksp", # 
                             spPL=c(0.01,0.01),
                             selection=FALSE, H.n=TRUE, gamma=1, aut.sp=TRUE, fp=FALSE,
                             pPen1 = NULL, pPen2 = NULL, 
                             rinit=1, rmax=100, fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     iterlimsp=50, pr.tolsp=1e-6, 
                             control.sp=list(maxit=50,tol=1e-6,step.half=25,rank.tol=sqrt(.Machine$double.eps)) ){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################

  uidf <- sp <- qu.mag <- sp.xi1 <- sp.xi2 <- NULL 
  
  opc <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","T")
  bb1 <- c("BB1.0","BB1.90","BB1.180","BB1.270")
  bb6 <- c("BB6.0","BB6.90","BB6.180","BB6.270")
  bb7 <- c("BB7.0","BB7.90","BB7.180","BB7.270")
  bb8 <- c("BB8.0","BB8.90","BB8.180","BB8.270")
  tpc <- c(bb1,bb6,bb7,bb8) 
  ppl <- c("P", "PP", "RPP", "SN")    
  pplf <- c("fixed","unpLik","pLik","pLiksp")
  scc <- c("C0","C180","J0","J180","G0","G180",
           "BB1.0","BB1.180","BB6.0","BB6.180",
           "BB7.0","BB7.180","BB8.0","BB8.180")

  #if(PL!="P") stop("Power link function approach not completed yet.")  
  if(!(PL %in% ppl)) stop("Error in parameter PL value. It should be one of: P, PP, RPP or SN.")
  if(BivD %in% tpc ) stop("Archimedean copulas with two parameters not ready yet.") 
  if(!(fitPL %in% pplf)) stop("Error in parameter fitPL value. It should be one of: fixed, unpLik, pLik or pLiksp.")
  if(PL!="P" && BivD %in% tpc ) stop("Models with power link functions and Archimedean copulas with two parameters not available yet.")
  if(!(BivD %in% c(opc,tpc))) stop("Error in parameter BivD value. It should be one of: N,C0,C90,C180,C270,J0,J90,J180,J270,G0,G90,G180,G270,F,T,BB1.0,BB1.90,BB1.180,BB1.270,BB6.0,BB6.90,BB6.180,BB6.270,BB7.0,BB7.90,BB7.180,BB7.270,BB8.0,BB8.90,BB8.180,BB8.270.")
  if(BivD!="N" && H.n==FALSE) stop("Bivariate copula models based on expected information not implemented.")
   
   if(PL!="P"){
      if(eqPL=="both")  {la1 <- valPL[1]; la2 <- valPL[2]; names(la1) <- "xi1.star"; names(la2) <- "xi2.star" }
      if(eqPL=="first") {la1 <- valPL[1];        names(la1) <- "xi1.star"}
      if(eqPL=="second"){la2 <- valPL[2];        names(la2) <- "xi2.star"}   
   } 
   
  if(is.null(start.v)){
       if(BivD %in% c("BB1.0", "BB1.180","BB6.0", "BB6.180","BB7.0", "BB7.180")) {delta.st <-  log(3); names(delta.st) <- "delta.star"} 
       if(BivD %in% c("BB1.90","BB1.270","BB6.90","BB6.270","BB7.90","BB7.270")) {delta.st <- -log(3); names(delta.st) <- "delta.star"}  
       if(BivD %in% bb8){ delta.st <- 0 ; names(delta.st) <- "delta.star"} 
  }
  
   ct <- data.frame( c(opc,tpc),
                     c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,2,7,27,17,37,8,28,18,38,9,29,19,39,10,30,20,40) 
                   )
   nC <- ct[which(ct[,1]==BivD),2]
   
   
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), gamma=gamma, weights=weights, data=data, paraPen=pPen1,control=list(keepData=TRUE)),list(weights=weights)))
  dim.od <- dim(gam1$data)[1]; if(!(is.null(weights))) orw <- weights; mw <- model.weights(gam1$model) 

  if(is.null(mw)){ weights <- rep(1,length(gam1$y))
                        if(dim.od!=length(weights)) weights <- rep(1,dim.od)
                      }

  if(!(is.null(mw))){ weights <- mw
                      if(dim.od!=length(weights)) weights <- orw 
                    }

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 


  if(selection==FALSE){

  startvSS <- NULL
  gam2  <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=weights, data=data, paraPen=pPen2),list(weights=weights)))
  y2 <- gam2$y 
  if(n!=length(y2)) stop("The dimensions of the two binary responses (and respective design matrices) do not match. This may be due to missing values in the dataset. Remove them before fitting the model.")
  mw2 <- model.weights(gam2$model)
  if(is.null(mw2)) weights <- rep(1,length(y2)) else weights <- mw2     

  X2 <- model.matrix(gam2); X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp); n.sel <- NULL
  gp1 <- gam1$nsdf; gp2 <- gam2$nsdf   
  
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE) sp <- c(gam1$sp,gam2$sp)
  if(l.sp1==0 && l.sp2!=0 && fp==FALSE) sp <- c(gam2$sp)
  if(l.sp1!=0 && l.sp2==0 && fp==FALSE) sp <- c(gam1$sp)
  
  

  if(is.null(start.v)){

  if(BivD %in% c("N","F","T")){i.rho <- polychor(y1,y2)
                               i.rho <- atanh( ifelse( abs(i.rho) > 0.90, sign(i.rho)*0.90, i.rho) )
                               if(BivD %in% c("N","T")) names(i.rho) <- "athrho" else names(i.rho) <- "theta.star"
                              }else{ if(BivD %in% scc) i.rho <-  log(3) else i.rho <- -log(3)
                                                       names(i.rho) <- "theta.star"
                                   }


                 if(PL=="P"){
                 if(BivD %in% tpc) start.v <- c(coef(gam1),coef(gam2),i.rho,delta.st) else start.v <- c(coef(gam1),coef(gam2),i.rho) 
                            } else{
                            
                            
                              if(fitPL=="fixed") start.v <- c(coef(gam1),coef(gam2),i.rho) else{  
                            
                 		if(eqPL=="both")   start.v <- c(coef(gam1),coef(gam2),i.rho,la1,la2)         
                 		if(eqPL=="first")  start.v <- c(coef(gam1),coef(gam2),i.rho,la1)           
                 		if(eqPL=="second") start.v <- c(coef(gam1),coef(gam2),i.rho,la2)  
                 		
                 		                                                                }
                 		
                 		
                      		  }


                      }

  cy1 <- NULL                   
  y1.y2 <- y1*y2
  y1.cy2 <- y1*(1-y2)
  cy1.y2 <- (1-y1)*y2
  cy1.cy2 <- (1-y1)*(1-y2)
  
  	  if(PL=="P"){ if(BivD %in% opc ) func.opt <- bprobgHs else func.opt <- bprobgHsBB   }                       
  	  if(PL!="P") func.opt <- bprobgHsPL
    

  }


  if(selection==TRUE){
  inde <- gam1$y > 0
  gam2 <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=weights, data=data, subset=inde, paraPen=pPen2),list(weights=weights,inde=inde)))  
  if(table(inde)[[2]]!=length(gam2$y)) stop("The length of the outcome variable does not match that of the selected observations.")
  X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2) 
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y; n.sel <- length(gam2$y)
  l.sp2 <- length(gam2$sp)
  gp1 <- gam1$nsdf; gp2 <- gam2$nsdf 

  if(is.null(start.v)){ startvSS <- try(startSS(gam1, gam2, formula.eq2, data, gamma, weights, inde, l.sp1, l.sp2, pPen2, fp),silent=TRUE)
                       if(class(startvSS)=="try-error"){ i.rho <- 0.5; 
                                                         names(i.rho) <- "athrho" 
                                                         start.v <- c(coef(gam1),coef(gam2),i.rho)} else start.v <- startvSS$start.v

  if( !(BivD %in% c("N", "T")) ){ if(BivD %in% scc) start.v[length(start.v)] <- log(3) else start.v[length(start.v)] <- -log(3) 
                                  names(start.v)[length(start.v)] <- "theta.star"
                                }

  if(PL=="P"){ if(BivD %in% tpc) start.v <- c(start.v,delta.st) else start.v <- start.v          
             }else{
             
                if(fitPL=="fixed") start.v <- c(coef(gam1),coef(gam2),i.rho) else{ 
                
  		if(eqPL=="both")   start.v <- c(start.v,la1,la2)  
  		if(eqPL=="first")  start.v <- c(start.v,la1)           
  		if(eqPL=="second") start.v <- c(start.v,la2)   
  		                                                                  }
                   }

		      }

  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){if(class(startvSS)=="try-error") sp <- c(gam1$sp,gam2$sp) else sp <- c(gam1$sp,startvSS$gam2.1$sp)}
  if(l.sp1==0 && l.sp2!=0 && fp==FALSE){if(class(startvSS)=="try-error") sp <- c(gam2$sp)         else sp <- c(startvSS$gam2.1$sp)}
  if(l.sp1!=0 && l.sp2==0 && fp==FALSE) sp <- c(gam1$sp)


  cy1 <- (1-y1)
  y1.y2 <- y1*y2; y1.cy2 <- y1*(1-y2)
  cy1.y2 <- cy1.cy2 <- NULL
  
  if(PL=="P"){ if(BivD %in% opc ) func.opt <- bprobgHsSS else func.opt <- bprobgHsSSBB }
  if(PL!="P") func.opt <- bprobgHsSSPL 
  
  
  }


if(PL!="P"){

  if(fitPL=="fixed" || fitPL=="unpLik") sp.xi1 <- sp.xi2 <- 0 else{ sp.xi1 <- spPL[1]; sp.xi2 <- spPL[2]} 

  names(sp.xi1) <- "xi1" 
  names(sp.xi2) <- "xi2"

  if( !(is.null(sp)) ){
    if(eqPL=="both")   sp <- c(sp, sp.xi1, sp.xi2)
    if(eqPL=="first")  sp <- c(sp, sp.xi1)
    if(eqPL=="second") sp <- c(sp, sp.xi2)
                      }else{
    if(eqPL=="both")   sp <- c(sp.xi1, sp.xi2)
    if(eqPL=="first")  sp <- sp.xi1
    if(eqPL=="second") sp <- sp.xi2
                           }
 }


  if( (l.sp1!=0 || l.sp2!=0) ) qu.mag <- S.m(gam1,gam2,l.sp1,l.sp2) 
                   if(PL!="P") qu.mag <- S.mPL(qu.mag,eqPL,start.v) 

  


  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################

  SemiParFit <- SemiParBIVProbit.fit(func.opt=func.opt, start.v=start.v, rinit=rinit, rmax=rmax, BivD=BivD, nC=nC, nu=nu,
                                     PL=PL, eqPL=eqPL, valPL=valPL, fitPL=fitPL, H.n=H.n, 
                                     y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                     X1=X1, X2=X2,  
                                     X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, 
                                     fp=fp,
                                     aut.sp=aut.sp, iterlimsp=iterlimsp, l.sp1=l.sp1, l.sp2=l.sp2, pr.tolsp=pr.tolsp,
                                     weights=weights, fterm=fterm, mterm=mterm, iterlim = 1e+4, gamma=gamma,
                                     control.sp=control.sp) 

  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParBIVProbit.fit.post(SemiParFit=SemiParFit, formula.eq2=formula.eq2, data=data, selection=selection, 
                                            BivD=BivD, nC=nC, nu=nu, PL=PL, eqPL=eqPL, valPL=valPL, fitPL=fitPL, 
                                            y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                            X1=X1, X2=X2,
                                            X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, qu.mag=qu.mag, 
                                            gam1=gam1, gam2=gam2, gp1=gp1, gp2=gp2, 
                                            fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, 
                                            weights=weights)
  SemiParFit <- SemiParFit.p$SemiParFit
  ##########################################################################################################################

L <- list(fit=SemiParFit$fit, gam1=gam1, gam2=gam2, gam2.1=startvSS$gam2.1, 
          coefficients=SemiParFit$fit$argument, weights=weights, 
          sp=SemiParFit.p$sp, iter.sp=SemiParFit$iter.sp, iter.if=SemiParFit$iter.if, iter.inner=SemiParFit$iter.inner,
          rho=SemiParFit.p$rho, theta=SemiParFit.p$theta, delta=SemiParFit.p$delta, KeT=SemiParFit.p$KeT,   
          n=n, n.sel=n.sel, 
          X1=SemiParFit.p$X1, X2=SemiParFit.p$X2, X1.d2=X1.d2, X2.d2=X2.d2, 
          l.sp1=l.sp1, l.sp2=l.sp2, He=SemiParFit.p$He, HeSh=SemiParFit.p$HeSh, Vb=SemiParFit.p$Vb, F=SemiParFit.p$F, 
          fp=fp, aut.sp=aut.sp,
          t.edf=SemiParFit.p$t.edf, edf1=SemiParFit.p$edf1, edf2=SemiParFit.p$edf2, bs.mgfit=SemiParFit$bs.mgfit, conv.sp=SemiParFit$conv.sp, 
          wor.c=SemiParFit$wor.c,
          p11=SemiParFit$fit$p11, p10=SemiParFit$fit$p10, p01=SemiParFit$fit$p01, p00=SemiParFit$fit$p00, p0=SemiParFit$fit$p0,  
          eta1=SemiParFit$fit$eta1, eta2=SemiParFit$fit$eta2, 
          y1=y1,y2=y2,sel=selection, BivD=BivD, nu=nu, 
          PL=PL, eqPL=eqPL, valPL=valPL, fitPL=fitPL, spPL=spPL, xi1=SemiParFit.p$xi1, xi2=SemiParFit.p$xi2, logLik=SemiParFit.p$logLik,
          nC=SemiParFit.p$nC,H.n=H.n,pPen1 = pPen1, pPen2 = pPen2, 
          good=SemiParFit$fit$good,
          y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,qu.mag=qu.mag, gp1=gp1, gp2=gp2, 
          X2s=SemiParFit.p$X2s)

class(L) <- "SemiParBIVProbit"

L

}
