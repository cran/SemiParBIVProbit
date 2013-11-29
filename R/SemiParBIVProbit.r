SemiParBIVProbit <- function(formula.eq1, formula.eq2, data=list(), weights=NULL,  
                             start.v=NULL, BivD="N", nu=3, PL="P", eqPL="both", xi1=1, xi2=1,
                             selection=FALSE, H.n=TRUE, gamma=1, aut.sp=TRUE, fp=FALSE,
                             RE=FALSE, RE.type="N", NGQ=10, K=2, id=NULL, e.npRE=TRUE, pPen1 = NULL, pPen2 = NULL, 
                             rinit=1, rmax=100, fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     iterlimsp=50, pr.tolsp=1e-6, 
                             control.sp=list(maxit=50,tol=1e-6,step.half=25,rank.tol=sqrt(.Machine$double.eps)) ){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  if(!(PL %in% c("P", "PP", "RPP"))) stop("Error in parameter PL value. It should be one of: P, PP or RPP.")
  if(BivD %in% c("BB1.0","BB1.90","BB1.180","BB1.270",
                 "BB6.0","BB6.90","BB6.180","BB6.270",
                 "BB7.0","BB7.90","BB7.180","BB7.270",
                 "BB8.0","BB8.90","BB8.180","BB8.270" ) ) stop("Archimedean copulas with two parameters not ready yet.") 
  if(RE==TRUE && RE.type=="N") stop("Models with bivariate normal random effects not ready yet.")
  if(RE==TRUE && PL!="P") stop("Models with random effects and power link functions not available.")
  if(RE==TRUE && selection=="TRUE") stop("Sample selection models with random effects not available.")
  if(!(RE.type %in% c("N", "NP"))) stop("Error in parameter RE.type value. It should be one of: NP or N.")
  if(!(BivD %in% c("N","C0","C90","C180","C270",
                       "J0","J90","J180","J270",
                       "G0","G90","G180","G270","F","T",
                       "BB1.0","BB1.90","BB1.180","BB1.270",
                       "BB6.0","BB6.90","BB6.180","BB6.270",
                       "BB7.0","BB7.90","BB7.180","BB7.270",
                       "BB8.0","BB8.90","BB8.180","BB8.270" ))) stop("Error in parameter BivD value. It should be one of: N,C0,C90,C180,C270,J0,J90,J180,J270,G0,G90,G180,G270,F,T,
                       BB1.0,BB1.90,BB1.180,BB1.270,BB6.0,BB6.90,BB6.180,BB6.270,BB7.0,BB7.90,BB7.180,BB7.270,BB8.0,BB8.90,BB8.180,BB8.270.")
  if(BivD!="N" && H.n==FALSE) stop("Bivariate copula models based on expected information not implemented.")
  if(BivD %in% c("BB1.0","BB1.90","BB1.180","BB1.270",
                 "BB6.0","BB6.90","BB6.180","BB6.270",
                 "BB7.0","BB7.90","BB7.180","BB7.270",
                 "BB8.0","BB8.90","BB8.180","BB8.270" ) && RE==TRUE) stop("Archimedean copulas with two parameters for the case of models with random effects not implemented.") 
  

   ct <- data.frame( c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","T",
                       "BB1.0","BB1.90","BB1.180","BB1.270",
                       "BB6.0","BB6.90","BB6.180","BB6.270",
                       "BB7.0","BB7.90","BB7.180","BB7.270",
                       "BB8.0","BB8.90","BB8.180","BB8.270"),
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

  X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]; l.sp1 <- length(gam1$sp); y1 <- gam1$y; n <- length(y1) 
  uidf <- sp <- NULL 

  #if(dim(gam1$data)[1]!=n) weights <- weights[-c(gam1$na.action)]              
  #if(!(is.null(gam1$na.action))) weights <- weights[-c(gam1$na.action)] 

  qu.mag <- NULL
  
  if(is.null(start.v)){
    if(BivD %in% c("BB1.0","BB1.180","BB6.0","BB6.180","BB7.0","BB7.180")) {delta.st <- log(3); names(delta.st) <- "delta.star"} 
    if(BivD %in% c("BB1.90","BB1.270","BB6.90","BB6.270","BB7.90","BB7.270")) {delta.st <- -log(3); names(delta.st) <- "delta.star"}  
    if(BivD %in% c("BB8.0","BB8.180","BB8.90","BB8.270")){ delta.st <- 0 ; names(delta.st) <- "delta.star"} 
    #if(BivD %in% c("BB8.0","BB8.180","BB8.90","BB8.270")){ delta.st <- 0 ; names(delta.st) <- "delta.star"} 
  }

  if(selection==FALSE){

  startvSS <- NULL
  gam2  <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=weights, data=data, paraPen=pPen2),list(weights=weights)))
  y2 <- gam2$y 
  if(n!=length(y2)) stop("The dimensions of the two binary responses (and respective design matrices) do not match. This may be due to missing values in the dataset. Remove them before fitting the model.")
  mw2 <- model.weights(gam2$model)
  if(is.null(mw2)) weights <- rep(1,length(y2)) else weights <- mw2     

  if(RE==TRUE && RE.type=="N") weights <- rep(weights,each=NGQ^2)
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
                              }
  else{ if(BivD %in% c("C0", "C180","J0", "J180","G0", "G180",
                       "BB1.0","BB1.180","BB6.0","BB6.180",
                       "BB7.0","BB7.180","BB8.0","BB8.180")) i.rho <-  log(3) else i.rho <- -log(3)
                                                                names(i.rho) <- "theta.star"
      }


  if(RE==FALSE){ #if(PL=="P"){

                     if(BivD %in% c("BB1.0","BB1.90","BB1.180","BB1.270",
                        	    "BB6.0","BB6.90","BB6.180","BB6.270",
                         	    "BB7.0","BB7.90","BB7.180","BB7.270",
                                    "BB8.0","BB8.90","BB8.180","BB8.270" )) start.v <- c(coef(gam1),coef(gam2),i.rho,delta.st) else start.v <- c(coef(gam1),coef(gam2),i.rho) 
                              
                            #} 
  
                 
  
                 #else {

                 #if(eqPL=="both"){
                 #la1 <- la2 <- 0; names(la1) <- "lambda1.star"; names(la2) <- "lambda2.star"
                 #start.v <- c(coef(gam1),coef(gam2),i.rho,la1,la2)
                 #}
                 #if(eqPL=="first"){
                 #la1 <- 0; names(la1) <- "lambda1.star"
                 #start.v <- c(coef(gam1),coef(gam2),i.rho,la1)
                 #}
                 #if(eqPL=="second"){
                 #la2 <- 0; names(la2) <- "lambda2.star"
                 #start.v <- c(coef(gam1),coef(gam2),i.rho,la2)
                 #}                   
                 
                 
                 #}
 
                 
               }

                      }

  cy1 <- NULL
  # if(RE==FALSE || (RE==TRUE && RE.type=="NP")) 
  y1.y2 <- y1*y2; y1.cy2 <- y1*(1-y2); cy1.y2 <- (1-y1)*y2; cy1.cy2 <- (1-y1)*(1-y2)
  
  if(RE==TRUE && RE.type=="N") func.opt <- bprobgHs.NRE 
  else{ 
  
  if(PL=="P"){ if(BivD %in% c("N","C0","C90","C180","C270",
                              "J0","J90","J180","J270",
                              "G0","G90","G180","G270","F","T") ) func.opt <- bprobgHs else func.opt <- bprobgHsBB   }
                              
  if(PL!="P" && BivD %in% c("N","C0","C90","C180","C270",
                              "J0","J90","J180","J270",
                              "G0","G90","G180","G270","F","T") ) func.opt <- bprobgHsPL
  
  
  
  }  

       

    if(RE==TRUE){
     
     if(dim(X1)[1]!=length(id)) stop("The length of your id does not match the number of rows of your data frame.")

      if(RE.type=="NP"){ 

           	if( (l.sp1!=0 || l.sp2!=0) && fp==FALSE) qu.mag <- S.m(gam1,gam2,l.sp1,l.sp2,K,RE=FALSE,RE.type)
      
           	if(is.null(start.v)){ 
           		start.v <- c(coef(gam1),coef(gam2),i.rho)
           		fit     <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, BivD=BivD, nC=nC, nu=nu, xi1=xi1, xi2=xi2, PL=PL, eqPL=eqPL, H.n=H.n, 
                                         y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                         X1=X1, X2=X2,
                		         X1.d2=X1.d2, X2.d2=X2.d2, sp=sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, fp=fp, 
                		         l.sp1=l.sp1, l.sp2=l.sp2, blather=TRUE, 
                        		 fterm=fterm, mterm=mterm, iterlim = 1e+4, weights=weights,
                         		K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses, NGQ=NGQ, dat1all=dat1all, dat2all=dat2all, W=W) 
           		start.v <- fit$argument }
                            
                        X1.d2 <- X1.d2-1; X2.d2 <- X2.d2-1 # THIS LINE MUST BE HERE     
        		start.v <- start.v[-c(1,X1.d2+2)]
         		start.v <- c(sqrt(2)*gauss.quad(K,"hermite")$nodes,start.v[1:X1.d2],
                      		     sqrt(2)*gauss.quad(K,"hermite")$nodes,start.v[(X1.d2+1):(X2.d2+X1.d2)],
                      	     	sign(start.v[X2.d2+X1.d2+1])*i.rho)
                	for(i in 1:K) {names(start.v)[i] <- paste("re",i,sep="")
                        	       names(start.v)[(X1.d2+K+i)] <- paste("re",i,sep="") }	              
                              
        	uid <- unique(id) 
      		N <- length(uid) 
      		uidf <- array(0,N) 
      		for (j in 1:N) uidf[j] <- sum(id==unique(id)[j]) 
      		cuid <- c(0,cumsum(uidf))

      		masses <- rep(1/K,K) 
      		func.opt <- bprobgHs.NPRE

                xx1 <- as.matrix(X1[,-1]); xx2 <- as.matrix(X2[,-1])
          	X1 <- function(u) cbind(matrix(replace(rep(0,K),u,1),ncol=K,nrow=n,byrow=TRUE),xx1)
          	X2 <- function(u) cbind(matrix(replace(rep(0,K),u,1),ncol=K,nrow=n,byrow=TRUE),xx2)

	        	if(e.npRE==TRUE){
		     		start.ve <- extraiterNP(paramNP=start.v, BivD, nC, nu, xi1, xi2, PL, eqPL,
                                                        y1.y2, y1.cy2, cy1.y2, cy1.cy2, cy1,
                                                        X1, X2, X1.d2, X2.d2, pPen1, pPen2, sp, qu.mag, gp1, gp2,
		     		                        fp, l.sp1, l.sp2, weights, K, n, N, cuid, uidf, masses, pr.tolsp, NGQ, 
		     		                        dat1all, dat2all, W)                
             	     		start.v <- start.ve$paramNP; masses <- start.ve$masses} 
                     }


      if(RE.type=="N"){ 

                uidf <- table(id)
                
                cuid <- c(0,cumsum(uidf))
                
                N    <- length(uidf)

                log.sigma1 <- log(1.1); names(log.sigma1) <- "log.sigma1"
                athrho.u  <- atanh(0.5); names(athrho.u) <- "athrho.u"
                log.sigma2 <- log(1.1); names(log.sigma2) <- "log.sigma2"

           	if(is.null(start.v)) start.v <- c(coef(gam1),log.sigma1,coef(gam2),athrho.u,log.sigma2,i.rho)

                Z1 <- Z2 <- sqrt(2)*gauss.quad(NGQ,"hermite")$nodes
                Z <- expand.grid(Z1,Z2)
  
                W1 <- W2 <- gauss.quad(NGQ,"hermite")$weights/sqrt(pi)
                W <- apply(expand.grid(W1,W2),1,prod)

  
                Zexp <- matrix(rep(unlist(t(Z)),n),ncol=2,byrow=TRUE)
  
                dat1exp <- matrix(rep(X1,each=NGQ^2),ncol=X1.d2)
                dat2exp <- matrix(rep(X2,each=NGQ^2),ncol=X2.d2)
                dat1all <- cbind(dat1exp,Zexp[,1])
                dat2all <- cbind(dat2exp,Zexp)

                y1e <- rep(y1,each=NGQ^2)
                y2e <- rep(y2,each=NGQ^2)
                          
                y1.y2 <- y1e*y2e; y1.cy2 <- y1e*(1-y2e); cy1.y2 <- (1-y1e)*y2e; cy1.cy2 <- (1-y1e)*(1-y2e)

    
                       }

     }



  }


  if(selection==TRUE){
  if(RE==TRUE) stop("Routine for bivariate probit sample selection modelling with random effects not implemented yet.")
  inde <- gam1$y > 0
  gam2 <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=weights, data=data, subset=inde, paraPen=pPen2),list(weights=weights,inde=inde)))  
  if(table(inde)[[2]]!=length(gam2$y)) stop("The length of the outcome variable does not match that of the selected observations.")
  X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2) 
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y; n.sel <- length(gam2$y)
  l.sp2 <- length(gam2$sp)

  if(is.null(start.v)){startvSS <- startSS(gam1, gam2, formula.eq2, data, gamma, weights, inde, l.sp1, l.sp2, pPen2, fp)
                       start.v <- startvSS$start.v

  


  if( !(BivD %in% c("N", "T")) ){ if(BivD %in% c("C0","C180","J0","J180","G0","G180",
                                                 "BB1.0","BB1.180","BB6.0","BB6.180",
                                                 "BB7.0","BB7.180","BB8.0","BB8.180")) start.v[length(start.v)] <- log(3) else start.v[length(start.v)] <- -log(3) 
                                  names(start.v)[length(start.v)] <- "theta.star"
                                }

  #if(PL=="P"){ 
               if(BivD %in% c("BB1.0","BB1.90","BB1.180","BB1.270",
                              "BB6.0","BB6.90","BB6.180","BB6.270",
                              "BB7.0","BB7.90","BB7.180","BB7.270",
                              "BB8.0","BB8.90","BB8.180","BB8.270" )) start.v <- c(start.v,delta.st) else start.v <- start.v 
                                
                            #}else{
  
  #if(eqPL=="both"){
  #la1 <- la2 <- 0; names(la1) <- "lambda1.star"; names(la2) <- "lambda2.star";  
  #start.v <- c(start.v,la1,la2)
  #}
  #  if(eqPL=="first"){
  #  la1 <- 0; names(la1) <- "lambda1.star"  
  #  start.v <- c(start.v,la1)
  #}
  #    if(eqPL=="second"){
  #la2 <- 0; names(la2) <- "lambda2.star";  
  #start.v <- c(start.v,la2)
  #}
  
  
  #}


        gp1 <- gam1$nsdf; gp2 <- gam2$nsdf 
  	if(l.sp1!=0 && l.sp2!=0 && fp==FALSE) sp <- c(gam1$sp,startvSS$gam2.1$sp)
  	if(l.sp1==0 && l.sp2!=0 && fp==FALSE) sp <- c(startvSS$gam2.1$sp)
        if(l.sp1!=0 && l.sp2==0 && fp==FALSE) sp <- c(gam1$sp)

			}

  cy1 <- (1-y1)
  y1.y2 <- y1*y2; y1.cy2 <- y1*(1-y2)
  cy1.y2 <- cy1.cy2 <- NULL
  
  if(PL=="P"){ if(BivD %in% c("N","C0","C90","C180","C270",
                                "J0","J90","J180","J270",
                                "G0","G90","G180","G270","F","T") ) func.opt <- bprobgHsSS else func.opt <- bprobgHsSSBB }
  
  
  if(PL!="P" && BivD %in% c("N","C0","C90","C180","C270",
                                "J0","J90","J180","J270",
                                "G0","G90","G180","G270","F","T") ) func.opt <- bprobgHsSSPL 
  
  }



  if( (l.sp1!=0 || l.sp2!=0) && fp==FALSE) qu.mag <- S.m(gam1,gam2,l.sp1,l.sp2,K,RE,RE.type) 


  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################

  SemiParFit <- SemiParBIVProbit.fit(func.opt=func.opt, start.v=start.v, rinit=rinit, rmax=rmax, BivD=BivD, nC=nC, nu=nu,
                                     xi1=xi1, xi2=xi2, PL=PL, eqPL=eqPL, H.n=H.n, 
                                     y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                     X1=X1, X2=X2,  
                                     X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, sp=sp, qu.mag=qu.mag, gp1=gp1, gp2=gp2, 
                                     fp=fp, RE=RE, RE.type=RE.type,
                                     aut.sp=aut.sp, iterlimsp=iterlimsp, l.sp1=l.sp1, l.sp2=l.sp2, pr.tolsp=pr.tolsp,
                                     weights=weights, fterm=fterm, mterm=mterm, iterlim = 1e+4, e.npRE=e.npRE, gamma=gamma,
                                     K=K, n=n, N=N, cuid=cuid, uidf=uidf, masses=masses, control.sp=control.sp,NGQ=NGQ,
                                     dat1all=dat1all, dat2all=dat2all, W=W) 

  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParBIVProbit.fit.post(SemiParFit=SemiParFit, formula.eq2=formula.eq2, data=data, selection=selection, 
                                            RE=RE, RE.type=RE.type, BivD=BivD, nC=nC, nu=nu, xi1=xi1, xi2=xi2, PL=PL, eqPL=eqPL, 
                                            y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,
                                            X1=X1, X2=X2,
                                            X1.d2=X1.d2, X2.d2=X2.d2, pPen1=pPen1, pPen2=pPen2, qu.mag=qu.mag, gam1=gam1, gam2=gam2, gp1=gp1, gp2=gp2, 
                                            fp=fp, l.sp1=l.sp1, l.sp2=l.sp2, 
                                            weights=weights, K=K, n=n, N=N, cuid=cuid, uidf=uidf)

  ##########################################################################################################################

  SemiParFit <- SemiParFit.p$SemiParFit # needed for NP re stuff

                               

L <- list(fit=SemiParFit$fit, gam1=gam1, gam2=gam2, gam2.1=startvSS$gam2.1, 
          coefficients=SemiParFit$fit$argument, weights=weights, 
          sp=SemiParFit.p$sp, iter.sp=SemiParFit$iter.sp, iter.if=SemiParFit$iter.if, iter.inner=SemiParFit$iter.inner,
          rho=SemiParFit.p$rho, theta=SemiParFit.p$theta, delta=SemiParFit.p$delta, KeT=SemiParFit.p$KeT,   
          sigma1=SemiParFit.p$sigma1, rho.u=SemiParFit.p$rho.u, sigma2=SemiParFit.p$sigma2, n=n, n.sel=n.sel, 
          X1=SemiParFit.p$X1, X2=SemiParFit.p$X2, X1.d2=X1.d2, X2.d2=X2.d2, 
          l.sp1=l.sp1, l.sp2=l.sp2, He=SemiParFit.p$He, HeSh=SemiParFit.p$HeSh, Vb=SemiParFit.p$Vb, F=SemiParFit.p$F, 
          fp=fp, aut.sp=aut.sp,
          t.edf=SemiParFit.p$t.edf, edf1=SemiParFit.p$edf1, edf2=SemiParFit.p$edf2, bs.mgfit=SemiParFit$bs.mgfit, conv.sp=SemiParFit$conv.sp, 
          wor.c=SemiParFit$wor.c,
          p11=SemiParFit$fit$p11, p10=SemiParFit$fit$p10, p01=SemiParFit$fit$p01, p00=SemiParFit$fit$p00, p0=SemiParFit$fit$p0,  
          eta1=SemiParFit$fit$eta1, eta2=SemiParFit$fit$eta2, 
          y1=y1,y2=y2,sel=selection, K=SemiParFit.p$K, masses=SemiParFit$fit$masses, RE=RE, RE.type=RE.type, BivD=BivD, nu=nu, 
          PL=PL, eqPL=eqPL, xi1=SemiParFit.p$xi1, xi2=SemiParFit.p$xi2, logLik=SemiParFit.p$logLik,
          eb.u1=SemiParFit.p$NP.qu.int$eb.u1, eb.u2=SemiParFit.p$NP.qu.int$eb.u2, 
          Eb.u1=SemiParFit.p$NP.qu.int$Eb.u1, Eb.u2=SemiParFit.p$NP.qu.int$Eb.u2, 
          id=id, uidf=uidf, T.sv=SemiParFit.p$NP.qu.int$T.sv,
          eta1S=SemiParFit.p$eta1S,eta2S=SemiParFit.p$eta2S,ass.pS=SemiParFit.p$ass.pS,
          #lambda1S=SemiParFit.p$lambda1s, lambda2S=SemiParFit.p$lambda2s, 
          deltaS=SemiParFit.p$deltaS,
          nC=SemiParFit.p$nC,H.n=H.n,pPen1 = pPen1, pPen2 = pPen2, 
          good=SemiParFit$fit$good,
          y1.y2=y1.y2, y1.cy2=y1.cy2, cy1.y2=cy1.y2, cy1.cy2=cy1.cy2, cy1=cy1,qu.mag=qu.mag, gp1=gp1, gp2=gp2)

class(L) <- "SemiParBIVProbit"

L

}
