SemiParBIVProbit <- function(formula, data = list(), weights = NULL, subset = NULL, start.v = NULL, 
                             Model = "B", method = "trwlF", BivD = "N", nu = 3, fp = FALSE,
                             PL = "P", eqPL = "both", valPL = c(0,0), fitPL = "pLiksp", spPL = c(0.01,0.01),
                             hess = TRUE, infl.fac = 1, pPen1 = NULL, pPen2 = NULL, 
                             rinit = 1, rmax = 100, iterlimsp = 50, pr.tolsp = 1e-6,
                             gc.l = FALSE, awlm = FALSE, parscale, extra.regI = "t"){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  sp <- qu.mag <- sp.xi1 <- sp.xi2 <- startvSS <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- NULL  
  end <- 0
  selection <- FALSE; if(Model=="BSS") selection <- TRUE
  
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","T")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180")
  ppl  <- c("P", "PP", "RPP", "SN")    
  pplf <- c("fixed","unpLik","pLik","pLiksp")
  mb   <- c("B", "BSS", "BPO")

  #if(Model == "BPO") stop("Check next release for final tested version of this model.")

  if(!(Model %in% mb)) stop("Error in parameter Model value. It should be one of: B, BSS or BPO.")
  if(!(PL %in% ppl)) stop("Error in parameter PL value. It should be one of: P, PP, RPP or SN.")
  if(!(fitPL %in% pplf)) stop("Error in parameter fitPL value. It should be one of: fixed, unpLik, pLik or pLiksp.")
  if(!(BivD %in% opc)) stop("Error in parameter BivD value. It should be one of: N, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270, F, T")
  if(PL!="P" && Model %in% c("BSS","BPO")) stop("Sample selection models with asymmetric links not implemented.")
  if(!(method %in% c("trwl","trwlF"))) stop("Error in parameter method value. It should be one of: trwl or trwlF.")
  if(!(extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of: t, pC or sED.")
  
  ig <- interpret.gam(formula)
  mf <- match.call(expand.dots = FALSE)

  pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names,ig[[2]]$response))
  fake.formula <- paste(ig[[1]]$response, "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(ig$fake.formula)
  mf$formula <- fake.formula 
  mf$start.v <- mf$Model <- mf$BivD <- mf$nu <- mf$fp <- mf$PL <- mf$eqPL <- mf$valPL <- mf$fitPL <- mf$spPL <- mf$hess <- mf$infl.fac <- mf$pPen1 <- mf$pPen2 <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$pr.tolsp <- mf$gc.l <- mf$awlm <- mf$method <- mf$parscale <- mf$extra.regI <- NULL                           
  mf$drop.unused.levels <- TRUE 
  if(Model=="BSS") mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()
        
  if(Model=="BSS"){ 
     indS <- as.logical(data[,ig[[1]]$response])==FALSE 
     indS <- ifelse( is.na(indS), FALSE, indS) 
     data[indS, ig[[2]]$response] <- ifelse( is.na(data[indS, ig[[2]]$response]), 0, data[indS, ig[[2]]$response]) 
     data <- na.omit(data)
     }
  
  if(is.null(weights)) weights <- rep(1,dim(data)[1]) else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  
  
  if(Model=="B"){
  if(ig[[1]]$response %in% ig[[2]]$pred.names ) end <- 1
  if(ig[[2]]$response %in% ig[[1]]$pred.names ) end <- 2
  }

  ct <- data.frame( c(opc),
                    c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,2) 
                     )
  nC <- ct[which(ct[,1]==BivD),2]
  
  
   if(PL!="P" && selection == FALSE){
      if(eqPL=="both")  {la1 <- valPL[1]; la2 <- valPL[2]; names(la1) <- "xi1.star"; names(la2) <- "xi2.star"}
      if(eqPL=="first") {la1 <- valPL[1];                  names(la1) <- "xi1.star"}
      if(eqPL=="second"){la2 <- valPL[2];                                            names(la2) <- "xi2.star"}   
   } 
   
   

  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), gamma=infl.fac, weights=weights, 
                              data=data, paraPen=pPen1),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 



  if(Model=="B" || Model=="BPO"){
  
  gam2  <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=infl.fac, weights=weights, 
                           data=data, paraPen=pPen2),list(weights=weights))) # check at later stage the need of eval and substitute

  X2 <- model.matrix(gam2)
  X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp)
  y2 <- gam2$y 
  
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE) sp <- c(gam1$sp,gam2$sp)
  if(l.sp1==0 && l.sp2!=0 && fp==FALSE) sp <- c(gam2$sp)
  if(l.sp1!=0 && l.sp2==0 && fp==FALSE) sp <- c(gam1$sp)
  
  

  if(is.null(start.v)){

  if(BivD %in% c("N","F","T")) i.rho <- 0.1 else{ if(BivD %in% scc) i.rho <-  log(3) else i.rho <- -log(3)}              
  if(BivD %in% c("N","T")) names(i.rho) <- "athrho" else names(i.rho) <- "theta.star"               


  if(PL=="P") start.v <- c(coef(gam1), coef(gam2), i.rho) else{
                            
                            
                            if(fitPL=="fixed") start.v <- c(coef(gam1),coef(gam2),i.rho) else{  
                            
                 		if(eqPL=="both")   start.v <- c(coef(gam1),coef(gam2),i.rho,la1,la2)         
                 		if(eqPL=="first")  start.v <- c(coef(gam1),coef(gam2),i.rho,la1)           
                 		if(eqPL=="second") start.v <- c(coef(gam1),coef(gam2),i.rho,la2)  
                 		
                 		                                                              }

                      		                               }
                      }
   
   
   
  if(Model=="B"){ 
   
  y1.y2 <- y1*y2
  y1.cy2 <- y1*(1-y2)
  cy1.y2 <- (1-y1)*y2
  cy1.cy2 <- (1-y1)*(1-y2)

  }

  
  if(Model=="BPO" ){ 

      cy <- 1 - y1
  
  }
  
  
  

  } # end big if





  if(Model=="BSS"){

  inde <- y1 > 0

  gam2 <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=infl.fac, weights=weights, 
                              data=data, subset=inde, paraPen=pPen2),list(weights=weights,inde=inde)))  
                              
  X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2) 
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y   # ; n.sel <- length(gam2$y)
  l.sp2 <- length(gam2$sp)


  if(is.null(start.v)){ startvSS <- try(startSS(gam1, gam2, formula.eq2, data, infl.fac, weights, inde, l.sp1, l.sp2, pPen2, fp),silent=TRUE)
                       if(class(startvSS)=="try-error"){ i.rho <- 0.5; names(i.rho) <- "athrho" 
                                                         start.v <- c(coef(gam1),coef(gam2),i.rho)} else start.v <- startvSS$start.v

  if( !(BivD %in% c("N", "T")) ){ if(BivD %in% scc) start.v[length(start.v)] <- log(3) else start.v[length(start.v)] <- -log(3) 
                                  names(start.v)[length(start.v)] <- "theta.star"
                                }

		      }
		      
		      
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){if(class(startvSS)=="try-error") sp <- c(gam1$sp,gam2$sp) else sp <- c(gam1$sp,startvSS$gam2.1$sp)}
  if(l.sp1==0 && l.sp2!=0 && fp==FALSE){if(class(startvSS)=="try-error") sp <- c(gam2$sp)         else sp <- c(startvSS$gam2.1$sp)}
  if(l.sp1!=0 && l.sp2==0 && fp==FALSE) sp <- c(gam1$sp)


  cy1 <- (1-y1)
  y1.y2 <- y1*y2
  y1.cy2 <- y1*(1-y2)
  
  }


  gp1 <- gam1$nsdf 
  gp2 <- gam2$nsdf   

  if( (l.sp1!=0 || l.sp2!=0) )  qu.mag <- S.m(gam1,gam2,l.sp1,l.sp2) 


if(PL!="P" && selection == FALSE){

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
                           
 qu.mag <- S.mPL(qu.mag,eqPL,start.v)                           
                           
}


if(missing(parscale)) parscale <- rep(1, length(start.v))



  respvec <- list(y1 = y1,
                  y2 = y2,
                  y1.y2 = y1.y2, 
                  y1.cy2 = y1.cy2, 
                  cy1.y2 = cy1.y2, 
                  cy1.cy2 = cy1.cy2, 
                  cy1 = cy1,
                  cy = cy)
  
  VC <- list(X1 = X1, 
             X2 = X2, 
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             gp1 = gp1, 
             gp2 = gp2,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             infl.fac = infl.fac,
             weights = weights,
             fp = fp,
             hess = hess,
             pPen1 = pPen1,
             pPen2 = pPen2,
             Model = Model,
             PL = PL, end = end,
             BivD = BivD, nu = nu,
             nC = nC, gc.l = gc.l, awlm = awlm, n = n, extra.regI = extra.regI,
             parscale = parscale) # original n only needed in SemiParBIVProbit.fit
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################

  if(Model=="B"){ 
  if(PL=="P") func.opt <- bprobgHs   
  if(PL!="P") func.opt <- bprobgHsPL
                }
  if(Model=="BPO") func.opt <- bprobgHsPO 
  if(Model=="BSS") func.opt <- bprobgHsSS 

  if(method == "trwl")  fit.func <- SemiParBIVProbit.fit
  if(method == "trwlF") fit.func <- SemiParBIVProbit.fit1

  SemiParFit <- fit.func(func.opt = func.opt, start.v = start.v, 
                                     rinit = rinit, rmax = rmax, iterlim = 1e+4, iterlimsp = iterlimsp, pr.tolsp = pr.tolsp,
                                     PL = PL, eqPL = eqPL, valPL = valPL, fitPL = fitPL, 
                                     respvec = respvec, VC = VC, 
                                     sp = sp, qu.mag = qu.mag) 
  
  n <- sum(as.numeric(SemiParFit$fit$good))  
  if(Model=="BSS") n.sel <- sum(as.numeric(inde[SemiParFit$fit$good]))
  
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParBIVProbit.fit.post(SemiParFit = SemiParFit, formula.eq2 = formula.eq2, data = data, 
                                            Model = Model, VC = VC,  
                                            PL = PL, eqPL = eqPL, valPL = valPL, fitPL = fitPL, 
                                            qu.mag = qu.mag, gam1 = gam1, gam2 = gam2)
                                            
  SemiParFit <- SemiParFit.p$SemiParFit # useful for SS models, eta2 calculatons etc.
  ##########################################################################################################################

rm(data)
if(gc.l == TRUE) gc()


L <- list(fit = SemiParFit$fit, 
          gam1 = gam1, gam2 = gam2, gam2.1 = startvSS$gam2.1, 
          coefficients = SemiParFit$fit$argument, 
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, l.sp1 = l.sp1, l.sp2 = l.sp2, fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          rho = SemiParFit.p$rho, theta = SemiParFit.p$theta, OR = SemiParFit.p$OR, GM = SemiParFit.p$GM, #####  KeT = SemiParFit.p$KeT,   
          n = n, n.sel = n.sel, 
          X1 = X1, X2 = X2, X1.d2 = X1.d2, X2.d2 = X2.d2, 
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, F = SemiParFit.p$F, 
          t.edf = SemiParFit.p$t.edf, edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, 
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,
          p11 = SemiParFit$fit$p11, p10 = SemiParFit$fit$p10, p01 = SemiParFit$fit$p01, p00 = SemiParFit$fit$p00, p0 = SemiParFit$fit$p0,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
          y1 = y1, y2 = y2, 
          sel = selection, 
          BivD = BivD, nu = nu, 
          PL = PL, eqPL = eqPL, valPL = valPL, fitPL = fitPL, spPL = spPL, xi1 = SemiParFit.p$xi1, xi2 = SemiParFit.p$xi2, 
          logLik = SemiParFit.p$logLik,
          nC = nC, hess = hess, pPen1 = pPen1, pPen2 = pPen2, 
          good = SemiParFit$fit$good,
          respvec = respvec,
          qu.mag = qu.mag, 
          gp1 = gp1, gp2 = gp2, 
          X2s = SemiParFit.p$X2s,
          VC = VC, Model = Model, ig = ig, method = method, magpp = SemiParFit$magpp)

class(L) <- "SemiParBIVProbit"

L

}
