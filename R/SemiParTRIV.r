SemiParTRIV <- function(formula, data = list(), weights = NULL, subset = NULL,
                             Model = "T", penCor = "unpen", sp.penCor = 3, approx = FALSE, 
                             infl.fac = 1, gamma = 1, w.alasso = NULL, 
                             rinit = 1, rmax = 100, 
                             iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t"){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  i.rho <- sp <- qu.mag <- qu.mag1 <- qu.mag2 <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss <- inde <- spgamlss <- n.sel1 <- n.sel2 <- NULL  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- 0
  ngc <- 2
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- NULL
  
  fp <- FALSE
  
  X2s <- X3s <- NULL
    
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- NULL  
  sp7 <- gp7 <- gam7 <- X7 <- NULL  
  spCor <- NULL

  y1.y2.y3    <- NULL
  y1.y2.cy3   <- NULL
  cy1.y2.y3   <- NULL 
  cy1.y2.cy3  <- NULL
  cy1.cy2.cy3 <- NULL
  cy1.cy2.y3  <- NULL
  y1.cy2.cy3  <- NULL
  y1.cy2.y3   <- NULL

  cy1         <- NULL
  y1.cy2      <- NULL
  y1.y2.cy3   <- NULL
  y1.y2.y3    <- NULL


  if(penCor != "unpen"){ spCor <- sp.penCor; names(spCor) <- "spCor"} 
  

  
  opc  <- NULL 
  scc  <- NULL 
  sccn <- NULL 
  mb   <- NULL 
  m2   <- m2d <- NULL 
  m3   <- m3d <- NULL 
  
  bl   <- NULL 
  
  hess <- TRUE
    
    
    
  mb   <- c("T", "TSS")
  if(!(Model %in% mb)) stop("Error in parameter Model value. It should be one of: T, TSS.")
    
    
  if(!(penCor %in% c("unpen", "ridge", "lasso", "alasso"))) stop("Error in parameter penCor value. It should be one of:\nunpen, ridge, lasso, alasso.")

  if(is.null(w.alasso) && penCor == "alasso") stop("You must provide a vector of three weights when using alasso.")

  if(length(w.alasso)!=3 && penCor == "alasso") stop("You must supply a vector made up of three weights.")


  l.flist <- length(formula)
  if(!is.list(formula)) stop("You must specify a list of equations.")
  
  if(!(extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of:\nt, pC or sED.")
  
  if(l.flist > 3) stop("You are only allowed to specify three equations.") 
  
  
  if(Model == "TSS") stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")

 #######################################################################################  

 
  mf  <- match.call(expand.dots = FALSE)
  ig <- interpret.gam(formula)
  
  if( l.flist == 3 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- all.vars(as.formula(formula[[3]]))[1]
    v3 <- c(v3, ig[[3]]$pred.names)    
    
    pred.n <- union(v1,c(v2,v3))
                    } 
  
  
                    
                    
  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$infl.fac <- mf$rinit <- mf$approx <- mf$gamma <- mf$w.alasso <- mf$rmax <- mf$Model <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$penCor <- mf$sp.penCor <- NULL                           
  mf$drop.unused.levels <- TRUE 
  if(Model=="TSS") mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  
  
  if(Model=="TSS"){  
  
     data[is.na(data[, v1[1]]), v1[1]] <- 0
     indS1 <- data[, v1[1]] 
     indS2 <- data[, v2[1]]     
     indS1[is.na(indS1)] <- 0   
     indS2[is.na(indS2)] <- 0       
     indS1 <- as.logical(indS1)
     indS2 <- as.logical(indS2) 
     
     data[indS1 == FALSE, v2[1]] <- 0
     data[indS1 == FALSE, v3[1]] <- 0
     data[indS2 == FALSE, v3[1]] <- 0     
     
     data <- na.omit(data)
     
                   }
  

  if(is.null(weights)) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  formula.eq3 <- formula[[3]]  
  
  
 # if(Model=="B"){
 #   if(v1[1] %in% v2[-1]) end <- 1
 #   if(v2[1] %in% v1[-1]) end <- 2
 #               }


  ct <- cta <- nC  <- nCa <- NULL  
  
  
 ##############################################################  
 # Equations 1, 2 and 3
 ##############################################################  
 
  
  
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link = "probit"), gamma=infl.fac, weights=weights, data=data),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 
  if(l.sp1 != 0) sp1 <- gam1$sp
  
  inde1 <- inde2 <- inde2.1 <- rep(TRUE, n) # useful for double ss
  
  if(Model == "TSS") inde1 <- inde2 <- as.logical(y1)
  
  

  gam2 <- eval(substitute(gam(formula.eq2, binomial(link = "probit"), gamma=infl.fac, weights=weights, data=data, subset=inde1),list(weights=weights,inde1=inde1))) 

  if(Model == "TSS"){
  
  X2s <- try(predict.gam(gam2, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
  if(class(X2s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?SemiParTRIV for more information.")    
  
  }
  
  X2 <- model.matrix(gam2)
  X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp)
  y2 <- gam2$y 
  if(l.sp2 != 0) sp2 <- gam2$sp
  

 
  
  if(Model == "TSS"){ inde2[inde1] <- as.logical(gam2$y); inde2.1 <- inde2[inde1]}
  

  gam3 <- eval(substitute(gam(formula.eq3, binomial(link = "probit"), gamma=infl.fac, weights=weights, data=data, subset=inde2),list(weights=weights,inde2=inde2))) 


  if(Model == "TSS"){
  
  X3s <- try(predict.gam(gam3, newdata = data[,-dim(data)[2]], type = "lpmatrix"), silent = TRUE)
  if(class(X3s)=="try-error") stop("Check that the numbers of factor variables' levels\nin the selected sample are the same as those in the complete dataset.\nRead the Details section in ?SemiParTRIV for more information.")    
    
  }  
    
  X3 <- model.matrix(gam3)
  X3.d2 <- dim(X3)[2]
  l.sp3 <- length(gam3$sp)
  y3 <- gam3$y 
  if(l.sp3 != 0) sp3 <- gam3$sp  


 ##############################################################

if(Model == "T"){

  y1.y2.y3    <- y1*y2*y3 
  y1.y2.cy3   <- y1*y2*(1-y3)
  cy1.y2.y3   <- (1-y1)*y2*y3 
  cy1.y2.cy3  <- (1-y1)*y2*(1-y3) 
  cy1.cy2.cy3 <- (1-y1)*(1-y2)*(1-y3)
  cy1.cy2.y3  <- (1-y1)*(1-y2)*y3
  y1.cy2.cy3  <- y1*(1-y2)*(1-y3)
  y1.cy2.y3   <- y1*(1-y2)*y3
  
} 

if(Model == "TSS"){
 
  cy1       <- (1-y1)
  y1.cy2    <- y1[inde1]*(1-y2)
  y1.y2.cy3 <- y1[inde2]*y2[inde2.1]*(1-y3) 
  y1.y2.y3  <- y1[inde2]*y2[inde2.1]*y3  
  n.sel1 <- table(inde1)[2]
  n.sel2 <- table(inde2)[2]

}


  gp1 <- gam1$nsdf 
  gp2 <- gam2$nsdf
  gp3 <- gam3$nsdf

##############################################################
# Starting values
##############################################################

if(Model == "T"){

res1 <- residuals(gam1)
res2 <- residuals(gam2)
res3 <- residuals(gam3)

theta12 <- atanh(cor(res1, res2))
theta13 <- atanh(cor(res1, res3))
theta23 <- atanh(cor(res2, res3))

}


if(Model == "TSS") theta12 <- theta13 <- theta23 <- atanh(0.01)
  
names(theta12) <- "theta12.st"
names(theta13) <- "theta13.st"
names(theta23) <- "theta23.st"
  
start.v <- c(coef(gam1), coef(gam2), coef(gam3), theta12, theta13, theta23)

##############################################################

if(Model == "T")   func.opt <- triprobgHs
if(Model == "TSS") func.opt <- triprobgHsSS
  
l.gam1 <- length(coef(gam1))
l.gam2 <- length(coef(gam2))
l.gam3 <- length(coef(gam3))
l.gam4 <- length(coef(gam4))
l.gam5 <- length(coef(gam5))
l.gam6 <- length(coef(gam6))
l.gam7 <- length(coef(gam7))


if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0) && fp==FALSE ){ 
  
                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7)
                 
                 qu.mag <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, gam7, 
                               l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, l.sp6, l.sp7, 
                               l.gam1, l.gam2, l.gam3, l.gam4, l.gam5, l.gam6, l.gam7) 
       
  }  
  


if(!(penCor %in% c("unpen")) && fp==FALSE){
    
    l.sp4  <- 1
    sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, spCor)  
  
########################################### 
# no penalty for corrs they are added after  
########################################### 
 
# here I set up rank and off as I can not do it later
# then I will incorporate the penalty for the corrs
 
if(l.sp1==0 && l.sp2==0 && l.sp3==0) qu.mag <- list(rank = 3, off = X1.d2 + X2.d2 + X3.d2 + 1, Ss = NULL) 
if(l.sp1!=0 || l.sp2!=0 || l.sp3!=0) qu.mag <- list(rank = c(qu.mag$rank, 3), off = c(qu.mag$off, X1.d2 + X2.d2 + X3.d2 + 1), Ss = qu.mag$Ss)

}


##########################################################


if(missing(parscale)) parscale <- 1   

  respvec <- list(y1 = y1,
                  y2 = y2,
                  y3 = y3,
  		  y1.y2.y3    = y1.y2.y3   , 
  		  y1.y2.cy3   = y1.y2.cy3  ,
  		  cy1.y2.y3   = cy1.y2.y3  , 
  		  cy1.y2.cy3  = cy1.y2.cy3 , 
       		  cy1.cy2.cy3 = cy1.cy2.cy3,
  		  cy1.cy2.y3  = cy1.cy2.y3 , 
                  y1.cy2.cy3  = y1.cy2.cy3 , 
                  y1.cy2.y3   = y1.cy2.y3,                    
                  cy1       = cy1,
                  y1.cy2    = y1.cy2,
                  y1.y2.cy3 = y1.y2.cy3, univ = 0)
 
   
  VC <- list(X1 = X1, inde1 = inde1, inde2 = inde2, inde2.1 = inde2.1,
             X2 = X2, 
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,
             X7 = X7,
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             X3.d2 = X3.d2,
             X4.d2 = X4.d2,
             X5.d2 = X5.d2,
             X6.d2 = X6.d2,   
             X7.d2 = X7.d2, 
             gp1 = gp1, 
             gp2 = gp2,
             gp3 = gp3,
             gp4 = gp4, 
             gp5 = gp5,
             gp6 = gp6,  
             gp7 = gp7,  
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6,  
             l.sp7 = l.sp7,
             infl.fac = infl.fac,
             weights = weights,
             fp = fp,
             hess = hess, nCa = nCa,
             Model = Model, gamlssfit = TRUE,
             end = NULL,
             BivD = "N", penCor = penCor,
             nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = c("probit","probit","probit"),
             Cont = "NO", ccss = "no", m2 = m2, m3 = m3, m2d = m2d, m3d = m3d, bl = bl, triv = TRUE,
             X2s = X2s, X3s = X3s,
             approx = approx, gamma = gamma, wc = w.alasso, qu.mag = qu.mag)
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################

  SemiParFit <- SemiParBIVProbit.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = VC$qu.mag) # there is a reason for calling qu.mag from VC
                                            
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParTRIV.fit.post(SemiParFit = SemiParFit, 
                                            VC = VC, Model = Model,
                                            gam1 = gam1, gam2 = gam2, gam3 = gam3,
                                            gam4 = gam4, gam5 = gam5, gam6 = gam6)
    
  SemiParFit <- SemiParFit.p$SemiParFit # useful for SS models, eta2 calculatons etc.    
    
  ##########################################################################################################################
                                            
                                            
if(gc.l == TRUE) gc()

e.v <- min(eigen(SemiParFit$fit$hessian, symmetric=TRUE, only.values = TRUE)$values)
gradi <- round(max(abs(SemiParFit$fit$gradient)),1)

me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?SemiParTRIProbit."

if(gradi > 10 && e.v <= 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 10 && e.v > 0)   warning(paste(me1,"\n",me3), call. = FALSE)
if(gradi < 10 && e.v <= 0)  warning(paste(me2,"\n",me3), call. = FALSE)

  ##########################################################################################################################

L <- list(fit = SemiParFit$fit, formula = formula, Model = Model,
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7,  
          coefficients = SemiParFit$fit$argument,  iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,l.sp7 = l.sp7, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta12   = SemiParFit.p$theta12, 
          theta12.a = SemiParFit.p$theta12.a,
          theta13   = SemiParFit.p$theta13, 
          theta13.a = SemiParFit.p$theta13.a,
          theta23   = SemiParFit.p$theta23, 
          theta23.a = SemiParFit.p$theta23.a,           
          n = n, 
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4 = X4, X5 = X5, X6 = X6, X7 = X7, X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2,           
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf11 = SemiParFit.p$edf11,   # what is this for?
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,
          p111 = SemiParFit$fit$p111,
          p011 = SemiParFit$fit$p011,
          p101 = SemiParFit$fit$p101,
          p110 = SemiParFit$fit$p110,
          p100 = SemiParFit$fit$p100,
          p010 = SemiParFit$fit$p010,
          p001 = SemiParFit$fit$p001,
          p000 = SemiParFit$fit$p000,
          eta1 = SemiParFit$fit$eta1, 
          eta2 = SemiParFit$fit$eta2,
          eta3 = SemiParFit$fit$eta3,
          y1 = y1, y2 = y2, y3 = y3,  n.sel1 = n.sel1, n.sel2 = n.sel2, 
          logLik = SemiParFit.p$logLik,
          nC = nC, hess = hess, 
          respvec = respvec, 
          qu.mag = qu.mag, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, 
          VC = VC, magpp = SemiParFit$magpp,
          gamlss = gamlss, gamlssfit = TRUE, Cont = "NO", triv = TRUE,  
          l.flist = l.flist, margins = VC$margins,
          inde2 = inde2, X2s = X2s, X3s = X3s,
          p1n = SemiParFit.p$p1n, p2n = SemiParFit.p$p2n, p3n = SemiParFit.p$p3n)

class(L) <- c("SemiParTRIV","SemiParBIVProbit")

L

}

