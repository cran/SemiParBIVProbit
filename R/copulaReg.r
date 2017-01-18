copulaReg <- function(formula, data = list(), weights = NULL, subset = NULL,  
                             BivD = "N", margins = c("N", "N"), dof = 3, 
                             gamlssfit = FALSE, fp = FALSE, infl.fac = 1, 
                             rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t"){
  
  bl   <- c("probit", "logit", "cloglog", "cauchit")  
  
  if(margins[1] %in% bl){
                          
  L <- eval(substitute(SemiParBIVProbit(formula, data, weights, subset,
                               Model = "B", BivD, margins, dof, gamlssfit,
                               fp, hess = TRUE, infl.fac, 
                               rinit, rmax, iterlimsp, tolsp,
                               gc.l, parscale, extra.regI, intf = TRUE, 
                               theta.fx = NULL),list(weights=weights)))                               
  
  }
  
  if(!(margins[1] %in% bl)){
    
  ##########################################################################################################################
  # preamble
  ##########################################################################################################################  
  robust <- FALSE; t.c = 3
  sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <-  gam1 <- gam2 <- y1m <- y2m <- NULL  
  i.rho <- log.sig2.2 <- log.nu.2 <- log.nu.1 <- log.sig2.1 <- dof.st <- NULL
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- NULL    
  sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- NULL   
  
  BivD2 <- c("C0C90","C0C270","C180C90","C180C270",
             "J0J90","J0J270","J180J90","J180J270",
             "G0G90","G0G270","G180G90","G180G270")
             
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM","T")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180", BivD2)
  sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")
  m2   <- c("N","N2","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK")
  m3   <- c("DAGUM","SM")
  m1d  <- c("PO", "ZTP")
  m2d  <- c("NBI", "NBII","NBIa", "NBIIa","PIG")
  m3d  <- c("DEL","SICHEL")
  

  
  M    <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, BivD = BivD, 
               robust = robust, opc = opc, extra.regI = extra.regI, margins = margins, BivD2 = BivD2, dof = dof)

  ct  <- data.frame( c(opc), c(1:14,55,56,57) )
  cta <- data.frame( c(opc), c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56,2) )     
  
  
  
  if(BivD %in% BivD2){
  
  if(BivD %in% BivD2[1:4])  BivDt <- "C0" 
  if(BivD %in% BivD2[5:12]) BivDt <- "J0"
  
  nC  <-  ct[which( ct[,1]==BivDt),2]
  nCa <- cta[which(cta[,1]==BivDt),2]     
  
  }
  
  
  if(!(BivD %in% BivD2)){
    
  nC  <-  ct[which( ct[,1]==BivD),2]
  nCa <- cta[which(cta[,1]==BivD),2]     
    
  }
  
  


 #######################################################################################  
  if(!is.list(formula)) stop("You must specify a list of equations.")
  l.flist <- length(formula)
  pream.wm(formula, margins, M, l.flist)
  form.check(formula, l.flist) 
       
  mf <- match.call(expand.dots = FALSE)
            
  pred.varR <- pred.var(formula, l.flist) 
   
  v1     <- pred.varR$v1  
  v2     <- pred.varR$v2
  pred.n <- pred.varR$pred.n  

  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$BivD <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
  mf$drop.unused.levels <- TRUE 
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  
 
  n <- dim(data)[1]
        
  if(is.null(weights)) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
    
 ##############################################################  
 # Equation 1
 ##############################################################  
   
 form.eq12R <- form.eq12(formula.eq1, data, v1, margins[1], m1d, m2d)   
 
 formula.eq1  <- form.eq12R$formula.eq1
 formula.eq1r <- form.eq12R$formula.eq1r
 y1           <- form.eq12R$y1
 y1.test      <- form.eq12R$y1.test 
 y1m          <- form.eq12R$y1m

 gam1         <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
 gam1$formula <- formula.eq1r  
 names(gam1$model)[1] <- as.character(formula.eq1r[2])
  
 y1 <- y1.test 
 if( margins[1] %in% c("LN") ) y1 <- log(y1) 
 
 X1 <- model.matrix(gam1)
 X1.d2 <- dim(X1)[2]
 l.sp1 <- length(gam1$sp)
 if(l.sp1 != 0) sp1 <- gam1$sp
 gp1 <- gam1$nsdf 
    
 ##############################################################
 # Equation 2 
 ##############################################################  

 form.eq12R <- form.eq12(formula.eq2, data, v2, margins[2], m1d, m2d)   
 
 formula.eq2  <- form.eq12R$formula.eq1
 formula.eq2r <- form.eq12R$formula.eq1r
 y2           <- form.eq12R$y1
 y2.test      <- form.eq12R$y1.test 
 y2m          <- form.eq12R$y1m

 gam2         <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
 gam2$formula <- formula.eq2r  
 names(gam2$model)[1] <- as.character(formula.eq2r[2])    
 
 y2 <- y2.test 
 if( margins[2] %in% c("LN") ) y2 <- log(y2) 
  
 X2 <- model.matrix(gam2)
 X2.d2 <- dim(X2)[2]
 l.sp2 <- length(gam2$sp)
 if(l.sp2 != 0) sp2 <- gam2$sp 
 gp2 <- gam2$nsdf    
    
#################################################################
# Starting value for dependence parameter (and dof for T if used)
#################################################################

res1 <- residuals(gam1)
res2 <- residuals(gam2)
                  
ass.s <- cor(res1, res2, method = "kendall")
ass.s <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))

i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)

dof.st <- log(dof - 2) 
names(dof.st) <- "dof.star"   
                           
##############################################################
# Other starting values + overall
##############################################################
           
if( !(margins[1] %in% c(m1d)) ){

start.snR <- startsn(margins[1], y1)
    
log.sig2.1 <- start.snR$log.sig2.1; names(log.sig2.1) <- "sigma2.1.star"
if( margins[1] %in% c(m3) ){ log.nu.1   <- start.snR$log.nu.1;   names(log.nu.1)   <- "nu.1.star"}     

}

if( !(margins[2] %in% c(m1d)) ){

start.snR <- startsn(margins[2], y2)
    
log.sig2.2 <- start.snR$log.sig2.1; names(log.sig2.2) <- "sigma2.2.star"
if( margins[2] %in% c(m3) ){ log.nu.2   <- start.snR$log.nu.1;   names(log.nu.2)   <- "nu.2.star"}     

}

vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, log.sig2.2 = log.sig2.2, 
           log.nu.2 = log.nu.2, log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, dof.st = dof.st, n = n )
start.v <- overall.sv(margins, M, vo)
   			
##############################################################  
# starting values for case of predictors on all parameters
##############################################################  
  
    if(l.flist > 2){
    
    overall.svGR <- overall.svG(formula, data, ngc = 2, margins, M, vo, gam1, gam2)
    
    start.v = overall.svGR$start.v 
    X3 = overall.svGR$X3; X4 = overall.svGR$X4; X5 = overall.svGR$X5
    X6 = overall.svGR$X6; X7 = overall.svGR$X7; X8 = overall.svGR$X8
    X3.d2 = overall.svGR$X3.d2; X4.d2 = overall.svGR$X4.d2; X5.d2 = overall.svGR$X5.d2
    X6.d2 = overall.svGR$X6.d2; X7.d2 = overall.svGR$X7.d2; X8.d2 = overall.svGR$X8.d2
    gp3 = overall.svGR$gp3; gp4 = overall.svGR$gp4; gp5 = overall.svGR$gp5
    gp6 = overall.svGR$gp6; gp7 = overall.svGR$gp7; gp8 = overall.svGR$gp8
    gam3 = overall.svGR$gam3; gam4 = overall.svGR$gam4; gam5 = overall.svGR$gam5
    gam6 = overall.svGR$gam6; gam7 = overall.svGR$gam7; gam8 = overall.svGR$gam8
    l.sp3 = overall.svGR$l.sp3; l.sp4 = overall.svGR$l.sp4; l.sp5 = overall.svGR$l.sp5
    l.sp6 = overall.svGR$l.sp6; l.sp7 = overall.svGR$l.sp7; l.sp8 = overall.svGR$l.sp8
    sp3 = overall.svGR$sp3; sp4 = overall.svGR$sp4; sp5 = overall.svGR$sp5
    sp6 = overall.svGR$sp6; sp7 = overall.svGR$sp7; sp8 = overall.svGR$sp8
    
    }
    

##########################################################
# SPs and penalties
##########################################################
  

GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8)   


if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0) && fp==FALSE ){ 

L.GAM <- list(l.gam1 = length(coef(gam1)), l.gam2 = length(coef(gam2)), l.gam3 = length(coef(gam3)), l.gam4 = length(coef(gam4)),
              l.gam5 = length(coef(gam5)), l.gam6 = length(coef(gam6)), l.gam7 = length(coef(gam7)), l.gam8 = length(coef(gam8)))

L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
             l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8)

                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
                 qu.mag <- S.m(GAM, L.SP, L.GAM)                             
                                                        }
  

##########################################################
# general lists
##########################################################

if(missing(parscale)) parscale <- 1   

  respvec <- respvec2 <- respvec3 <- list(y1 = y1, y2 = y2,
                  y1.y2 = NULL, y1.cy2 = NULL, 
                  cy1.y2 = NULL, cy1.cy2 = NULL, 
                  cy1 = NULL, cy = NULL, univ = 0)
 
  my.env <- new.env()
  my.env$signind <- 1


  VC <- list(X1 = X1, 
             X2 = X2, 
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,  
             X7 = X7, 
             X8 = X8,
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             X3.d2 = X3.d2,
             X4.d2 = X4.d2,
             X5.d2 = X5.d2,
             X6.d2 = X6.d2,
             X7.d2 = X7.d2,
             X8.d2 = X8.d2,
             gp1 = gp1, 
             gp2 = gp2,
             gp3 = gp3,
             gp4 = gp4, 
             gp5 = gp5,
             gp6 = gp6, 
             gp7 = gp7,
             gp8 = gp8,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6, 
             l.sp7 = l.sp7, 
             l.sp8 = l.sp8, my.env = my.env,
             infl.fac = infl.fac,
             weights = weights,
             fp = fp, 
             gamlssfit = gamlssfit,
             hess = NULL,
             Model = "CC", univ.gamls = FALSE,
             end = end,
             BivD = BivD, nCa = nCa,
             nC = nC, gc.l = gc.l, 
             n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "YES", ccss = "no", m2 = m2, m3 = m3, 
             m1d = m1d, m2d = m2d, m3d = m3d, 
             bl = bl, triv = FALSE,
             y1m = y1m, y2m = y2m, 
             tc = t.c,
             i.rho = i.rho, dof = dof,
             dof.st = dof.st, BivD2 = BivD2, cta = cta, ct = ct,
             zerov = -10) # original n only needed in SemiParBIVProbit.fit
  

  
  
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  ##########################################################################################################################
  # GAMLSS fit
  ##########################################################################################################################
  ##########################################################################################################################

if(gamlssfit == TRUE){ 

  form.gamlR <- form.gaml(formula, l.flist, M)

  gamlss1 <- eval(substitute(gamlss(form.gamlR$formula.gamlss1, data = data, weights = weights, subset = subset,  
                   margin = margins[1], infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI), list(weights=weights)))

  gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = subset,  
                   margin = margins[2], infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI), list(weights=weights)))   
                      
  # updated starting values   

  SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)
  gamls.upsvR <- gamls.upsv(gamlss1, gamlss2, margins, M, l.flist, nstv = names(start.v), VC, GAM, SP)
  sp <- gamls.upsvR$sp
  start.v <- gamls.upsvR$start.v 

}

  ##########################################################################################################################
  ##########################################################################################################################
  
  func.opt <- func.OPT(margins, M)  
  
  SemiParFit <- SemiParBIVProbit.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
    
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- copulaReg.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)                                     
 
  y1.m <- y1; if(margins[1] == "LN") y1.m <- exp(y1) 
  y2.m <- y2; if(margins[2] == "LN") y2.m <- exp(y2)

  SemiParFit <- SemiParFit.p$SemiParFit  

  if(gc.l == TRUE) gc()

  ##########################################################################################################################


e.v <- min(eigen(SemiParFit$fit$hessian, symmetric=TRUE, only.values = TRUE)$values)
gradi <- round(max(abs(SemiParFit$fit$gradient)),1)

me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?copulaReg."

if(gradi > 10 && e.v <= 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 10 && e.v > 0)   warning(paste(me1,"\n",me3), call. = FALSE)
if(gradi < 10 && e.v <= 0)  warning(paste(me2,"\n",me3), call. = FALSE)


  ##########################################################################################################################


L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, gamlss1 = gamlss1, gamlss2 = gamlss2, formula = formula,        
          edf11 = SemiParFit.p$edf11,   
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,  
          sigma21 = SemiParFit.p$sigma21, sigma22 = SemiParFit.p$sigma22, 
          sigma21.a = SemiParFit.p$sigma21.a, sigma22.a = SemiParFit.p$sigma22.a,
          nu1 = SemiParFit.p$nu1, nu2 = SemiParFit.p$nu2, 
          nu1.a = SemiParFit.p$nu1.a, nu2.a = SemiParFit.p$nu2.a,
          dof.a = SemiParFit.p$dof.a, dof = SemiParFit.p$dof,
          X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8,
          X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2,            
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, edf7 = SemiParFit.p$edf7,
          edf8 = SemiParFit.p$edf8,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, edf1.7 = SemiParFit.p$edf1.7, 
          edf1.8 = SemiParFit.p$edf1.8, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
          etad=SemiParFit$fit$etad, etas1 = SemiParFit$fit$etas1, etas2 = SemiParFit$fit$etas2,
          y1 = y1.m, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, 
          respvec = respvec, hess = TRUE,
          qu.mag = qu.mag, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          VC = VC, magpp = SemiParFit$magpp,
          gamlssfit = gamlssfit, Cont = "YES",
          tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE,
          BivD2 = BivD2)
  
if(BivD %in% BivD2){       

L$teta1     <- SemiParFit$fit$teta1
L$teta.ind1 <- SemiParFit$fit$teta.ind1   
L$teta2     <- SemiParFit$fit$teta2
L$teta.ind2 <- SemiParFit$fit$teta.ind2   
L$Cop1      <- SemiParFit$fit$Cop1
L$Cop2      <- SemiParFit$fit$Cop2

}          

class(L) <- c("copulaReg","SemiParBIVProbit")


}

L


}

