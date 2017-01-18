gamlss <- function(formula, data = list(), weights = NULL, subset = NULL,  
                   margin = "N", robust = FALSE, rc = 3, rsim = 100, infl.fac = 1, 
                   rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                   gc.l = FALSE, parscale, extra.regI = "t"){
  
  ##########################################################################################################################
  # preamble 
  ##########################################################################################################################
  
  
  
  
  i.rho <- sp <- qu.mag <- qu.mag1 <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- spgamlss1 <- NULL  
  end <- X2.d2 <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- y1m <- y2m <- NULL
  fp <- FALSE
  gp2 <- gp3 <- 1
  sp1 <- sp2 <- gam2 <- X2 <- sp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- NULL     

  m2   <- c("N","N2","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK")
  m3   <- c("DAGUM","SM")
  m1d  <- c("PO", "ZTP")
  m2d  <- c("NBI", "NBII","NBIa", "NBIIa","PIG")
  m3d  <- c("DEL","SICHEL")
  M    <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, robust = robust, extra.regI = extra.regI, margin = margin)

  ##########################################################################################################################
  if(!is.list(formula)) stop("You must specify a list of equations.")
  l.flist <- length(formula)
  pream.wm(formula, margins = NULL, M, l.flist, type = "gamls")
  form.check(formula, l.flist, gamlss = TRUE) 
  
  mf <- match.call(expand.dots = FALSE)

  pred.varR <- pred.var(formula, l.flist, gaml = TRUE) 
   
  v1     <- pred.varR$v1  
  v2     <- pred.varR$v2
  pred.n <- pred.varR$pred.n  

  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  mf$rsim <- mf$robust <- mf$rc <- mf$margin <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- NULL                           
  mf$drop.unused.levels <- TRUE 
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
  if(gc.l == TRUE) gc()  
 
  n <- dim(data)[1]
        
  if(is.null(weights)) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  
 ##############################################################  
 ##############################################################  
   
 form.eq12R <- form.eq12(formula.eq1, data, v1, margin, m1d, m2d)   
 
 formula.eq1  <- form.eq12R$formula.eq1
 formula.eq1r <- form.eq12R$formula.eq1r
 y1           <- form.eq12R$y1
 y1.test      <- form.eq12R$y1.test  
 y1m          <- form.eq12R$y1m
   
 gam1         <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
 gam1$formula <- formula.eq1r  
 names(gam1$model)[1] <- as.character(formula.eq1r[2])
     
 y1 <- y1.test 
 if( margin %in% c("LN") ) y1 <- log(y1) 
    
 X1 <- model.matrix(gam1)
 X1.d2 <- dim(X1)[2]
 l.sp1 <- length(gam1$sp)
 if(l.sp1 != 0) sp1 <- gam1$sp
 gp1 <- gam1$nsdf 
                      
############

log.nu.1 <- log.sig2.1 <- NULL 

if( !(margin %in% c(m1d)) ){

start.snR <- startsn(margin, y1)
    
log.sig2.1 <- start.snR$log.sig2.1; names(log.sig2.1) <- "sigma2.star"
if( margin %in% c(m3) ){ log.nu.1   <- start.snR$log.nu.1;   names(log.nu.1)   <- "nu.star"}     

}


if(margin %in% c(m1d)    ) start.v1 <- c( coef(gam1)                       )
if(margin %in% c(m2,m2d) ) start.v1 <- c( coef(gam1), log.sig2.1           ) 
if(margin %in% c(m3,m3d) ) start.v1 <- c( coef(gam1), log.sig2.1, log.nu.1 ) 

##############################################################  
##############################################################  
  
if(l.flist > 1){
    
    vo <- list(log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, n = n)
    overall.svGR <- overall.svG(formula, data, ngc = 2, margin, M, vo, gam1, gam2, type = "gaml")
    
    start.v1 <- overall.svGR$start.v 
    X2 <- overall.svGR$X2
    X3 <- overall.svGR$X3  
    X2.d2 <- overall.svGR$X2.d2
    X3.d2 <- overall.svGR$X3.d2
    gp2 <- overall.svGR$gp2
    gp3 <- overall.svGR$gp3
    gam2 <- overall.svGR$gam2
    gam3 <- overall.svGR$gam3
    l.sp2 <- overall.svGR$l.sp2
    l.sp3 <- overall.svGR$l.sp3
    sp2 <- overall.svGR$sp2
    sp3 <- overall.svGR$sp3
    
}
  

##########################################################
# SPs and penalties
##########################################################

spgamlss1 <- c(sp1, sp2, sp3)
GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8) 


if(l.sp1 !=0 || l.sp2 !=0 || l.sp3 !=0) { ##


L.GAM <- list(l.gam1 = length(coef(gam1)), l.gam2 = length(coef(gam2)), 
              l.gam3 = length(coef(gam3)), l.gam4 = 0, l.gam5 = 0, 
              l.gam6 = 0, l.gam7 = 0, l.gam8 = 0)                                 
  
L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = 0, 
             l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0)                               
                      
qu.mag1 <- S.m(GAM, L.SP, L.GAM)  

} ##


##########################################################
##########################################################

if(missing(parscale)) parscale <- 1   

respvec2 <- list(y1 = y1, univ = 2)               
    
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
             l.sp8 = l.sp8, 
             infl.fac = infl.fac, rsim = rsim,
             weights = weights,
             fp = fp, 
             hess = NULL,
             Model = "CC", univ.gamls = TRUE,
             gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = c(margin, margin),
             Cont = "YES", ccss = "no", m2 = m2, m3 = m3, m1d = m1d, 
             m2d = m2d, m3d = m3d, bl = NULL, triv = FALSE,
             y1m = y1m, y2m = y2m, robust = robust, rc = rc) 
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################

  if(margin %in% c(m1d, m2d, m2) ) func.opt1 <- bprobgHsContUniv else func.opt1 <- bprobgHsContUniv3


    SemiParFit <- SemiParBIVProbit.fit(func.opt = func.opt1, start.v = start.v1, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec2, VC = VC, sp = spgamlss1, qu.mag = qu.mag1) 
                                  
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- gamlss.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)  
                                                          
  y1.m <- y1; if(margin == "LN") y1.m <- exp(y1) 

  SemiParFit <- SemiParFit.p$SemiParFit  

  ##########################################################################################################################

if(gc.l == TRUE) gc()

  ##########################################################################################################################


e.v <- min(eigen(SemiParFit$fit$hessian, symmetric=TRUE, only.values = TRUE)$values)
gradi <- round(max(abs(SemiParFit$fit$gradient)),1)

me1 <- "Largest absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?gamlss."

if(gradi > 10 && e.v <= 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 10 && e.v > 0)   warning(paste(me1,"\n",me3), call. = FALSE)
if(gradi < 10 && e.v <= 0)  warning(paste(me2,"\n",me3), call. = FALSE)


  ##########################################################################################################################



L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, formula = formula,        
          edf11 = SemiParFit.p$edf11,   ## this is for RE
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, 
          gam6 = gam6, gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, iterlimsp = iterlimsp,
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
          l.sp7 = l.sp7, l.sp8 = l.sp8,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,  
          sigma2 = SemiParFit.p$sigma2,  
          sigma2.a = SemiParFit.p$sigma2.a, 
          nu = SemiParFit.p$nu, 
          nu.a = SemiParFit.p$nu.a,
          X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, 
          X6 = X6, X7 = X7, X8 = X8,
          X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, 
          X7.d2 = X7.d2, X8.d2 = X8.d2,             
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, 
          edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
          edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$etas1, 
          eta3 = SemiParFit$fit$etan1, 
          y1 = y1.m, 
          margins = c(margin, margin),   
          logLik = SemiParFit.p$logLik,
          hess = TRUE,
          qu.mag = qu.mag1, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, 
          gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          VC = VC, magpp = SemiParFit$magpp,
          Cont = "YES",
          l.flist = l.flist, triv = FALSE, univar.gamlss = TRUE)

class(L) <- c("gamlss","SemiParBIVProbit")


L


}

