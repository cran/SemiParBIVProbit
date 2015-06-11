SemiParBIVProbit <- function(formula, data = list(), weights = NULL, subset = NULL, start.v = NULL, 
                             Model = "B", BivD = "N", margins = c("probit","probit"), gamlssfit = FALSE,
                             fp = FALSE, hess = TRUE, infl.fac = 1, 
                             rinit = 1, rmax = 100, iterlimsp = 50, pr.tolsp = 1e-6,
                             gc.l = FALSE, parscale, extra.regI = "t"){
  
  ##########################################################################################################################
  # model set up and starting values
  ##########################################################################################################################
  
  i.rho <- sp <- qu.mag <- qu.mag1 <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss <- NULL  
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- 0
  sp3 <- gp3 <- gam3 <- X3 <- NULL  
  sp4 <- gp4 <- gam4 <- X4 <- NULL  
  sp5 <- gp5 <- gam5 <- X5 <- NULL   
  sp6 <- gp6 <- gam6 <- X6 <- NULL   
    
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180")
  sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")
  mb   <- c("B", "BSS", "BPO", "BPO0")
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA")
  

  if(margins[2]%in% m2) stop("Check next release for final tested version of this model.")
  
  if(Model == "BPO" && BivD != "N") stop("This model is not defined for copulae.")
  
  if(!(Model %in% mb)) stop("Error in parameter Model value. It should be one of: B, BSS, BPO, BPO0.")
  if(!(BivD %in% opc)) stop("Error in parameter BivD value. It should be one of: N, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270, F.")
  if(!(extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of: t, pC or sED.")
  
  if(margins[1] != "probit" ) stop("Error in first margin value. It can only be: probit.")
  if(!(margins[2] %in% c("probit",m2)) ) stop("Error in second margin value. It can be: probit, N, GU, rGU, LO, LN, WEI, iG, GA.")  
  if(margins[2] %in% m2 && (Model == "BPO" || Model == "BSS" || Model == "BPO0") ) stop("For continuous responses, only bivariate models are allowed for.")   
  
  if(length(formula) > 2 && margins[2] %in% m2){ if(length(formula)!=4) stop("You need to specify four equations.") }  
  if( length(formula) == 3 && Model == "BPO0") stop("The chosen model does not have a correlation parameter.")
  if( length(formula) > 3 && margins[2] == "probit") stop("The chosen model can not have more than three equations.")
  if( length(formula) > 4 && margins[2] != "probit") stop("The chosen model can not have more than four equations.")

 
 #######################################################################################  
 # formula check  
 #######################################################################################  
  
  
    if(length(formula) > 2){
    
    f3t <- try(formula[[3]][[3]], silent = TRUE)  
    if(class(f3t)!="try-error") stop("The third equation does not require a response.")
    
    	if(length(formula) > 3){
    	
    f4t <- try(formula[[4]][[3]], silent = TRUE)  
    if(class(f4t)!="try-error") stop("The fourth equation does not require a response.")  
    
    		if(length(formula) > 4){
    		
    f5t <- try(formula[[5]][[3]], silent = TRUE)  
    if(class(f5t)!="try-error") stop("The fifth equation does not require a response.")    
    
    			if(length(formula) > 5){
    			
    f6t <- try(formula[[6]][[3]], silent = TRUE)  
    if(class(f6t)!="try-error") stop("The sixth equation does not require a response.")      			
    			
    			                        }    
    					}
    				}    				
                            }

 #######################################################################################  
  
 
  ig <- interpret.gam(formula)
  mf <- match.call(expand.dots = FALSE)

  if( length(formula) == 2 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[2]]$response))
  if( length(formula) == 3 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[2]]$response))
  if( length(formula) == 4 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[4]]$pred.names, ig[[2]]$response))
  if( length(formula) == 5 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[4]]$pred.names, ig[[5]]$pred.names, ig[[2]]$response))
  if( length(formula) == 6 ) pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names, ig[[3]]$pred.names, ig[[4]]$pred.names, ig[[5]]$pred.names, ig[[6]]$pred.names, ig[[2]]$response))
  
  
  fake.formula <- paste(ig[[1]]$response, "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(ig$fake.formula)
  mf$formula <- fake.formula 
  mf$start.v <- mf$Model <- mf$BivD <- mf$margins <- mf$fp <- mf$hess <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$pr.tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
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
  
  if(is.null(weights)) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"} else weights <- data[,"(weights)"]    
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
  
  
  if(Model=="B"){
  if(ig[[1]]$response %in% ig[[2]]$pred.names ) end <- 1
  if(ig[[2]]$response %in% ig[[1]]$pred.names ) end <- 2
  }

  ct  <- data.frame( c(opc),
                    c(1:14) 
                     )
  cta <- data.frame( c(opc),
                     c(1,3,23,13,33,6,26,16,36,4,24,14,34,5) 
                     )                   
  nC  <- ct[which(ct[,1]==BivD),2]
  nCa <- cta[which(cta[,1]==BivD),2]  
  
  
 ##############################################################  
 # Equation 1
 ##############################################################  
   
  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), gamma=infl.fac, weights=weights, 
                              data=data),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 
  if(l.sp1 != 0) sp1 <- gam1$sp else sp1 <- NULL 

 ##############################################################
 # Equation 2 for BPO and binary B
 ##############################################################  

  if((Model=="B" || Model=="BPO" || Model=="BPO0") && margins[2]=="probit"){
  
  gam2  <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=infl.fac, weights=weights, 
                           data=data),list(weights=weights))) # check at later stage the need of eval and substitute

  X2 <- model.matrix(gam2)
  X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp)
  y2 <- gam2$y 
  if(l.sp2 != 0) sp2 <- gam2$sp else sp2 <- NULL 
     
  if(Model=="B"){  
  	y1.y2 <- y1*y2
  	y1.cy2 <- y1*(1-y2)
  	cy1.y2 <- (1-y1)*y2
  	cy1.cy2 <- (1-y1)*(1-y2)
                }

  if(Model=="BPO" || Model=="BPO0" ) cy <- 1 - y1
  
  } 
  
 ##############################################################
 # Equation 2 for B and continuous response 
 ##############################################################  

  if(Model=="B" && margins[2] != "probit" ){
  
  
  if(margins[2] %in% c("N","GU","rGU","LO") ){
  
    gam2  <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
      
  }
  
  if(margins[2] == "LN"){
  
    formula.eq2r <- formula.eq2 
    formula.eq2  <- update(formula.eq2, log(.) ~ . )
    gam2         <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data),list(weights=weights)))
    gam2$formula <- formula.eq2r 
    
  } 

  if(margins[2] %in% c("WEI", "iG", "GA", "iGA") ){
  
  gam2  <- eval(substitute(gam(formula.eq2, family = Gamma(link = "log"), gamma=infl.fac, weights=weights, data=data),list(weights=weights)))

  }
  
    
    y2 <- gam2$y
    X2 <- model.matrix(gam2)
    X2.d2 <- dim(X2)[2]
    l.sp2 <- length(gam2$sp)
    if(l.sp2 != 0) sp2 <- gam2$sp else sp2 <- NULL 
    
    cy <- 1 - y1
    
  } 

 ##############################################################
 # Equation 2 for BSS 
 ##############################################################  

  if(Model=="BSS"){

  inde <- y1 > 0

  gam2 <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=infl.fac, weights=weights, 
                              data=data, subset=inde),list(weights=weights,inde=inde)))  
                              
  X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2) 
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y   # ; n.sel <- length(gam2$y)
  l.sp2 <- length(gam2$sp)
  if(l.sp2 != 0) sp2 <- gam2$sp else sp2 <- NULL 
  
  cy1 <- (1-y1)
  y1.y2 <- y1*y2
  y1.cy2 <- y1*(1-y2)

  }
  
##############################################################  
  
  gp1 <- gam1$nsdf 
  gp2 <- gam2$nsdf
  
##############################################################
# Starting values for dependence parameter
##############################################################

# starting value of association
    
if(is.null(start.v)){    
    
if( !(Model %in% c("BPO","BPO0")) ){    
    
if(Model=="B")    res1 <- residuals(gam1)
if(Model=="BSS")  res1 <- residuals(gam1)[inde]
                  res2 <- residuals(gam2)
ass.s <- cor(res1,res2) + 0.0000001

if(BivD %in% scc)  ass.s <-  abs(ass.s)   
if(BivD %in% sccn) ass.s <- -abs(ass.s) 

i.rho <- BiCopTau2Par(family = nCa, tau = ass.s)

if(BivD == "N") i.rho <- atanh( i.rho ) 
if(BivD == "F") i.rho <- i.rho + 0.0000001
if(!(BivD %in% c("N","F"))) i.rho <- abs(i.rho)

if(BivD %in% c("C0","C180","C90","C270"))                            i.rho <-  log(i.rho)   
if(BivD %in% c("J0","J180","G0","G180","J90","J270","G90","G270"))   i.rho <-  log(i.rho - 1)   

}


if(Model=="BPO"){ 

ass.s <- 0.1

if(BivD %in% scc)  ass.s <-  abs(ass.s)   
if(BivD %in% sccn) ass.s <- -abs(ass.s) 

i.rho <- BiCopTau2Par(family = nCa, tau = ass.s)

if(BivD == "N") i.rho <- atanh( i.rho ) 
if(BivD == "F") i.rho <- i.rho + 0.0000001
if(!(BivD %in% c("N","F"))) i.rho <- abs(i.rho)

if(BivD %in% c("C0","C180","C90","C270"))                            i.rho <-  log(i.rho)   
if(BivD %in% c("J0","J180","G0","G180","J90","J270","G90","G270"))   i.rho <-  log(i.rho - 1)   

}

if(Model=="BPO0")  i.rho <- 0

names(i.rho) <- "theta.star"   


}

                           
##############################################################
# Starting values for whole parameter vector
##############################################################
                      
    if(length(formula) == 2 && margins[1] == "probit" && margins[2] == "probit" && Model != "BPO0" ){
    
       if(is.null(start.v)) start.v <- c(coef(gam1), coef(gam2), i.rho) 
       
    } 
    
    if(margins[1] == "probit" && margins[2] != "probit" ){
    
       log.sig2 <- log(gam2$sig2); names(log.sig2) <- "sigma2.star"
       if(is.null(start.v)) start.v <- c(coef(gam1), coef(gam2), log.sig2, i.rho) 
       
                           start.v1 <- c(            coef(gam2), log.sig2       ) 
       
       
    } 
    
    if(Model == "BPO0"){
    
       if(is.null(start.v)) start.v <- c(coef(gam1), coef(gam2)) 
       
    }     
    
    
##############################################################  
  
    if(length(formula)==3){
    
    formula.eq3 <- formula[[3]] 
    
    set.seed(0)
    theta <- rnorm(length(y2)) # seq(-4, 4, length.out = length(y2))
    
    nad <- "theta"  
    
    formula.eq3 <- as.formula( paste(nad,"~",formula.eq3[2],sep="") ) 

    gam3 <- eval(substitute(gam(formula.eq3, weights=weights, data=data))) 
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]

    l.sp3 <- length(gam3$sp)
    if(l.sp3 != 0) sp3 <- gam3$sp else sp3 <- NULL 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
  
 
    names(i.rho) <- names(coef(gam3))[1]
    g3c <- rep(0, length(coef(gam3))-1 )
    names(g3c) <- names(coef(gam3))[-1]
    
    
    start.v <- c( coef(gam1), coef(gam2), i.rho, g3c )
    
  }
  
  
    if(length(formula)==4){
    
    formula.eq3 <- formula[[3]] 
    formula.eq4 <- formula[[4]] 
     
    set.seed(0) 
    sigma2 <- rnorm(length(y2)) # seq(-4, 4, length.out = length(y2)); nad1 <- "sigma2" 
    theta  <- rnorm(length(y2)) # seq(-4, 4, length.out = length(y2)); nad2 <- "theta"  
    
    nad1 <- "sigma2" 
    nad2 <- "theta" 
    
    formula.eq3 <- as.formula( paste(nad1,"~",formula.eq3[2],sep="") ) 
    formula.eq4 <- as.formula( paste(nad2,"~",formula.eq4[2],sep="") ) 

    gam3 <- eval(substitute(gam(formula.eq3, weights=weights, data=data))) 
    gam4 <- eval(substitute(gam(formula.eq4, weights=weights, data=data)))     
    X3 <- model.matrix(gam3)
    X3.d2 <- dim(X3)[2]
    X4 <- model.matrix(gam4)
    X4.d2 <- dim(X4)[2]    

    l.sp3 <- length(gam3$sp)
    if(l.sp3 != 0) sp3 <- gam3$sp else sp3 <- NULL 
    environment(gam3$formula) <- environment(gam2$formula)
    gp3 <- gam3$nsdf 
    
    l.sp4 <- length(gam4$sp)
    if(l.sp4 != 0) sp4 <- gam4$sp else sp4 <- NULL 
    environment(gam4$formula) <- environment(gam2$formula)
    gp4 <- gam4$nsdf     
  
  
    names(log.sig2) <- names(coef(gam3))[1]
    g3c <- rep(0, length(coef(gam3))-1 )
    names(g3c) <- names(coef(gam3))[-1]  
 
    names(i.rho) <- names(coef(gam4))[1]
    g4c <- rep(0, length(coef(gam4))-1 )
    names(g4c) <- names(coef(gam4))[-1]
    
    
    start.v  <- c(coef(gam1), coef(gam2), log.sig2, g3c , i.rho, g4c )
    start.v1 <- c(            coef(gam2), log.sig2, g3c              ) # introduced for univariate contin outcome model
    
  }  
  
  
  
  
##########################################################
  
  
  if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0) && fp==FALSE ){ 
  
                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6)
                 qu.mag <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, l.sp6) 
       
       
  if(margins[2] %in% m2 && gamlssfit == TRUE ) {qu.mag1 <- S.m(gam1, gam2, gam3, gam4, gam5, gam6, 
                                          l.sp1, l.sp2, l.sp3, l.sp4, l.sp5, l.sp6, 
                                          eq1 = "no")                                                                       
                           sp1 <- c(sp2, sp3, sp4, sp5, sp6)                                                          
                          }
                 
                                                                                           }

##########################################################


if(missing(parscale)) parscale <- 1   


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
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,             
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             X3.d2 = X3.d2,
             X4.d2 = X4.d2,
             X5.d2 = X5.d2,
             X6.d2 = X6.d2,             
             gp1 = gp1, 
             gp2 = gp2,
             gp3 = gp3,
             gp4 = gp4, 
             gp5 = gp5,
             gp6 = gp6,             
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6,             
             infl.fac = infl.fac,
             weights = weights,
             fp = fp,
             hess = hess,
             Model = Model,
             end = end,
             BivD = BivD,
             nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins) # original n only needed in SemiParBIVProbit.fit
             
  if(gc.l == TRUE) gc()           
             
  ##########################################################################################################################
  # model fitting
  ##########################################################################################################################

  if(Model=="B" && margins[2]=="probit")   func.opt <- bprobgHs   
  if(Model=="B" && margins[2]%in% m2)      func.opt <- bprobgHsCont   
  
  if(Model=="BPO")                         func.opt <- bprobgHsPO 
  if(Model=="BPO0")                        func.opt <- bprobgHsPO0   
  if(Model=="BSS")                         func.opt <- bprobgHsSS 



  if(margins[2] %in% m2 && gamlssfit == TRUE){ 
  
  gamlss <- SemiParBIVProbit.fit(func.opt = bprobgHsContUniv, start.v = start.v1, 
                         rinit = rinit, rmax = rmax, iterlim = 1e+4, iterlimsp = iterlimsp, pr.tolsp = pr.tolsp,
                         respvec = respvec, VC = VC, sp = sp1, qu.mag = qu.mag1, naive = TRUE)
                         
                         
                         
  # new starting values                       
                         
  if( length(formula) == 2 ) start.v <- c(coef(gam1), gamlss$fit$argument, i.rho)
  if( length(formula) == 4 ) start.v <- c(coef(gam1), gamlss$fit$argument, i.rho, g4c)
  if(l.sp2 != 0) sp2 <- gamlss$sp[1:l.sp2] 
  if(l.sp3 != 0) sp3 <- gamlss$sp[l.sp2 + (1:l.sp3)] 
  sp <- c(sp1, sp2, sp3, sp4, sp5, sp6)
  
  
  # change gam2 output
  #if(margins[2]=="N")   {fam <- "Gaussian";       ; lin <- "identity"    }
  #if(margins[2]=="GU")  {fam <- "Gumbel"          ; lin <- "identity"    }
  #if(margins[2]=="rGU") {fam <- "reverse Gumbel"  ; lin <- "identity"    }
  #if(margins[2]=="LO")  {fam <- "logistic"        ; lin <- "identity"    }
  #if(margins[2]=="LN")  {fam <- "log-normal"      ; lin <- "log"         }
  #if(margins[2]=="WEI") {fam <- "Weibull"         ; lin <- "log"         }
  #if(margins[2]=="iG")  {fam <- "inverse Gaussian"; lin <- "log"         }   
  #
  #gam2$family[[1]] <- fam 
  #gam2$family[[2]] <- lin
  
  
  gam2$coefficients <- gamlss$fit$argument[1:X2.d2]
  gam2$Vp <- gamlss$magpp$Vb[1:X2.d2,1:X2.d2]
  gam2$sig2 <- 1
  gam2$edf <- gamlss$magpp$edf[1:X2.d2]
                                          
                                          
      if(length(formula) > 2){

          gam3$coefficients <- gamlss$fit$argument[X2.d2 + (1:X3.d2)]
  	  gam3$Vp <- gamlss$magpp$Vb[X2.d2 + (1:X3.d2), X2.d2 + (1:X3.d2)]
  	  gam3$sig2 <- 1
          gam3$edf <- gamlss$magpp$edf[X2.d2 + (1:X3.d2)]  
                                                                            
                           } 
  
  }  
  
  
  
  SemiParFit <- SemiParBIVProbit.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 1e+4, iterlimsp = iterlimsp, pr.tolsp = pr.tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag, naive = FALSE) 
                         
                        
  
  if(margins[2]=="probit") n <- sum(as.numeric(SemiParFit$fit$good))  
  if(Model=="BSS") n.sel <- sum(as.numeric(inde[SemiParFit$fit$good]))
  
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- SemiParBIVProbit.fit.post(SemiParFit = SemiParFit, formula.eq2 = formula.eq2, data = data, 
                                            Model = Model, VC = VC, 
                                            qu.mag = qu.mag, gam1 = gam1, gam2 = gam2, gam3 = gam3,
                                            gam4 = gam4, gam5 = gam5, gam6 = gam6)
                                            
  SemiParFit <- SemiParFit.p$SemiParFit # useful for SS models, eta2 calculatons etc.
 
  y2.m <- y2  
  if(margins[2] == "LN")  y2.m <- exp(y2)

  ##########################################################################################################################

rm(data)
if(gc.l == TRUE) gc()


L <- list(fit = SemiParFit$fit, 
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6,  
          coefficients = SemiParFit$fit$argument, 
          weights = weights, 
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,
          OR = SemiParFit.p$OR, GM = SemiParFit.p$GM,    
          n = n, n.sel = n.sel, 
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4 = X4, X5 = X5, X6 = X6, X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2,           
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, edf11 = SemiParFit.p$edf11,   
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,
          p11 = SemiParFit$fit$p11, p10 = SemiParFit$fit$p10, p01 = SemiParFit$fit$p01, p00 = SemiParFit$fit$p00, 
          p1 = SemiParFit$fit$p1, p2 = SemiParFit$fit$p2,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, etad = SemiParFit$fit$etad,
          etas = SemiParFit$fit$etas,
          y1 = y1, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, hess = hess, 
          good = SemiParFit$fit$good,
          respvec = respvec,
          qu.mag = qu.mag, sigma2 = SemiParFit.p$sigma2, sigma2.a = SemiParFit.p$sigma2.a,
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, 
          X2s = SemiParFit.p$X2s, p1n=SemiParFit.p$p1n , p2n = SemiParFit.p$p2n, 
          VC = VC, Model = Model, ig = ig, magpp = SemiParFit$magpp,
          gamlss = gamlss, gamlssfit = gamlssfit)

class(L) <- "SemiParBIVProbit"

L

}

