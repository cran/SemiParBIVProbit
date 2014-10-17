LM.bpm <- function(formula, data = list(), weights = NULL, subset = NULL, 
                   Model = "B", hess = FALSE, gamma = 1, pPen1 = NULL, pPen2 = NULL){

  sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- NULL  
  end <- 0
  BivD <- "N"; PL <- "P"
  fp <- FALSE

  if(Model == "BPO") stop("This test's implementation is not valid for bivariate models with partial observability.")

  ig <- interpret.gam(formula)
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- ig$fake.formula 
  mf$Model <- mf$hess <- mf$gamma <- mf$pPen1 <- mf$pPen2 <- NULL  
  mf$drop.unused.levels <- TRUE 
  if(Model=="BSS") mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
  
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




  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), gamma=gamma, weights=weights, 
                              data=data, paraPen=pPen1),list(weights=weights))) 

  X1 <- model.matrix(gam1)
  X1.d2 <- dim(X1)[2]
  l.sp1 <- length(gam1$sp)
  y1 <- gam1$y
  n <- length(y1) 


  if(Model=="B"){
  
  gam2  <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=weights, 
                           data=data, paraPen=pPen2),list(weights=weights))) 
  X2 <- model.matrix(gam2)
  X2.d2 <- dim(X2)[2]
  l.sp2 <- length(gam2$sp)
  y2 <- gam2$y 

  y1.y2 <- y1*y2
  y1.cy2 <- y1*(1-y2)
  cy1.y2 <- (1-y1)*y2
  cy1.cy2 <- (1-y1)*(1-y2)

  func.opt <- bprobgHs                       
  
  }
  

  
  if(Model=="BSS"){

  inde <- y1 > 0
  gam2 <- eval(substitute(gam(formula.eq2, binomial(link="probit"), gamma=gamma, weights=weights, 
                              data=data, subset=inde, paraPen=pPen2),list(weights=weights,inde=inde)))                              
  X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2) 
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y
  l.sp2 <- length(gam2$sp)

  cy1 <- (1-y1)
  y1.y2 <- y1*y2
  y1.cy2 <- y1*(1-y2)
  
  func.opt <- bprobgHsSS 

  }


  gp1 <- gam1$nsdf
  gp2 <- gam2$nsdf   
  
  if(l.sp1!=0 && l.sp2!=0) sp <- c(gam1$sp,gam2$sp)
  if(l.sp1==0 && l.sp2!=0) sp <- c(gam2$sp)
  if(l.sp1!=0 && l.sp2==0) sp <- c(gam1$sp)


  if( (l.sp1!=0 || l.sp2!=0) ) qu.mag <- S.m(gam1, gam2, l.sp1, l.sp2) 


  respvec <- list(y1 = y1,
                  y2 = y2,
                  y1.y2 = y1.y2, 
                  y1.cy2 = y1.cy2, 
                  cy1.y2 = cy1.y2, 
                  cy1.cy2 = cy1.cy2, 
                  cy1 = cy1)
  
  VC <- list(X1 = X1, 
             X2 = X2, 
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             gp1 = gp1, 
             gp2 = gp2,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             gamma = gamma,
             weights = weights,
             hess = hess,
             pPen1 = pPen1,
             pPen2 = pPen2,
             Model = Model,
             end = end, fp = fp,
             BivD = BivD, nC = 1, nu = 3)


params <- c(coef(gam1),coef(gam2),0)

resf <- func.opt(params, sp.xi1 = NULL, sp.xi2 = NULL, PL, eqPL = NULL, valPL = NULL, 
                 fitPL = NULL, respvec, VC, sp, qu.mag)

G   <- resf$gradient
var <- resf$hessian

var.eig <- eigen(var, symmetric=TRUE)   
if(min(var.eig$values) < .Machine$double.eps){ var <- as.matrix( nearPD( var, ensureSymmetry = FALSE )$mat )
                                               var.eig <- eigen(var, symmetric=TRUE) }
var <- var.eig$vec%*%tcrossprod(diag(1/var.eig$val),var.eig$vec)  

ev <- as.numeric(t(G)%*%var%*%G)
return(pchisq(ev,1,lower.tail=FALSE))

}






