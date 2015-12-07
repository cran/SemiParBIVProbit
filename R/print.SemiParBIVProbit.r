print.SemiParBIVProbit <- function(x, ...){

  if(x$BivD=="N" && x$Model!="BPO0")  {cop <- "BIVARIATE GAUSSIAN"                              ;lind <- "atanh(theta)"}
  if(x$BivD=="N" && x$Model=="BPO0")  {cop <- "INDEPENDENT GAUSSIAN"                              ;lind <- "atanh(theta)"}  
  if(x$BivD=="F")    {cop <- "BIVARIATE FRANK COPULA"                          ;lind <- "theta"}
  if(x$BivD=="C0")   {cop <- "BIVARIATE CLAYTON COPULA"                        ;lind <- "log(theta)"}
  if(x$BivD=="C90")  {cop <- "BIVARIATE 90\u00B0 CLAYTON COPULA"   ;lind <- "log(-theta)"}
  if(x$BivD=="C180") {cop <- "BIVARIATE 180\u00B0 CLAYTON COPULA"               ;lind <- "log(theta)"}
  if(x$BivD=="C270") {cop <- "BIVARIATE 270\u00B0 CLAYTON COPULA"  ;lind <- "log(-theta)"}  
  if(x$BivD=="J0")   {cop <- "BIVARIATE JOE COPULA"                            ;lind <- "log(theta-1)"}
  if(x$BivD=="J90")  {cop <- "BIVARIATE 90\u00B0 JOE COPULA"       ;lind <- "log(-theta-1)"}
  if(x$BivD=="J180") {cop <- "BIVARIATE 180\u00B0 JOE COPULA"                   ;lind <- "log(theta-1)"}
  if(x$BivD=="J270") {cop <- "BIVARIATE 270\u00B0 JOE COPULA"      ;lind <- "log(-theta-1)"}
  if(x$BivD=="G0")   {cop <- "BIVARIATE GUMBEL COPULA"                         ;lind <- "log(theta-1)"}
  if(x$BivD=="G90")  {cop <- "BIVARIATE 90\u00B0 GUMBEL COPULA"    ;lind <- "log(-theta-1)"}
  if(x$BivD=="G180") {cop <- "BIVARIATE 180\u00B0 GUMBEL COPULA"                ;lind <- "log(theta-1)"}
  if(x$BivD=="G270") {cop <- "BIVARIATE 270\u00B0 GUMBEL COPULA"   ;lind <- "log(-theta-1)"} 
    
  cp <- "  theta = "; as.p <- x$theta.a
  
  if(x$margins[2]=="probit")                     m2l <- "probit"
  if(x$margins[2] %in% c("N","GU","rGU","LO") )  m2l <- "identity"
  if(x$margins[2] %in% c("LN","WEI","WEI2","iG","GA","iGA","DAGUM") ) m2l <- "log" 
  

  
  cat("\nERRORS' DISTRIBUTION:",cop)

  cat("\n\nEQUATION 1")
  cat("\nFamily: Bernoulli") 
  cat("\nLink function: probit")
  cat("\nFormula: "); print(x$gam1$formula)

  cat("\nEQUATION 2")
  
  if(x$margins[2]=="probit") cat("\nFamily: Bernoulli") 
  if(x$margins[2]=="N")      cat("\nFamily: Gaussian")  
  if(x$margins[2]=="GU")     cat("\nFamily: Gumbel")    
  if(x$margins[2]=="rGU")    cat("\nFamily: reverse Gumbel")  
  if(x$margins[2]=="LO")     cat("\nFamily: logistic")   
  if(x$margins[2]=="LN")     cat("\nFamily: log-normal") 
  if(x$margins[2]=="WEI")    cat("\nFamily: Weibull") 
  if(x$margins[2]=="WEI2")   cat("\nFamily: Weibull (type 2)")   
  if(x$margins[2]=="iG")     cat("\nFamily: inverse Gaussian") 
  if(x$margins[2]=="GA")     cat("\nFamily: gamma")    
  if(x$margins[2]=="iGA")    cat("\nFamily: inverse gamma")    
  if(x$margins[2]=="DAGUM")  cat("\nFamily: DAGUM")  
  #if(x$margins[2]=="ZAGA")   cat("\nFamily: Zero adjusted gamma")  
  
  
 
  cat("\nLink function:",m2l,"\n")
  cat("Formula: "); print(x$gam2$formula)
  
  
  
  if(!is.null(x$X3) && is.null(x$X4) ){
  
  cat("\nEQUATION 3")
  cat("\nLink function:",lind,"\n") 
  cat("Formula: "); print(x$gam3$formula)
  
  }
  
  if(!is.null(x$X3) && !is.null(x$X4) &&  is.null(x$X5)){
  
  cat("\nEQUATION 3")
  cat("\nLink function:","log(sigma^2)","\n") 
  cat("Formula: "); print(x$gam3$formula)  
  
  cat("\nEQUATION 4")
  cat("\nLink function:",lind,"\n") 
  cat("Formula: "); print(x$gam4$formula)
  
  } 
  
  if(!is.null(x$X3) && !is.null(x$X4) && !is.null(x$X5)){
  
  cat("\nEQUATION 3")
  cat("\nLink function:","log(sigma^2)","\n") 
  cat("Formula: "); print(x$gam3$formula)  
  
  
  cat("\nEQUATION 4")
  if(x$margins[2]=="DAGUM") cat("\nLink function:","log(nu)","\n")  
  #if(x$margins[2]=="ZAGA")  cat("\nLink function:","logit(nu)","\n")  
  cat("Formula: "); print(x$gam4$formula)  
    

  cat("\nEQUATION 5")
  cat("\nLink function:",lind,"\n") 
  cat("Formula: "); print(x$gam5$formula)
  
  }  
  
  
  
  
  cont2par <- c("N","GU","rGU","LO","LN","WEI","WEI2","iG","GA","iGA")  
  cont3par <- c("DAGUM")  
  
  
  cat("\n")
            
  if(x$Model %in% c("B","BPO") && x$margins[2]=="probit") cat("n = ",x$n,cp,format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$Model == "BPO0")                                   cat("n = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

  if(x$Model=="BSS")                                      cat("n = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  
  if(x$Model=="B" && x$margins[2] %in% cont2par ) cat("n = ",x$n,"  sigma^2 = ",x$sigma2.a, cp, format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

  if(x$Model=="B" && x$margins[2] %in% cont3par ) cat("n = ",x$n,"  sigma^2 = ",x$sigma2.a, "  nu = ",x$nu.a, "\ntheta = ", format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")




invisible(x)

}

