print.summary.SemiParBIVProbit <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

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
  
  
  main.t <- "\nERRORS' DISTRIBUTION:" 
  cp <- "  theta = "; as.p <- x$theta

  if(x$margins[2]=="probit")                     m2l <- "probit"
  if(x$margins[2] %in% c("N","GU","rGU","LO") )  m2l <- "identity"
  if(x$margins[2] %in% c("LN","WEI","WEI2","iG","GA","iGA","DAGUM") ) m2l <- "log"   
  
  
  cat(main.t,cop) 
  cat("\n\nEQUATION 1")  
  cat("\nFamily: Bernoulli") 
  cat("\nLink function: probit\n")
  cat("Formula: "); print(x$formula1)
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){ 
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
    
  cat("\nEQUATION 2")
  
  if(x$margins[2]=="probit") cat("\nFamily: Bernoulli") 
  if(x$margins[2]=="N")      cat("\nFamily: Gaussian")  
  if(x$margins[2]=="GU")     cat("\nFamily: Gumbel")    
  if(x$margins[2]=="rGU")    cat("\nFamily: reverse Gumbel")  
  if(x$margins[2]=="LO")     cat("\nFamily: logistic")   
  if(x$margins[2]=="LN")     cat("\nFamily: log-normal") 
  if(x$margins[2]=="WEI")    cat("\nFamily: Weibull")  
  if(x$margins[2]=="WEI2")    cat("\nFamily: Weibull (type 2)")   
  if(x$margins[2]=="iG")     cat("\nFamily: inverse Gaussian")    
  if(x$margins[2]=="GA")     cat("\nFamily: gamma")    
  if(x$margins[2]=="iGA")    cat("\nFamily: inverse gamma")    
  if(x$margins[2]=="DAGUM")  cat("\nFamily: DAGUM")    
  
  
  cat("\nLink function:",m2l,"\n")
  cat("Formula: "); print(x$formula2)    
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    


if(!is.null(x$tableP3) && is.null(x$tableP4)  ){

  cat("\nEQUATION 3")
  cat("\nLink function:",lind,"\n") 
  cat("Formula: "); print(x$formula3)
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 


if(!is.null(x$tableP3) && !is.null(x$tableP4) && is.null(x$tableP5)  ){


  cat("\nEQUATION 3")
  cat("\nLink function:","log(sigma^2)","\n") 
  cat("Formula: "); print(x$formula3)
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }  


  cat("\nEQUATION 4")
  cat("\nLink function:",lind,"\n") 
  cat("Formula: "); print(x$formula4)
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 


if(!is.null(x$tableP3) && !is.null(x$tableP4) && !is.null(x$tableP5)  ){


  cat("\nEQUATION 3")
  cat("\nLink function:","log(sigma^2)","\n") 
  cat("Formula: "); print(x$formula3)
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }  
    
  cat("\nEQUATION 4")
  cat("\nLink function:","log(nu)","\n") 
  cat("Formula: "); print(x$formula4)
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }     


  cat("\nEQUATION 5")
  cat("\nLink function:",lind,"\n") 
  cat("Formula: "); print(x$formula5)
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp5!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 





  cont2par <- c("N","GU","rGU","LO","LN","WEI","WEI2","iG","GA","iGA")  
  cont3par <- c("DAGUM")  
  


  CIrs <- colMeans(x$CItheta, na.rm = TRUE)
  if(x$margins[2] %in% cont2par)  CIsig2 <- colMeans(x$CIsig2, na.rm = TRUE)
  if(x$margins[2] %in% cont3par) {CIsig2 <- colMeans(x$CIsig2, na.rm = TRUE) ; CInu <- colMeans(x$CInu, na.rm = TRUE)}




  nodi <- 3
  
  if( (x$Model=="B" || x$Model=="BPO") && x$margins[2]=="probit") cat("\nn = ",x$n,cp,format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

  if( x$Model == "BPO0" ) cat("\nn = ",x$n,"  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
    
    
  if(x$Model=="BSS") cat("\nn = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="") 
     
  if(x$Model=="B" && x$margins[2] %in% cont2par ) cat("\nn = ",x$n,cp,format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","  sigma^2 = ",format(x$sigma2,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       
  if(x$Model=="B" && x$margins[2] %in% cont3par ) cat("\nn = ",x$n,cp,format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","  sigma^2 = ",format(x$sigma2,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")","\nnu = ",format(x$nu,digits = nodi),"(",format(CInu[1],digits=nodi),",",format(CInu[2],digits=nodi),")","  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       
       
invisible(x)
                
}





















