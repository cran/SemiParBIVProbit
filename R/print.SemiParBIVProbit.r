print.SemiParBIVProbit <- function(x, ...){

  if(x$BivD=="N")    {cop <- "BIVARIATE GAUSSIAN"                              ;lind <- "atanh(rho)"}
  if(x$BivD=="F")    {cop <- "BIVARIATE FRANK COPULA"                          ;lind <- "theta"}
  if(x$BivD=="C0")   {cop <- "BIVARIATE CLAYTON COPULA"                        ;lind <- "log(theta)"}
  if(x$BivD=="C90")  {cop <- "BIVARIATE ROTATED CLAYTON COPULA (90 DEGREES)"   ;lind <- "log(-theta)"}
  if(x$BivD=="C180") {cop <- "BIVARIATE SURVIVAL CLAYTON COPULA"               ;lind <- "log(theta)"}
  if(x$BivD=="C270") {cop <- "BIVARIATE ROTATED CLAYTON COPULA (270 DEGREES)"  ;lind <- "log(-theta)"}  
  if(x$BivD=="J0")   {cop <- "BIVARIATE JOE COPULA"                            ;lind <- "log(theta-1)"}
  if(x$BivD=="J90")  {cop <- "BIVARIATE ROTATED JOE COPULA (90 DEGREES)"       ;lind <- "log(-theta-1)"}
  if(x$BivD=="J180") {cop <- "BIVARIATE SURVIVAL JOE COPULA"                   ;lind <- "log(theta-1)"}
  if(x$BivD=="J270") {cop <- "BIVARIATE ROTATED JOE COPULA (270 DEGREES)"      ;lind <- "log(-theta-1)"}
  if(x$BivD=="G0")   {cop <- "BIVARIATE GUMBEL COPULA"                         ;lind <- "log(theta-1)"}
  if(x$BivD=="G90")  {cop <- "BIVARIATE ROTATED GUMBEL COPULA (90 DEGREES)"    ;lind <- "log(-theta-1)"}
  if(x$BivD=="G180") {cop <- "BIVARIATE SURVIVAL GUMBEL COPULA"                ;lind <- "log(theta-1)"}
  if(x$BivD=="G270") {cop <- "BIVARIATE ROTATED GUMBEL COPULA (270 DEGREES)"   ;lind <- "log(-theta-1)"} 
  
  if(x$BivD %in% c("N") ) {cp <- "  rho = "; as.p <- x$rho.a} else{ cp <- "  theta = "; as.p <- x$theta.a}

  cat("\nERRORS' DISTRIBUTION:",cop)

  cat("\n\nEQUATION 1")
  if(x$PL=="P") cat("\nLink function: probit\n")

  if(x$PL=="PP" || x$PL=="RPP"){ 
    if(round(x$xi1,3)==1) cat("\nLink function: probit\n") else{ 
  if(x$PL=="PP") cat("\nLink function: power probit\n") 
  if(x$PL=="RPP") cat("\nLink function: reciprocal power probit\n") }
  }
  
    if(x$PL=="SN"){ 
      if(round(x$xi1,3)==0) cat("\nLink function: probit\n") else 
    cat("\nLink function: skew normal\n") 
  }


  print(x$gam1$formula)

  cat("\nEQUATION 2")
  if(x$PL=="P") cat("\nLink function: probit\n")
  
  if(x$PL=="PP" || x$PL=="RPP"){ 
    if(round(x$xi2,3)==1) cat("\nLink function: probit\n") else{  
  if(x$PL=="PP") cat("\nLink function: power probit\n") 
  if(x$PL=="RPP") cat("\nLink function: reciprocal power probit\n") }
}

    if(x$PL=="SN"){ 
      if(round(x$xi2,3)==0) cat("\nLink function: probit\n") else 
    cat("\nLink function: skew normal\n") 
  }

  print(x$gam2$formula)
  
  
  
  if(!is.null(x$X3)){
  
  cat("\nEQUATION 3")
  cat("\nLink function:",lind,"\n") 
  
  f <- x$gam3$formula; environment(f) <- environment(x$gam2$formula) 
  print(f)
  
  }
  
  
  cat("\n")
            
  if(x$sel==FALSE) cat("n = ",x$n,cp,format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  
  if(x$sel==TRUE)  cat("n = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  
  

invisible(x)

}

