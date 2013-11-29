print.SemiParBIVProbit <- function(x,...){

  if(x$BivD=="N")    cop <- "BIVARIATE NORMAL"
  if(x$BivD=="F")    cop <- "BIVARIATE FRANK COPULA"
  if(x$BivD=="T")    cop <- paste("BIVARIATE STUDENT-T COPULA (",x$nu," DEGREES OF FREEDOM)",sep="")
  if(x$BivD=="C0")   cop <- "BIVARIATE CLAYTON COPULA"
  if(x$BivD=="C90")  cop <- "BIVARIATE ROTATED CLAYTON COPULA (90 DEGREES)"
  if(x$BivD=="C180") cop <- "BIVARIATE SURVIVAL CLAYTON COPULA"
  if(x$BivD=="C270") cop <- "BIVARIATE ROTATED CLAYTON COPULA (270 DEGREES)"  
  if(x$BivD=="J0")   cop <- "BIVARIATE JOE COPULA"
  if(x$BivD=="J90")  cop <- "BIVARIATE ROTATED JOE COPULA (90 DEGREES)"
  if(x$BivD=="J180") cop <- "BIVARIATE SURVIVAL JOE COPULA"
  if(x$BivD=="J270") cop <- "BIVARIATE ROTATED JOE COPULA (270 DEGREES)"
  if(x$BivD=="G0")   cop <- "BIVARIATE GUMBEL COPULA"
  if(x$BivD=="G90")  cop <- "BIVARIATE ROTATED GUMBEL COPULA (90 DEGREES)"
  if(x$BivD=="G180") cop <- "BIVARIATE SURVIVAL GUMBEL COPULA"
  if(x$BivD=="G270") cop <- "BIVARIATE ROTATED GUMBEL COPULA (270 DEGREES)"   
  
  if(x$BivD=="BB1.0")   cop <- "BIVARIATE CLAYTON-GUMBEL COPULA" 
  if(x$BivD=="BB1.90")  cop <- "BIVARIATE ROTATED CLAYTON-GUMBEL COPULA (90 DEGREES)"  
  if(x$BivD=="BB1.180") cop <- "BIVARIATE ROTATED CLAYTON-GUMBEL COPULA (180 DEGREES)"  
  if(x$BivD=="BB1.270") cop <- "BIVARIATE ROTATED CLAYTON-GUMBEL COPULA (270 DEGREES)"  
  
  if(x$BivD=="BB6.0")   cop <- "BIVARIATE JOE-GUMBEL COPULA" 
  if(x$BivD=="BB6.90")  cop <- "BIVARIATE ROTATED JOE-GUMBEL COPULA (90 DEGREES)"  
  if(x$BivD=="BB6.180") cop <- "BIVARIATE ROTATED JOE-GUMBEL COPULA (180 DEGREES)"  
  if(x$BivD=="BB6.270") cop <- "BIVARIATE ROTATED JOE-GUMBEL COPULA (270 DEGREES)"  
  
  if(x$BivD=="BB7.0")   cop <- "BIVARIATE JOE-GUMBEL COPULA" 
  if(x$BivD=="BB7.90")  cop <- "BIVARIATE ROTATED JOE-CLAYTON COPULA (90 DEGREES)"  
  if(x$BivD=="BB7.180") cop <- "BIVARIATE ROTATED JOE-CLAYTON COPULA (180 DEGREES)"  
  if(x$BivD=="BB7.270") cop <- "BIVARIATE ROTATED JOE-CLAYTON COPULA (270 DEGREES)"  
  
  if(x$BivD=="BB8.0")   cop <- "BIVARIATE JOE-FRANK COPULA" 
  if(x$BivD=="BB8.90")  cop <- "BIVARIATE ROTATED JOE-FRANK COPULA (90 DEGREES)"  
  if(x$BivD=="BB8.180") cop <- "BIVARIATE ROTATED JOE-FRANK COPULA (180 DEGREES)"  
  if(x$BivD=="BB8.270") cop <- "BIVARIATE ROTATED JOE-FRANK COPULA (270 DEGREES)"  

  if(x$BivD %in% c("N","T") ) {cp <- "  rho = "; as.p <- x$rho} else{ cp <- "  theta = "; as.p <- x$theta}

  if(x$RE==TRUE){ if(x$RE.type=="NP") re <- "Bivariate Nonparametric Random Effects in linear predictors"; if(x$RE.type=="N") re <- "Bivariate Normal Random Effects in linear predictors"}

  cat("\nERRORS' DISTRIBUTION:",cop)
  if(x$RE==TRUE) cat("\n",re,sep ="")  

  cat("\n\nEQUATION 1")
  if(x$PL=="P") cat("\nLink function: probit\n")

  if(x$PL!="P"){ 
    if(x$xi1==1) cat("\nLink function: probit\n") else{ 
  if(x$PL=="PP") cat("\nLink function: power probit\n") 
  if(x$PL=="RPP") cat("\nLink function: reciprocal power probit\n") }
  }


  print(x$gam1$formula)

  cat("\nEQUATION 2")
  if(x$PL=="P") cat("\nLink function: probit\n")
  if(x$PL!="P"){ 
    if(x$xi2==1) cat("\nLink function: probit\n") else{  
  if(x$PL=="PP") cat("\nLink function: power probit\n") 
  if(x$PL=="RPP") cat("\nLink function: reciprocal power probit\n") }
}
  print(x$gam2$formula)
  cat("\n")

  del <- ""

  if(x$BivD %in% c("BB1.0","BB1.180","BB1.90","BB1.270",
                   "BB6.0","BB6.180","BB6.90","BB6.270",
                   "BB7.0","BB7.180","BB7.90","BB7.270",
                   "BB8.0","BB8.180","BB8.90","BB8.270") ) del <- paste("  delta = ",format(x$delta, digits=3)) 
                   
  if(x$sel==FALSE) cat("n = ",x$n,cp,format(as.p, digits=3),del,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  
  if(x$sel==TRUE)  cat("n = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p, digits=3),del,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  
  

invisible(x)

}

