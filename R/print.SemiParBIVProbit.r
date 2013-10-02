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

  if(x$BivD %in% c("N","T") ) {cp <- "  rho = "; as.p <- x$rho} else{ cp <- "  theta = "; as.p <- x$theta}

  if(x$RE==TRUE){ if(x$RE.type=="NP") re <- "Bivariate Nonparametric Random Effects in linear predictors"; if(x$RE.type=="N") re <- "Bivariate Normal Random Effects in linear predictors"}

  cat("\nERRORS' DISTRIBUTION:",cop)
  if(x$RE==TRUE) cat("\n",re,sep ="")  

  cat("\n\nEQUATION 1")
  if(x$gev.eq1==FALSE) cat("\nLink function: probit\n") else cat("\nLink function: gev\n")
  print(x$gam1$formula)

  cat("\nEQUATION 2")
  if(x$gev.eq2==FALSE) cat("\nLink function: probit\n") else cat("\nLink function: gev\n")
  print(x$gam2$formula)
  cat("\n")

  if(x$sel==FALSE) cat("n = ",x$n,cp,format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="") else cat("n = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

invisible(x)

}


















