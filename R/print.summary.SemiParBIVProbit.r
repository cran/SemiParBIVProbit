print.summary.SemiParBIVProbit <- function(x,digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

  if(x$BivD=="N")    cop <- "BIVARIATE PROBIT"
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

  re <- ""
  if(x$RE==TRUE && x$BivD=="N"){ if(x$RE.type=="NP") re <- "with NP Random Effects (RE)"; if(x$RE.type=="N") re <- "with N Random Effects"}
  if(x$sel==FALSE) cat("\nFamily:",cop,re,"\n\nEQUATION 1: ") else cat("\nFamily",cop,"\n\nSELECTION EQ.: ") 
  print(x$formula1)
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

  pP1 <- length(x$pPen1)
  pP2 <- length(x$pPen2)
  
    if(x$l.sc1!=0){ #|| (pP1!=0 && x$l.sp1>1)
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

  if(x$sel==FALSE) cat("\nEQUATION 2: ") else cat("\nOUTCOME EQ.: ")
  print(x$formula2)
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sc2!=0){
    cat("Smooth components' approximate significance:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
  if(x$RE==TRUE){
  cat("\nEstimated parameters of RE distribution:\n")
  printCoefmat(x$table.RE,digits = digits, signif.stars = signif.stars, na.print = "NA",...)
  cat("\n")
  }

  if(x$BivD %in% c("N","T")) {cp <- "  rho = "; as.p <- x$rho} else{ cp <- "  theta = "; as.p <- x$theta}


  if(x$sel==FALSE && x$RE==FALSE) cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  Kendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","\ntotal edf = ",format(x$t.edf,digits=3),"  MR = ",format(x$MR,digits=3),"%","  QPS1 = ",format(x$QPS1,digits=3),"  QPS2 = ",format(x$QPS2,digits=3),"\nCR1 = ",format(x$CR1,digits=3),"%  CR2 = ",format(x$CR2,digits=3),"%\n\n", sep="")  
     
  if(x$sel==TRUE  && x$RE==FALSE) cat("\nn = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","\nKendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

  if(x$RE==TRUE) cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","\nKendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

                       
invisible(x)
                
}

