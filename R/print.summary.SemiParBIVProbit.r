print.summary.SemiParBIVProbit <- function(x,digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

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
  main.t <- "\nERRORS' DISTRIBUTION:" 

  if(x$RE==TRUE){ if(x$RE.type=="NP") re <- "\nBivariate Nonparametric Random Effects (RE) in linear predictors"; if(x$RE.type=="N") re <- "\nBivariate Normal Random Effects in linear predictors"}
  cat(main.t,cop)
  if(x$RE==TRUE) cat(re)  
  if(x$sel==FALSE) cat("\n\nEQUATION 1") else cat("\n\nSELECTION EQ.") 
  if(x$PL=="P") cat("\nLink function: probit\n")

  #if(x$PL=="PP") cat("\nLink function: power probit, ",format(x$xi1,digits=3),"(",format(x$CIl1[1],digits=3),",",format(x$CIl1[2],digits=3),")\n", sep="") 
  #if(x$PL=="RPP") cat("\nLink function: reciprocal power probit, ",format(x$lambda1,digits=3),"(",format(x$CIl1[1],digits=3),",",format(x$CIl1[2],digits=3),")\n",sep="") 


  if(x$PL!="P"){ 
    if(x$xi1==1) cat("\nLink function: probit\n") else{ 
  	if(x$PL=="PP") cat("\nLink function: power probit, ",format(x$xi1,digits=3),"\n", sep="") 
  	if(x$PL=="RPP") cat("\nLink function: reciprocal power probit, ",format(x$xi1,digits=3),"\n", sep="") 
                                                          }
  }


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

  if(x$sel==FALSE) cat("\nEQUATION 2") else cat("\nOUTCOME EQ.")
  if(x$PL=="P") cat("\nLink function: probit\n")

  #if(x$PL=="PP") cat("\nLink function: power probit, ",format(x$xi2,digits=3),"(",format(x$CIl2[1],digits=3),",",format(x$CIl2[2],digits=3),")\n", sep="") 
  #if(x$PL=="RPP") cat("\nLink function: reciprocal power probit, ",format(x$xi2,digits=3),"(",format(x$CIl2[1],digits=3),",",format(x$CIl2[2],digits=3),")\n",sep="") 


  if(x$PL!="P"){ 
    if(x$xi2==1) cat("\nLink function: probit\n") else{  

  if(x$PL=="PP") cat("\nLink function: power probit, ",format(x$xi2,digits=3),"\n", sep="") 
  if(x$PL=="RPP") cat("\nLink function: reciprocal power probit, ",format(x$xi2,digits=3),"\n", sep="") }
}




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



  del <- ""

  if(x$BivD %in% c("BB1.0","BB1.180","BB1.90","BB1.270",
                   "BB6.0","BB6.180","BB6.90","BB6.270",
                   "BB7.0","BB7.180","BB7.90","BB7.270",
                   "BB8.0","BB8.180","BB8.90","BB8.270") ) del <- paste("\ndelta = ",format(x$delta, digits=3),"(",format(x$CId[1],digits=3),",",format(x$CId[2],digits=3),")", sep="") 
                   


  if(x$sel==FALSE && x$RE==FALSE) cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")",del,"  Kendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","\ntotal edf = ",format(x$t.edf,digits=3),"  MR = ",format(x$MR,digits=3),"%","  QPS1 = ",format(x$QPS1,digits=3),"  QPS2 = ",format(x$QPS2,digits=3),"\nCR1 = ",format(x$CR1,digits=3),"%  CR2 = ",format(x$CR2,digits=3),"%\n\n", sep="")  
     
  if(x$sel==TRUE  && x$RE==FALSE) cat("\nn = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")",del,"\nKendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

  if(x$RE==TRUE) cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","\nKendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

                       
invisible(x)
                
}

