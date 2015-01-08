print.summary.SemiParBIVProbit <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

  if(x$BivD=="N")    cop <- "BIVARIATE GAUSSIAN"
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
  
  
  main.t <- "\nERRORS' DISTRIBUTION:" 

  cat(main.t,cop)
  if(x$sel==FALSE) cat("\n\nEQUATION 1") else cat("\n\nSELECTION EQ.") 
  if(x$PL=="P") cat("\nLink function: probit\n")

  if(x$PL=="PP" || x$PL=="RPP"){ 
    if(round(x$xi1,3)==1) cat("\nLink function: probit\n") else{ 
  	if(x$PL=="PP") cat("\nLink function: power probit, ",format(x$xi1,digits=3),"\n", sep="") 
  	if(x$PL=="RPP") cat("\nLink function: reciprocal power probit, ",format(x$xi1,digits=3),"\n", sep="") 
                                                          }
  }
  
    if(x$PL=="SN"){ 
      if(round(x$xi1,3)==0) cat("\nLink function: probit\n") else 
      cat("\nLink function: skew normal, ",format(x$xi1,digits=3),"\n", sep="") 

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


  if(x$PL=="PP" || x$PL=="RPP"){ 
    if(round(x$xi2,3)==1) cat("\nLink function: probit\n") else{  

  if(x$PL=="PP") cat("\nLink function: power probit, ",format(x$xi2,digits=3),"\n", sep="") 
  if(x$PL=="RPP") cat("\nLink function: reciprocal power probit, ",format(x$xi2,digits=3),"\n", sep="") }
}

    if(x$PL=="SN"){ 
      if(round(x$xi2,3)==0) cat("\nLink function: probit\n") else 
      cat("\nLink function: skew normal, ",format(x$xi2,digits=3),"\n", sep="") 

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
    

  if(x$BivD %in% c("N","T")) {cp <- "  rho = "; as.p <- x$rho} else{ cp <- "  theta = "; as.p <- x$theta}


  #if(x$Model=="B" )  cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"\nMR = ",format(x$MR,digits=3),"%","  QPS1 = ",format(x$QPS1,digits=3),"  QPS2 = ",format(x$QPS2,digits=3),"\nCR1 = ",format(x$CR1,digits=3),"%  CR2 = ",format(x$CR2,digits=3),"%\n\n", sep="")     
  #if(x$Model=="BSS") cat("\nn = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","\ntotal edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 
  #if(x$Model=="BPO") cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

       
  
  #if(x$Model=="B" ) cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  Kendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","\ntotal edf = ",format(x$t.edf,digits=3),"  MR = ",format(x$MR,digits=3),"%","  QPS1 = ",format(x$QPS1,digits=3),"  QPS2 = ",format(x$QPS2,digits=3),"\nCR1 = ",format(x$CR1,digits=3),"%  CR2 = ",format(x$CR2,digits=3),"%\n\n", sep="")  
     
  #if(x$Model=="BSS") cat("\nn = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","\nKendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

  #if(x$Model=="BPO") cat("\nn = ",x$n,cp,format(as.p,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","\nKendall's Tau = ",round(x$KeT,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 
     
     
     
     
  nodi <- 3
  
  if(x$Model=="B" ) cat("\nn = ",x$n,cp,format(as.p,digits=nodi),"(",format(x$CIrs[1],digits=nodi),",",format(x$CIrs[2],digits=nodi),")","  Gamma measure = ",format(x$GM,digits=nodi),"(",format(x$CIgm[1],digits=nodi),",",format(x$CIgm[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"  MR = ",format(x$MR,digits=nodi),"%","  QPS1 = ",format(x$QPS1,digits=nodi),"  QPS2 = ",format(x$QPS2,digits=nodi),"\nCR1 = ",format(x$CR1,digits=nodi),"%  CR2 = ",format(x$CR2,digits=nodi),"%\n\n", sep="")  
     
  if(x$Model=="BSS") cat("\nn = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p,digits=nodi),"(",format(x$CIrs[1],digits=nodi),",",format(x$CIrs[2],digits=nodi),")","\nGamma measure = ",format(x$GM,digits=nodi),"(",format(x$CIgm[1],digits=nodi),",",format(x$CIgm[2],digits=nodi),")","  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="") 

  if(x$Model=="BPO") cat("\nn = ",x$n,cp,format(as.p,digits=nodi),"(",format(x$CIrs[1],digits=nodi),",",format(x$CIrs[2],digits=nodi),")","  Gamma measure = ",format(x$GM,digits=nodi),"(",format(x$CIgm[1],digits=nodi),",",format(x$CIgm[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="") 
            
       
       
invisible(x)
                
}

