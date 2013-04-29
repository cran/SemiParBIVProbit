print.summary.SemiParBIVProbit <- function(x,digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

  re <- ""
  if(x$RE==TRUE){ if(x$RE.type=="NP") re <- "with NP Random Effects (RE)"; if(x$RE.type=="N") re <- "with N Random Effects"}
  if(x$sel==FALSE) cat("\nFamily: BIVARIATE PROBIT",re,"\n\nEQUATION 1: ") else cat("\nFamily: BIVARIATE PROBIT\n\nSELECTION EQ.: ") 
  print(x$formula1)
  cat("\n") 
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){
    cat("Smooth terms' approximate significance:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

  if(x$sel==FALSE) cat("\nEQUATION 2: ") else cat("\nOUTCOME EQ.: ")
  print(x$formula2)
  cat("\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Smooth terms' approximate significance:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
  if(x$RE==TRUE){
  cat("\nEstimated parameters of RE distribution: ")
  cat("\n\n")
  printCoefmat(x$table.RE,digits = digits, signif.stars = signif.stars, na.print = "NA",...)
  cat("\n")
  }
  
  if(x$sel==FALSE && x$RE==FALSE) cat("\nn = ",x$n,"  rho = ",format(x$rho,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"  MR = ",format(x$MR,digits=3),"%\n","QPS1 = ",format(x$QPS1,digits=3),"  QPS2 = ",format(x$QPS2,digits=3),"  CR1 = ",format(x$CR1,digits=3),"%  CR2 = ",format(x$CR2,digits=3),"%\n\n", sep="") 

  if(x$sel==TRUE && x$RE==FALSE) cat("\nn = ",x$n,"  n.sel = ",x$n.sel,"  rho = ",format(x$rho,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

  if(x$RE==TRUE) cat("\nn = ",x$n,"  rho = ",format(x$rho,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") 

                       
invisible(x)
                
}














