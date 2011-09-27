print.summary.SemiParBIVProbit <- function(x,digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

  if(x$sel==FALSE) cat("\nFamily: BIVARIATE PROBIT\n\nEQUATION 1: ") else cat("\nFamily: BIVARIATE PROBIT\n\nSELECTION EQ.: ") 
  print(x$formula1)
  cat("\n") 
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0 && x$l.sp2!=0){
    cat("Smooth terms' approximate significance:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

  if(x$sel==FALSE) cat("\nEQUATION 2: ") else cat("\nOUTCOME EQ.: ")
  print(x$formula2)
  cat("\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",)
  cat("\n")

    if(x$l.sp1!=0 && x$l.sp2!=0){
    cat("Smooth terms' approximate significance:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

  if(x$sel==FALSE){ 
    if(x$l.sp1!=0 && x$l.sp2!=0){ cat("n = ",x$n,"  rho = ",format(x$rho,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"  MR = ",format(x$MR,digits=3),"%\n","QPS1 = ",format(x$QPS1,digits=3),"  QPS2 = ",format(x$QPS2,digits=3),"  CR1 = ",format(x$CR1,digits=3),"%  CR2 = ",format(x$CR2,digits=3),"%\n\n", sep="") }
    else{ cat("n = ",x$n,"  rho = ",format(x$rho,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"  MR = ",format(x$MR,digits=3),"%\n","QPS1 = ",format(x$QPS1,digits=3),"  QPS2 = ",format(x$QPS2,digits=3),"  CR1 = ",format(x$CR1,digits=3),"%  CR2 = ",format(x$CR2,digits=3),"%\n\n",sep="")  }
                  }else{
    if(x$l.sp1!=0 && x$l.sp2!=0){ cat("n = ",x$n,"  n.sel = ",x$n.sel,"  rho = ",format(x$rho,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n", sep="") }
    else{ cat("n = ",x$n,"  n.sel = ",x$n.sel,"  rho = ",format(x$rho,digits=3),"(",format(x$CIrs[1],digits=3),",",format(x$CIrs[2],digits=3),")","  total edf = ",format(x$t.edf,digits=3),"\n\n",sep="")  }
                       }

                        

invisible(x)
                
}
















