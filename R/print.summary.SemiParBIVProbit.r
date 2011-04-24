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
    if(x$l.sp1!=0 && x$l.sp2!=0){ cat("n = ",x$n,"  rho = ",round(x$rho,3),"(",round(x$CIrs[1],3),",",round(x$CIrs[2],3),")","  total edf = ",round(x$t.edf,3),"  MR = ",round(x$MR,3),"%\n","QPS1 = ",round(x$QPS1,3),"  QPS2 = ",round(x$QPS2,3),"  CR1 = ",round(x$CR1,3),"%  CR2 = ",round(x$CR2,3),"%\n\n", sep="") }
    else{ cat("n = ",x$n,"  rho = ",round(x$rho,3),"(",round(x$CIrs[1],3),",",round(x$CIrs[2],3),")","  total edf = ",round(x$t.edf,3),"  MR = ",round(x$MR,3),"%\n","QPS1 = ",round(x$QPS1,3),"  QPS2 = ",round(x$QPS2,3),"  CR1 = ",round(x$CR1,3),"%  CR2 = ",round(x$CR2,3),"%\n\n",sep="")  }
                  }else{
    if(x$l.sp1!=0 && x$l.sp2!=0){ cat("n = ",x$n,"  rho = ",round(x$rho,3),"(",round(x$CIrs[1],3),",",round(x$CIrs[2],3),")","  total edf = ",round(x$t.edf,3),"\n\n", sep="") }
    else{ cat("n = ",x$n,"  rho = ",round(x$rho,3),"(",round(x$CIrs[1],3),",",round(x$CIrs[2],3),")","  total edf = ",round(x$t.edf,3),"\n\n",sep="")  }
                       }

                        

invisible(x)
                
}
















