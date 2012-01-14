print.SemiParBIVProbit <- function(x,...){

  if(x$sel==FALSE){ 
  re <- ""
  if(x$npRE==TRUE) re <- "with Random Effects"
  cat("\nFamily: BIVARIATE PROBIT",re,"\n\nFormula Eq. 1: ")
  print(x$gam1$formula)

  cat("Formula Eq. 2: ")
  print(x$gam2$formula)
  cat("\n")

  cat("n = ",x$n,"  rho = ",format(x$rho, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n")

  }else{

  cat("\nFamily: BIVARIATE PROBIT\n\nSELECTION EQ.: ")
  print(x$gam1$formula)

  cat("OUTCOME   EQ.: ")
  print(x$gam2$formula)
  cat("\n")

  cat("n = ",x$n,"  n.sel = ",x$n.sel,"  rho = ",format(x$rho, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n")
  }


invisible(x)

}
