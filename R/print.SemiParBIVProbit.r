print.SemiParBIVProbit <- function(x,...){

  cat("\nFamily: BIVARIATE PROBIT\n\nFormula Eq. 1:\n")
  print(x$gam1$formula)

  cat("\nFormula Eq. 2:\n")
  print(x$gam2$formula)
  cat("\n")

  cat("n = ",x$n,"  rho = ",round(x$rho,3),"  total edf = ",round(x$t.edf,3),"\n\n")

invisible(x)

}
