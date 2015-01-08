print.mb <- function(x, ...){

if(x$Model == "B") nf <- nf1 <- "ATE"
if(x$Model == "BSS"){ nf <- "prevalence"; nf1 <- "Prevalence"}

cat("\nWorst-Case Manski's bounds for ",nf," (%):\n","Lower Bound: ",x$LB,"\n","Upper Bound: ",x$UB,"\n",sep="")
cat("\nImbens&Manski's CI for ",nf," (%): [",x$CI[1],",",x$CI[2],"]","\n",sep="")

cat("\n",nf1," assuming random assignment (%): ",x$av.p,"\n\n",sep="")

invisible(x)

}



