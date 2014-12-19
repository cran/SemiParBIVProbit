print.est.prev <- function(x, ...){

es <- format(x$res*100, digits = 3, trim=TRUE)

cat("\nEstimated prevalence (%) with corresponding ",(1-x$sig.lev)*100,"% Confidence Intervals:\n\n",sep="")

cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

invisible(x)

}



