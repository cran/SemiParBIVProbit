print.AT <- function(x, ...){

es <- format(x$res*100, digits = 3, trim=TRUE)

cat("\nAverage treatment effect (%) with ",(1-x$prob.lev)*100,"% confidence interval:\n\n",sep="")

cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

invisible(x)

}
