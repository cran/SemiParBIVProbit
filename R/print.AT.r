print.AT <- function(x, ...){

if(x$mar2=="probit") es <- format(x$res*100, digits = 3, trim=TRUE)

if(x$mar2!="probit") es <- format(x$res, digits = 3, trim=TRUE)

if(x$mar2=="probit") cat("\nAverage treatment effect (%) with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")

if(x$mar2!="probit") cat("\nAverage treatment effect with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")

cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

invisible(x)

}
