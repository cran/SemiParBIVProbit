print.AT <- function(x,...){

es <- format(x$res, digits = 3)

cat("\nAverage Effect with corresponding ",(1-x$sig.lev)*100,"% 'Confidence' Intervals:\n\n",sep="")
cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

invisible(x)

}
