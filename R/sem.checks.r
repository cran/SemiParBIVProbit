sem.checks <- function(x){

e.v <- eigen(x$He, symmetric=TRUE, only.values = TRUE)$values

n.drop <- x$n-sum(as.numeric(x$good==TRUE))

cat("\nLargest absolute gradient value:",max(abs(x$fit$gradient)))

if(x$hess==TRUE) mv <- "Observed" else mv <- "Expected" 

if (min(e.v) > 0) cat("\n",mv," information matrix is positive definite\n",sep="")
else cat("\nEigenvalue range of the information matrix: [",min(e.v),",",max(e.v),"]\n", sep = "")

if( ((x$l.sp1!=0 || x$l.sp2!=0) && x$fp==FALSE) || (x$PL!="P" && x$fitPL=="pLiksp" ) ){

cat("\nTrust region iterations before smoothing parameter estimation:",x$iter.if)
cat("\nLoops for smoothing parameter estimation:",x$iter.sp) #,"\n\n")
cat("\nTrust region iterations within smoothing loops:",x$iter.inner,"\n\n")
}
                        else{cat("\nTrust region iterations:",x$iter.if,"\n\n")}

if(n.drop > 0) cat("Number of observations dropped during fitting:",n.drop,"\n\n")


}
