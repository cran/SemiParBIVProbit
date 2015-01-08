sem.checks <- function(x){

e.v <- eigen(x$fit$hessian, symmetric=TRUE, only.values = TRUE)$values

n.drop <- x$n - sum(as.numeric(x$good==TRUE))

cat("\nLargest absolute gradient value:",max(abs(x$fit$gradient)))

if(x$hess==TRUE) mv <- "Observed" else mv <- "Expected" 

if (min(e.v) > 0) cat("\n",mv," information matrix is positive definite\n",sep="") else cat("\n",mv," information matrix is NOT positive definite\n",sep="")

cat("Eigenvalue range: [",min(e.v),",",max(e.v),"]\n", sep = "")

if( ((x$l.sp1!=0 || x$l.sp2!=0) && x$fp==FALSE) || (x$PL!="P" && x$fitPL=="pLiksp" ) ){

cat("\nTrust region iterations before smoothing parameter estimation:",x$iter.if)
cat("\nLoops for smoothing parameter estimation:",x$iter.sp) 
cat("\nTrust region iterations within smoothing loops:",x$iter.inner,"\n\n")

}else{cat("\nTrust region iterations:",x$iter.if,"\n\n")}

if(!is.null(x$conv.sp)){
if(x$conv.sp == FALSE) cat("WARNING: Smoothing algorithm reached the max. number of iterations allowed.\n\n")
if(x$conv.sp == FALSE) eb <- "\n" else eb <-""
}

if(n.drop > 0) cat(eb,"WARNING: Number of observations dropped during fitting:",n.drop,"\n\n", sep = "")


}
