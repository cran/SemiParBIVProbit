mb <- function(treat, outc){

if(length(table(treat))!=2 ) stop("The treatment variable must be binary.")
if(length(table(outc))!=2  ) stop("The outcome variable must be binary.")

tab1 <- table(treat)
pT1  <- prop.table(tab1)[2] 
pT0  <- 1 - pT1 
tab2 <- table(treat, outc)

pY1cT1 <- prop.table(tab2,1)[4] 
pY1cT0 <- prop.table(tab2,1)[3] 

stuff  <- pY1cT1*pT1 - pY1cT0*pT0

UB <- stuff + pT0
LB <- stuff - pT1

avte <- (pY1cT1 - pY1cT0)*100
avte <- format(avte, digits = 3, trim=TRUE)

WCB <- c(LB, UB)*100 

WCB <- format(WCB, digits = 3, trim=TRUE)

cat("\nWORST-CASE MANSKI'S BOUND (%)\n\n","Lower Bound: ",WCB[1],"\n","Upper Bound: ",WCB[2],"\n\n",sep="")

cat("\nATE Assuming Random Assignment (%): ",avte,"\n\n",sep="")

}
