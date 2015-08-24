resp.check <- function(y, margin = "N", 
                           main = "Histogram and Density of Response",
                           xlab="Response", print.par = FALSE, plots = TRUE, ...){

m2 <- c("N","GU","rGU","LO","LN","WEI","iG","GA")
m3 <- c("DAGUM")

if(!(margin %in% c(m2,m3)) ) stop("Error in margin value. It can be: N, GU, rGU, LO, LN, WEI, iG, GA, DAGUM.") 

if(margin %in% c("LN","WEI","WEI2","iG","GA","DAGUM") && min(y, na.rm = TRUE) <= 0) stop("The response must be positive.")

y <- na.omit(y)

if(margin == "LN") y <- log(y)

margins <- c("probit", margin)

VC <- list(X2 = matrix(1, nrow = length(y), ncol = 1), X2.d2 = 1,
           X3 = matrix(1, nrow = length(y), ncol = 1), X3.d2 = 1,
           X4 = matrix(1, nrow = length(y), ncol = 1), X4.d2 = 1,
           l.sp2 = 0, l.sp3 = 0, l.sp4 = 0, 
           weights = 1, 
           margins = margins, fp = TRUE,
           extra.regI = "t")

respvec <- list(y2 = y)
           

if( margin %in% c("N","LN") )   st.v <- c( mean((y + mean(y))/2) ,           log( var(y) ) )  
if( margin %in% c("LO") )       st.v <- c( mean((y + mean(y))/2) ,           log(  3*var(y)/pi^2 ) )  
#if( margin %in% c("LN") )       st.v <- c( mean((log(y) + mean(log(y)))/2) , log(var(log(y))) )           # not sure log in front of mean   
if( margin %in% c("iG") )       st.v <- c( log( mean((y + mean(y))/2) ) , log( var(y)/mean(y)^3)  )    
if( margin %in% c("GU") )       st.v <- c( mean(y) + 0.57722*sqrt(var(y)/1.64493) ,  log(6*var(y)/pi^2) )    
if( margin %in% c("rGU") )      st.v <- c( mean(y) - 0.57722*sqrt(var(y)/1.64493) ,  log(6*var(y)/pi^2) )   
if( margin %in% c("WEI") )      st.v <- c( log( mean( exp(log(y) + 0.5772/(1.283/sqrt(var(log(y))))) )  ) , log( ( 1.283/sqrt(var(log(y))) )^2 ) ) 
if( margin %in% c("GA") )       st.v <- c( log(mean((y + mean(y))/2)), log(var(y)/mean(y)^2)  ) # log( 1^2 )             
if( margin %in% c("DAGUM") )    st.v <- c( log(mean((y + mean(y))/2)), log(sqrt(2)), log(1) )

if( margin %in% m2 ) names(st.v) <- c("eta2", "sigma2.star")
if( margin %in% m3 ) names(st.v) <- c("eta2", "sigma2.star","nu.star")


if(margin != "DAGUM") univfit <-  try(trust(bprobgHsContUniv, st.v, rinit = 1, rmax = 100, respvec = respvec, 
                                            VC = VC, sp = NULL, qu.mag = NULL, blather = TRUE), silent = TRUE)
                                            
if(margin == "DAGUM") univfit <-  try(trust(bprobgHsContUniv3, st.v, rinit = 1, rmax = 100, respvec = respvec, 
                                            VC = VC, sp = NULL, qu.mag = NULL, blather = TRUE), silent = TRUE)                                            
                 
                 
             
if(plots == TRUE){             
             
if(class(univfit) == "try-error") stop("The parameters of the chosen distribution could not be estimated. Try a different distribution.")                  

# univfit$argument

if(margin == "LN") y <- exp(y)

if(margin != "DAGUM") pp <- distrHsAT(y, univfit$argument[1], exp(univfit$argument[2]), 1, margin)
if(margin == "DAGUM") pp <- distrHsAT(y, univfit$argument[1], exp(univfit$argument[2]), exp(univfit$argument[3]), margin)


p <- pp$p2
d <- pp$pdf2

par(mfrow = c(1, 2))
hist(y, freq = FALSE, ylim=c(0, max(d, hist(y, plot = FALSE)$density) ),
     main=main,
     xlab=xlab, ...)

#xspline(sort(y),d[order(y)],lwd=2)

lines(sort(y),d[order(y)],lwd=2)

#lines(density(d[order(y)], adjust=1),lwd=2)
#lines(density(d, adjust=1),lwd=2)





qqnorm(qnorm(p))
abline(0, 1, col = "red")

if(print.par == TRUE) print(univfit$argument)


}

if(plots == FALSE) return( univfit$argument ) 



}

