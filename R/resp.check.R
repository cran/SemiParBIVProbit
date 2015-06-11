resp.check <- function(y, margin = "N", 
                           main = "Histogram and Density of Response",
                           xlab="Response", ...){

m2 <- c("N","GU","rGU","LO","LN","WEI","iG","GA")

if(!(margin %in% m2) ) stop("Error in margin value. It can be: N, GU, rGU, LO, LN, WEI, iG, GA.") 

if(margin %in% c("LN","WEI","iG","GA") && min(y) < 0) stop("The variable of interest must be positive.")

y <- na.omit(y)

if(margin == "LN") y <- log(y)

margins <- c("probit", margin)

VC <- list(X2 = matrix(1, nrow = length(y), ncol = 1), X2.d2 = 1,
           X3 = matrix(1, nrow = length(y), ncol = 1), X3.d2 = 1,
           l.sp2 = 0, l.sp3 = 0, 
           weights = 1, 
           margins = margins, fp = TRUE,
           extra.regI = "t")

respvec <- list(y2 = y)
           

univfit <- trust(bprobgHsContUniv, c(1,1), rinit = 1, rmax = 100, respvec = respvec, 
                 VC = VC, sp = NULL, qu.mag = NULL)

# univfit$argument

if(margin == "LN") y <- exp(y)

pp <- distrHsAT(y, univfit$argument[1], exp(univfit$argument[2]), margin)

p <- pp$p2
d <- pp$pdf2

par(mfrow = c(1, 2))
hist(y, freq = FALSE, ylim=c(0,max(d,hist(y, plot = FALSE)$density)),
     main=main,
     xlab=xlab, ...)
xspline(sort(y),d[order(y)],lwd=2)

qqnorm(qnorm(p))
abline(0, 1, col = "red")

}

