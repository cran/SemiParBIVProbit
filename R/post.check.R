post.check <- function(x, main = "Histogram and Density Estimate of Residuals", 
                          main2 = "Histogram and Density Estimate of Residuals",
                           xlab = "Quantile Residuals", xlab2 = "Quantile Residuals", ...){

if(x$Cont == "NO"){

p <- distrHsAT(x$y2, x$eta2, x$sigma2, x$nu, x$margins[2])$p2 
qqp <- qnorm(p)


par(mfrow = c(1, 2))

hist(qqp, freq = FALSE, #ylim=c(0, max(qqp, hist(qqp, plot = FALSE)$density) ),
     main=main,
     xlab=xlab, ylab = "Density", ...)
lines(density(qqp, adjust = 2),lwd=2)

qqnorm(qqp)
abline(0, 1, col = "red")


}




if(x$Cont == "YES"){

p1 <- distrHsAT(x$y1, x$eta1, x$sigma21, x$nu1, x$margins[1])$p2 
p2 <- distrHsAT(x$y2, x$eta2, x$sigma22, x$nu2, x$margins[2])$p2 

par(mfrow = c(2, 2))

qqp <- qnorm(p1)
hist(qqp, freq = FALSE, 
     main=main,
     xlab=xlab, ylab = "Density", ...)
lines(density(qqp, adjust = 2),lwd=2)

qqnorm(qqp)
abline(0, 1, col = "red")

qqp <- qnorm(p2)
hist(qqp, freq = FALSE, 
     main=main2,
     xlab=xlab2, ylab = "Density", ...)
lines(density(qqp, adjust = 2),lwd=2)

qqnorm(qqp)
abline(0, 1, col = "red")

}




}

